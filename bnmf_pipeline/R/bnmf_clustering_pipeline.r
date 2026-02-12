#' Production-grade bNMF Pipeline
#' @description Portable function that auto-detects sibling scripts and syncs data.
#' Production-grade bNMF Pipeline
#' @param script_path Explicit path to the directory containing run_bNMF_2025.R and utilities.R.
#' If NULL, the function attempts to auto-detect its own location.
#' Production-grade bNMF Pipeline
#' @description Portable function for genomic bNMF clustering. 
#' Handles Python integration and prevents package namespace conflicts.
run_genomic_bnmf_pipeline <- function(
    project_dir, 
    z_score_file, 
    sample_size_file, 
    main_gwas_id,
    n_reps = 100, 
    tolerance = 1e-6,
    corr_cutoff = 0.8,
    script_path = NULL
) {
  start_time <- Sys.time()
  
  # 1. Load Dependencies
  if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
  # Note: GenomicRanges often masks dplyr functions; we address this with explicit calls below
  pacman::p_load(tidyverse, data.table, readxl, magrittr, dplyr, strex, 
                 DT, kableExtra, GenomicRanges, rstudioapi, this.path, glue)

  # 2. Path Resolution (Python-style: Path(__file__).parent)
  if (is.null(script_path)) {
    script_dir <- tryCatch({
      this.path::this.dir()
    }, error = function(e) {
      if (rstudioapi::isAvailable()) {
        dirname(rstudioapi::getActiveDocumentContext()$path)
      } else {
        message("Warning: script_path not provided and auto-detection failed. Using getwd().")
        getwd()
      }
    })
  } else {
    script_dir <- script_path
  }

  # 3. Source Sibling Files
  source(file.path(script_dir, "run_bNMF_2025.r"))
  source(file.path(script_dir, "utilities.R"))
  
  # 4. Workspace Setup
  if (!dir.exists(project_dir)) dir.create(project_dir, recursive = TRUE)
  project_dir <- normalizePath(project_dir)
  setwd(project_dir)

  # 5. Data Ingestion & Alignment
  message("--- Loading and Syncing Matrices ---")
  z_raw <- data.table::fread(z_score_file)
  n_raw <- data.table::fread(sample_size_file)

  common_snps <- intersect(z_raw$chr_pos, n_raw$chr_pos)
  if(length(common_snps) == 0) stop("Zero overlapping SNPs between Z and N matrices.")

  # Helper: Explicitly using dplyr namespace to avoid GenomicRanges/S4Vectors conflicts
  process_df <- function(df, snps) {
    df %>%
      dplyr::filter(chr_pos %in% snps) %>%
      dplyr::arrange(chr_pos) %>%
      as.data.frame() %>%
      tibble::column_to_rownames("chr_pos") %>%
      dplyr::select(-dplyr::any_of("uniq_id"))
  }

  z_mat <- process_df(z_raw, common_snps)
  n_mat <- process_df(n_raw, common_snps)

  # 6. Matrix Preparation
  message("--- Preparing Z-matrix (Handling Signs) ---")
  prep_z_output <- prep_z_matrix(z_mat = z_mat, N_mat = n_mat, corr_cutoff = corr_cutoff)
  
  final_zscore_matrix <- as.matrix(prep_z_output$final_z_mat)
  storage.mode(final_zscore_matrix) <- "double"

  readr::write_csv(prep_z_output$df_traits, file.path(project_dir, "df_traits.csv"))
  if (file.exists("trait_cor_mat.txt")) {
    file.rename("trait_cor_mat.txt", file.path(project_dir, "trait_cor_mat.txt"))
  }

  # 7. Execute bNMF
  message(sprintf("--- Running %i bNMF Iterations ---", n_reps))
  bnmf_reps <- run_bNMF(final_zscore_matrix, n_reps = n_reps, tolerance = tolerance)
  
  # 8. Stability Analysis & Results Saving
  summarize_bNMF(bnmf_reps, dir_save = project_dir)
  save_bnmf_results(bnmf_reps, main_gwas_id, project_dir)

  # Calculate stability summary
  final_active <- sapply(bnmf_reps, function(x) as.numeric(tail(unlist(x$n.active), 1)))
  k_summary <- as.data.frame(table(final_active)) %>%
    dplyr::rename(K = final_active, n_runs = Freq) %>%
    dplyr::mutate(K = as.numeric(as.character(K)))
  
  top3_ks <- k_summary %>% 
    dplyr::arrange(dplyr::desc(n_runs)) %>% 
    dplyr::slice_head(n = 3) %>% 
    dplyr::pull(K)

  # 9. Report Generation
  rmd_path <- file.path(script_dir, "format_bNMF_results.Rmd") 
  out_dir <- file.path(project_dir, "reports")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  for (k in top3_ks) {
    message(sprintf("Rendering Report for K=%i", k))
    tryCatch({
      rmarkdown::render(
        input       = rmd_path,
        output_file = sprintf("results_for_K_%i.html", k),
        output_dir  = out_dir,
        params      = list(main_dir = project_dir, k = k, loci_file = "query", GTEx = FALSE, my_traits = data.frame()),
        quiet       = TRUE
      )
    }, error = function(e) message(sprintf("Failed to render K=%i: %s", k, e$message)))
  }

  # 10. Final Export
  save.image(file = file.path(project_dir, "pipeline_complete.RData"))
  message("--- Pipeline Finished ---")
  print(Sys.time() - start_time)
  
  return(list(summary = k_summary, best_ks = top3_ks))
}




# ## need to load two R files dynamically 
# source("/Users/JJOHN41/Downloads/bnmf-clustering/scripts/run_bNMF_2025.R")
# source("/Users/JJOHN41/Downloads/bnmf-clustering/scripts/utilities.R")

# # USER INPUTS!!!
# project_dir = '/mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe/' # path to where you want results saved
# user_token = '901021003036' # token for LDlinkR api
# PVCUTOFF = 5e-8
# setwd(project_dir)

# main_gwas_id='daner_PGC_SCZ_w3_90_0418b_ukbbdedupe'
# z_score_file='daner_PGC_SCZ_w3_90_0418b_ukbbdedupe_zscore_index.csv'
# sample_size_file='daner_PGC_SCZ_w3_90_0418b_ukbbdedupe_sample_size_index.csv'
# script_path='/mnt/disks/sdd/bnmf-clustering/bnmf_R_scripts/'



# run_genomic_bnmf_pipeline(
#     project_dir=project_dir, 
#     z_score_file=z_score_file, 
#     sample_size_file=sample_size_file, 
#     main_gwas_id=main_gwas_id,
#     n_reps = 100, 
#     tolerance = 1e-6,
#     corr_cutoff = 0.8,
#     script_path = script_path
# )