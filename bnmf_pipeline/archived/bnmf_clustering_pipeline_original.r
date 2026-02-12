pacman::p_load(tidyverse, data.table, readxl, magrittr, dplyr, strex,
               DT, kableExtra, GenomicRanges,rstudioapi,here)




## need to load two R files dynamically 
source("/mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/bnmf-pipeline/bnmf_pipeline/R/run_bNMF_2025.r")
source("/mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/bnmf-pipeline/bnmf_pipeline/R/utilities.R")
source("/mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/bnmf-pipeline/bnmf_pipeline/R/utilities.R")
rmd_path='/mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/bnmf-pipeline/bnmf_pipeline/R/format_bNMF_results.Rmd'
# USER INPUTS!!!
project_dir = '/mnt/disks/sdd/bnmf-clustering/Neuroimage/SCZ/Neuroimage_main_cluster/test_run/' # path to where you want results saved
user_token = '901021003036' # token for LDlinkR api
PVCUTOFF = 5e-8
setwd(project_dir)

main_gwas_id='daner_PGC_SCZ_w3_90_0418b_ukbbdedupe'
z_score_file='daner_PGC_SCZ_w3_90_0418b_ukbbdedupe_zscore_index.csv'
sample_size_file='daner_PGC_SCZ_w3_90_0418b_ukbbdedupe_sample_size_index.csv'
#----




# create project folder 
dir.create(project_dir)

z_mat<-fread(z_score_file)
N_mat<-fread(sample_size_file)

## Please need to check both files column names and row names are identical or not 


z_mat <- as.data.frame(z_mat)
z_mat$uniq_id<-NULL
rownames(z_mat) <-z_mat$chr_pos
z_mat$chr_pos<-NULL


N_mat <- as.data.frame(N_mat)
N_mat$uniq_id<-NULL
rownames(N_mat) <- N_mat$chr_pos
N_mat$chr_pos<-NULL



colSums(is.na(z_mat))
colSums(is.na(N_mat))








prep_z_output <- prep_z_matrix(z_mat = z_mat,
                               N_mat = N_mat,
                               corr_cutoff = 0.8)

# prep_z_output has two outputs:
#   1.) The scaled, non-negative z-score matrix
final_zscore_matrix <- prep_z_output$final_z_mat

#   2.) Results from the trait filtering
df_traits_filtered <- prep_z_output$df_traits
write_csv(x = df_traits_filtered,
          file = file.path(project_dir,"df_traits.csv"))

# prep_z_matrix also save trait correlation matrix to working dir, so move to project dir
system(sprintf("mv trait_cor_mat.txt %s", project_dir))

print(sprintf("Final matrix: %i SNPs x %i traits",
      nrow(final_zscore_matrix),
      ncol(final_zscore_matrix)/2))

save.image(file = file.path(project_dir, "pipeline_data.RData"))

#----

final_zscore_matrix <- as.matrix(final_zscore_matrix)
storage.mode(final_zscore_matrix) <- "double"


# Section 11.) Run bNMF 
bnmf_reps <- run_bNMF(final_zscore_matrix,
                      n_reps=100,
                      tolerance = 1e-6)

final_active <- sapply(bnmf_reps, function(x) {
  as.numeric(tail(unlist(x$n.active), 1))
  })





summarize_bNMF(bnmf_reps, dir_save=project_dir)

save.image(file = file.path(project_dir, "pipeline_data.RData"))

end=Sys.time()
print("Total pipeline runtime:")
print(end-start)

#----



  # Calculate stability summary
  final_active <- sapply(bnmf_reps, function(x) as.numeric(tail(unlist(x$n.active), 1)))
  k_summary <- as.data.frame(table(final_active)) %>%
    dplyr::rename(K = final_active, n_runs = Freq) %>%
    dplyr::mutate(K = as.numeric(as.character(K)))
  
  top3_k_vector <- k_summary %>% 
    dplyr::arrange(dplyr::desc(n_runs)) %>% 
    dplyr::slice_head(n = 3) %>% 
    dplyr::pull(K)



# format results


for (k in top3_k_vector ) {

  if (is.null(k)){
    html_filename <- "results_for_maxK.html"
  } else {
    html_filename <- sprintf("results_for_K_%i.html", k)
  }

  out_dir <- glue::glue("{project_dir}/reports")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


  rmarkdown::render(
    input       = rmd_path,
    output_file = html_filename,
    output_dir  = out_dir,
    params = list(
      main_dir  = project_dir,
      k         = k,
      loci_file = "query",
      GTEx      = FALSE,
      my_traits = data.frame()
    )
  )



}


