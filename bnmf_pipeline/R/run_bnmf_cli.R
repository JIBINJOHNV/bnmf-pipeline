#!/usr/bin/env Rscript

# Load library quietly to keep output clean
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-d", "--project_dir"), type="character", default=NULL, 
              help="Path to the project directory", metavar="path"),
  make_option(c("-z", "--z_score_file"), type="character", default=NULL, 
              help="Path to the synchronized Z-score matrix", metavar="path"),
  make_option(c("-s", "--sample_size_file"), type="character", default=NULL, 
              help="Path to the effective sample size matrix", metavar="path"),
  make_option(c("-i", "--main_gwas_id"), type="character", default="GWAS_Project", 
              help="Identifier for the main GWAS study [default %default]", metavar="string"),
  make_option(c("-r", "--n_reps"), type="integer", default=100, 
              help="Number of bNMF repetitions [default %default]", metavar="int"),
  make_option(c("-t", "--tolerance"), type="double", default=1e-06, 
              help="Convergence tolerance [default %default]", metavar="num"),
  make_option(c("-c", "--corr_cutoff"), type="double", default=0.8, 
              help="Correlation cutoff for redundant traits [default %default]", metavar="num"),
  make_option(c("-p", "--script_path"), type="character", default=NULL, 
              help="Path to the core bNMF R source scripts", metavar="path"),
  make_option(c("-k", "--maximum_k"), type = "integer", default = 30,
              help = "Maximum number of latent factors (K) allowed. Default: %default.")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "\nBNMF Genomic Clustering Step 3: R Core Pipeline"
)

opt <- parse_args(opt_parser)

# --- Validation Logic ---
required_args <- c("project_dir", "z_score_file", "sample_size_file", "script_path")
missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]))]

if (length(missing_args) > 0) {
  print_help(opt_parser)
  cat("\n[!] ERROR: The following required arguments are missing:\n")
  cat(paste("  -", missing_args, collapse = "\n"), "\n\n")
  quit(status = 1)
}

# --- Execute ---
core_script <- file.path(opt$script_path, "bnmf_clustering_pipeline.r")

if (!file.exists(core_script)) {
  stop(paste("Could not find core script at:", core_script), call. = FALSE)
}

# Source the function definition
source(core_script)

message(paste(">>> Starting bNMF pipeline for:", opt$main_gwas_id))

# --- ARGUMENT FILTERING (The Fix) ---
# 1. Get names of arguments that the function actually expects
valid_args <- names(formals(run_genomic_bnmf_pipeline))

# 2. Filter 'opt' to remove 'help' and other CLI-only artifacts
# This ensures only recognized parameters are passed
pipeline_params <- opt[names(opt) %in% valid_args]

# 3. Call the function safely
tryCatch({
    do.call(run_genomic_bnmf_pipeline, pipeline_params)
    message(">>> bNMF pipeline completed successfully.")
}, error = function(e) {
    message(paste("\n[!] Pipeline Error:", e$message))
    quit(status = 1)
})