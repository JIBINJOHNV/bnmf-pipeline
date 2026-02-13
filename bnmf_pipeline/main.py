import os
import sys
from pathlib import Path
import rpy2.robjects as robjects
import subprocess

# Relative path to your R scripts inside the package
INTERNAL_R_DIR = Path(__file__).parent / "R"

def run_full_pipeline(configs):
    """
    Orchestrates GWAS processing and bNMF clustering.
    """
    # 1. Output Directory Logic
    output_folder = configs.get('output_folder')
    final_output_path = os.path.abspath(output_folder) if output_folder else os.getcwd()
    os.makedirs(final_output_path, exist_ok=True)
    
    main_id = configs['main_gwas_id']

    # 3. Step 1: Main GWAS Summary Stat Processing
    print(f"--- [PYTHON] Starting Main GWAS Processing ---")
    from .process_main_sumstat import process_main_gwas
    
    # This returns a dictionary: {"main_sumstat_tsv": "/path/to/file.tsv"}
    step1_output = process_main_gwas(
        main_gwas_vcf_path=configs['main_gwas_vcf_path'],
        outdir=final_output_path,
        main_gwas_id=main_id,
        ld_folder=configs['ld_folder'],
        bcftools_bin=configs['bcftools_bin'],
        n_threads=configs['n_threads'],
        pop=configs['pop'],
        lead_p=configs['lead_p'],
        r2_clump=configs['r2_clump'],
        snp_to_include=configs['snp_to_include']
    )

    # EXTRACT THE PATH STRING FROM THE DICTIONARY
    # This prevents the Polars .read() error in Step 2
    main_gwas_loci = step1_output["main_sumstat_tsv"]

    # 2. Step 1: Trait Preparation (Polars)
    print(f"--- [PYTHON] Starting Trait Processing in {final_output_path} ---")
    from .process_traits_sumstat import process_traits_gwas_polars
    
    trait_outpit = process_traits_gwas_polars(
        main_gwas_id=main_id,
        main_gwas_loci=main_gwas_loci, # Now a proper string path
        ld_folder=configs['ld_folder'],
        trait_input_file=configs['trait_input_file'],
        dbsnp_folder=configs['dbsnp_folder'],
        output_folder=final_output_path,
        r2_threshold=configs['r2_threshold'],
        window=configs['window'],
        bcftools_bin=configs['bcftools_bin'],
        variant_missing_rate=configs['variant_missing_rate'],
        sample_missing_rate=configs['sample_missing_rate']
    )

    # 4. Step 3: bNMF Clustering (R via rpy2)
    print(f"--- [R via rpy2] Initializing bNMF Clustering ---")
    
    # FIX: Ensure we are capturing string paths. 
    # Based on your process_traits_gwas_polars return statement, 
    # we use 'zscore_file' and 'ss_final' keys.
    z_file = trait_outpit['zscore_file']
    n_file = trait_outpit['ss_final']
    
    # Validation: Ensure these are strings before calling os.path.exists
    if not isinstance(z_file, str) or not isinstance(n_file, str):
        # This catch is for safety if the R-script return wasn't updated yet
        raise TypeError(f"Expected file paths as strings, but got {type(z_file)} and {type(n_file)}")

    if not (os.path.exists(z_file) and os.path.exists(n_file)):
        raise FileNotFoundError(f"Input CSVs missing in {final_output_path}. Ensure Step 2 saved files correctly.")


    print(f"--- [Rscript] Launching bNMF Clustering CLI ---")

    r_cli_script = str(INTERNAL_R_DIR / "run_bnmf_cli.R")

    cmd = [
        "Rscript",
        r_cli_script,
        "--project_dir", final_output_path,
        "--z_score_file", z_file,
        "--sample_size_file", n_file,
        "--main_gwas_id", main_id,
        "--n_reps", str(configs["n_reps"]),
        "--tolerance", str(configs["tolerance"]),
        "--corr_cutoff", str(configs["corr_cutoff"]),
        "--maximum_k", str(configs["maximum_k"]),
        "--script_path", str(INTERNAL_R_DIR)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError("bNMF R pipeline failed")

    print(result.stdout)
    print(f"Success! Results saved in: {final_output_path}")


    # Source the specific R pipeline script
    # r_pipeline = str(INTERNAL_R_DIR / "bnmf_clustering_pipeline.r")
    # robjects.r.source(r_pipeline)
    # bnmf_func = robjects.globalenv['run_genomic_bnmf_pipeline']
    
    # bnmf_func(
    #     project_dir      = final_output_path,
    #     z_score_file     = z_file,
    #     sample_size_file = n_file,
    #     main_gwas_id     = main_id,
    #     n_reps           = configs['n_reps'],
    #     tolerance        = configs['tolerance'],
    #     corr_cutoff      = configs['corr_cutoff'],
    #     maximum_k        = configs['maximum_k'],
    #     script_path      = str(INTERNAL_R_DIR)
    # )
    # print(f"Success! Results saved in: {final_output_path}")
    
