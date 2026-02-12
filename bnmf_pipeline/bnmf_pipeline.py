import argparse
import os
import sys
from pathlib import Path
import rpy2.robjects as robjects

# Assuming these are in your module:
# from your_module import process_traits_gwas_polars, process_main_gwas

def main():
    parser = argparse.ArgumentParser(
        description="Integrated GWAS & bNMF Clustering Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --- Global / Shared Arguments ---
    parser.add_argument("--main_gwas_id", default="daner_PGC_SCZ_w3_90_0418b_ukbbdedupe", help="Primary GWAS identifier")
    parser.add_argument("--ld_folder", default="/mnt/disks/sdd/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/", help="Path to LD reference panel")
    parser.add_argument("--bcftools_bin", default="/usr/bin/bcftools", help="Path to bcftools executable")
    parser.add_argument("--n_threads", type=int, default=10, help="Number of threads for processing")
    parser.add_argument("--output_folder", default=None, help="Output path (Defaults to Current Working Directory)")

    # --- Trait Processing Arguments ---
    trait_group = parser.add_argument_group("Trait Processing Parameters")
    trait_group.add_argument("--main_gwas_loci", default="/mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe_independent_markers.tsv")
    trait_group.add_argument("--trait_input_file", default="/mnt/disks/sdd/bnmf-clustering/filtered_vcf_files/trait_gwas.csv")
    trait_group.add_argument("--dbsnp_folder", default="/mnt/disks/sdd/resourses/postgwas/gwas2vcf/GRCh37/dbSNP/vcf_files/")
    trait_group.add_argument("--r2_threshold", type=float, default=0.8)
    trait_group.add_argument("--window", type=int, default=1000000)

    # --- Main GWAS Processing Arguments ---
    gwas_group = parser.add_argument_group("Main GWAS Processing Parameters")
    gwas_group.add_argument("--vcf_path", default="/mnt/disks/sdd/bnmf-clustering/filtered_vcf_files/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe_GRCh37_filtered.vcf.gz")
    gwas_group.add_argument("--snp_to_include", default="/mnt/disks/sdd/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/EUR.chr1_22_Variants_to_include.txt")
    gwas_group.add_argument("--pop", default="EUR")
    gwas_group.add_argument("--lead_p", type=float, default=5e-8)
    gwas_group.add_argument("--r2_clump", type=float, default=0.05)

    # --- bNMF Clustering Arguments ---
    bnmf_group = parser.add_argument_group("bNMF Clustering Parameters")
    bnmf_group.add_argument("--n_reps", type=int, default=100, help="Number of bNMF repetitions")
    bnmf_group.add_argument("--tolerance", type=float, default=1e-6)
    bnmf_group.add_argument("--corr_cutoff", type=float, default=0.8)

    # --- Execution Flags ---
    parser.add_argument("--run_traits", action="store_true", help="Step 1: Process Traits (Polars)")
    parser.add_argument("--run_main", action="store_true", help="Step 2: Process Main GWAS")
    parser.add_argument("--run_bnmf", action="store_true", help="Step 3: Run bNMF Clustering (R via rpy2)")

    args = parser.parse_args()

    # --- Path Resolution ---
    script_dir = Path(__file__).parent.resolve()
    
    # Logic: Default to CWD if no output folder is provided
    if args.output_folder is None:
        final_output_path = os.getcwd()
    else:
        final_output_path = os.path.abspath(args.output_folder)

    os.makedirs(final_output_path, exist_ok=True)

    # --- STEP 1: Process Traits ---
    if args.run_traits:
        print(f"--- [PYTHON] Starting Trait Processing in {final_output_path} ---")
        process_traits_gwas_polars(
            main_gwas_id=args.main_gwas_id,
            main_gwas_loci=args.main_gwas_loci,
            ld_folder=args.ld_folder,
            trait_input_file=args.trait_input_file,
            dbsnp_folder=args.dbsnp_folder,
            output_folder=final_output_path,
            r2_threshold=args.r2_threshold,
            window=args.window,
            bcftools_bin=args.bcftools_bin
        )

    # --- STEP 2: Process Main GWAS ---
    if args.run_main:
        print(f"--- [PYTHON] Starting Main GWAS Processing ---")
        process_main_gwas(
            main_gwas_vcf_path=args.vcf_path,
            outdir=final_output_path,
            main_gwas_id=args.main_gwas_id,
            ld_folder=args.ld_folder,
            bcftools_bin=args.bcftools_bin,
            n_threads=args.n_threads,
            pop=args.pop,
            lead_p=args.lead_p,
            r2_clump=args.r2_clump,
            snp_to_include=args.snp_to_include
        )

    # --- STEP 3: Run bNMF Clustering (via rpy2) ---
    if args.run_bnmf:
        print(f"--- [R via rpy2] Initializing bNMF Clustering ---")
        
        # Define expected input files from previous steps
        z_file = os.path.join(final_output_path, f"{args.main_gwas_id}_zscore_index.csv")
        n_file = os.path.join(final_output_path, f"{args.main_gwas_id}_sample_size_index.csv")

        # Sanity Check
        if not (os.path.exists(z_file) and os.path.exists(n_file)):
            sys.exit(f"Error: Required input files not found in {final_output_path}.\nEnsure Step 1 & 2 completed successfully.")

        # Load and call the R function
        r_source_path = str(script_dir / "run_genomic_bnmf_pipeline.R")
        robjects.r.source(r_source_path)
        bnmf_pipeline_r = robjects.globalenv['run_genomic_bnmf_pipeline']
        
        results = bnmf_pipeline_r(
            project_dir      = final_output_path,
            z_score_file     = z_file,
            sample_size_file = n_file,
            main_gwas_id     = args.main_gwas_id,
            n_reps           = args.n_reps,
            tolerance        = args.tolerance,
            corr_cutoff      = args.corr_cutoff,
            script_path      = str(script_dir)
        )
        print(f"Pipeline Complete. Results saved in: {final_output_path}")

    if not any([args.run_traits, args.run_main, args.run_bnmf]):
        print("Done. No execution flags provided. Use --help for usage details.")

if __name__ == "__main__":
    main()