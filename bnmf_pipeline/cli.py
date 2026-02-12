import argparse
import sys
import textwrap
from rich_argparse import RichHelpFormatter
from .main import run_full_pipeline

def main():
    # Setup RichHelpFormatter styling for a professional look
    RichHelpFormatter.styles["argparse.groups"] = "bold gold1"
    RichHelpFormatter.styles["argparse.args"] = "cyan"
    RichHelpFormatter.styles["argparse.metavar"] = "italic underline grey70"
    RichHelpFormatter.styles["argparse.default"] = "bold dim white"  # Highlighted defaults

    parser = argparse.ArgumentParser(
        prog="bnmf-pipeline",
        description=textwrap.dedent('''
            [bold cyan]BNMF Genomic Clustering Pipeline[/bold cyan]
            [hr]
            This tool orchestrates the processing of GWAS summary statistics and 
            performs Bayesian Non-negative Matrix Factorization (bNMF) clustering.
            
            [italic]Example Usage:[/italic]
              [green]bnmf-pipeline --run_traits --run_main --run_bnmf[/green]
              [green]bnmf-pipeline --run_bnmf --n_reps 250 --main_gwas_id MY_GWAS_001[/green]
        '''),
        formatter_class=RichHelpFormatter
    )

    # --- Global Configuration ---
    global_group = parser.add_argument_group("Global Configuration")

    global_group.add_argument("--ld_folder", default="/mnt/disks/sdd/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/", metavar="",
                             help="Path to 1000G LD reference panels. (Default: %(default)s)")
    global_group.add_argument("--bcftools_bin", default="/usr/bin/bcftools", metavar="",
                             help="Path to bcftools executable. (Default: %(default)s)")
    global_group.add_argument("--n_threads", type=int, default=10, metavar="",
                             help="Number of CPU threads. (Default: %(default)s)")
    global_group.add_argument("--output_folder", default=None, metavar="",
                             help="Output directory. (Default: Current Working Directory)")

    # --- Step 2: GWAS Summary Stat Processing ---
    gwas_group = parser.add_argument_group("Step 2: GWAS Summary Stat Processing")
    gwas_group.add_argument("--main_gwas_id", default="daner_PGC_SCZ_w3_90_0418b_ukbbdedupe", metavar="",
                             help="Unique identifier for the primary GWAS study. (Default: %(default)s)")
    gwas_group.add_argument("--main_gwas_vcf_path", 
                             default="/mnt/disks/sdd/bnmf-clustering/filtered_vcf_files/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe_GRCh37_filtered.vcf.gz",
                             metavar="", help="Path to filtered main GWAS VCF. (Default: %(default)s)")
    gwas_group.add_argument("--snp_to_include", 
                             default="/mnt/disks/sdd/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/EUR.chr1_22_Variants_to_include.txt",
                             metavar="", help="List of high-quality SNPs. (Default: %(default)s)")
    gwas_group.add_argument("--pop", default="EUR", 
                             metavar="", help="Population code (e.g., EUR). (Default: %(default)s)")
    gwas_group.add_argument("--lead_p", type=float, default=5e-8, 
                             metavar="", help="P-value threshold for lead SNPs. (Default: %(default)s)")
    gwas_group.add_argument("--r2_clump", type=float, default=0.05, 
                             metavar="", help="R-squared threshold for clumping. (Default: %(default)s)")


    # --- Step 1: Trait Processing ---
    trait_group = parser.add_argument_group("Step 1: Trait Processing")
    # trait_group.add_argument("--main_gwas_loci", 
    #                          default="/mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe_independent_markers.tsv",
    #                          metavar="", help="Independent markers for main GWAS. (Default: %(default)s)")
    trait_group.add_argument("--trait_input_file", 
                             default="/mnt/disks/sdd/bnmf-clustering/filtered_vcf_files/trait_gwas.csv",
                             metavar="", help="CSV with paths to trait GWAS files. (Default: %(default)s)")
    trait_group.add_argument("--dbsnp_folder", 
                             default="/mnt/disks/sdd/resourses/postgwas/gwas2vcf/GRCh37/dbSNP/vcf_files/",
                             metavar="", help="Path to dbSNP VCF files. (Default: %(default)s)")
    trait_group.add_argument("--r2_threshold", type=float, default=0.8, 
                             metavar="", help="R-squared threshold for variant correlation. (Default: %(default)s)")
    trait_group.add_argument("--window", type=int, default=1000000, 
                             metavar="", help="Window size (bp) for trait extraction. (Default: %(default)s)")

    # --- Step 3: bNMF Parameters ---
    bnmf_group = parser.add_argument_group("Step 3: bNMF Parameters")
    bnmf_group.add_argument("--n_reps", type=int, default=100, 
                             metavar="", help="Number of repetitions for stability analysis. (Default: %(default)s)")
    bnmf_group.add_argument("--tolerance", type=float, default=1e-6, 
                             metavar="", help="Convergence tolerance for bNMF. (Default: %(default)s)")
    bnmf_group.add_argument("--corr_cutoff", type=float, default=0.8, 
                             metavar="", help="Correlation cutoff for redundant traits. (Default: %(default)s)")

    # Show help and exit if no arguments provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    # Validation: Ensure at least one step is requested
    if not any([args.run_traits, args.run_main, args.run_bnmf]):
        from rich.console import Console
        console = Console()
        console.print("\n[bold red][!] ERROR:[/bold red] You must specify at least one action flag: --run_traits, --run_main, or --run_bnmf.\n")
        parser.print_help(sys.stderr)
        sys.exit(1)

    run_full_pipeline(vars(args))

if __name__ == "__main__":
    main()