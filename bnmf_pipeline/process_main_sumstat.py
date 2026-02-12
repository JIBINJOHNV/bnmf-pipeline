import polars as pl
import psutil
import os
import math
import subprocess
from pathlib import Path
import textwrap
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
from packaging import version
import pandas as pd


def vcf_to_standered_ldclump(sumstat_vcf: str, output_folder: str, sample_name: str, bcftools_path="bcftools", n_threads=4):
    """
    Converts VCF to TSV using specified bcftools path with multithreading support 
    and detailed command logging.
    """
    vcf_path = Path(sumstat_vcf)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{sample_name}_formatted.tsv"
    log_file = output_dir / f"{sample_name}_vcf_to_standered_ldclump_conversion.log"
    # Constructing the bash command
    cmd = textwrap.dedent(f"""
        {{
            printf "chrcol\\tposcol\\tneacol\\teacol\\trsIDcol\\tpcol\\tbecol\\tsecol\\teafcol\\n"
            {bcftools_path} view --threads {n_threads} --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
            {bcftools_path} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t[%LP]\\t[%ES]\\t[%SE]\\t[%AF]\\n' | \
            sed 's|:|_|g' | \
            awk -F '\\t' 'BEGIN {{ OFS="\\t" }} {{
                raw_p = exp(-$6 * log(10))
                $6 = sprintf("%.6g", raw_p)
                if ($7 == "" || $7 == ".") $7 = "NA"
                if ($8 == "" || $8 == ".") $8 = "NA"
                print
            }}'
        }} > "{output_file}"
    """)
    try:
        with open(log_file, "w") as lf:
            # Writing metadata and the exact command to the log
            lf.write(f"--- VCF Conversion Log ---\n")
            lf.write(f"Sample: {sample_name}\n")
            lf.write(f"Threads: {n_threads}\n")
            lf.write(f"Command used:\n{cmd}\n")
            lf.write(f"--- Execution Output ---\n")
            lf.flush() # Ensure the command is written before the process starts
            # Execute the command
            subprocess.run(cmd, shell=True, executable="/bin/bash", check=True, stdout=lf, stderr=lf)
            lf.write(f"\n[STAMP] Conversion completed successfully.\n")
    except subprocess.CalledProcessError:
        raise RuntimeError(f"VCF conversion failed for {sample_name}. Check log for details: {log_file}")
    return str(output_file)

def get_ld_partners(chrom, pos, rsid, ld_file_path, r2_threshold):
    # Initial state: variant only paired with itself
    partners_r2 = {rsid: 1.0}
    
    start = max(0, pos - 1_000_000)
    end   = pos + 1_000_000
    region = f"{chrom}:{start}-{end}"
    # Version check for Polars
    import polars as pl
    from packaging import version
    pl_version = version.parse(pl.__version__)
    csv_kwargs = {"separator": "\t"} if pl_version >= version.parse("0.19.2") else {"sep": "\t"}
    
    try:
        # 1. Run Tabix
        result = subprocess.run(
            ["tabix", ld_file_path, region],
            capture_output=True,
            text=True,
            check=True
        )
        # STAGE 1: If tabix returns nothing for the region
        if not result.stdout.strip():
            return partners_r2, "NOT_IN_REFERENCE"
        ld_df = pl.read_csv(
                StringIO(result.stdout),
                has_header=False,
                new_columns=["chr_a", "pos_a", "snp_a",
                            "chr_b", "pos_b", "snp_b", "r2"],
                schema_overrides={
                    "chr_a": pl.Utf8,
                    "pos_a": pl.Int64,
                    "snp_a": pl.Utf8,
                    "chr_b": pl.Utf8,
                    "pos_b": pl.Int64,
                    "snp_b": pl.Utf8,
                    "r2": pl.Float64,
                },
                **csv_kwargs
            )
        # STAGE 2: Check if the rsID exists in the reference data at all
        exists_in_ref = ld_df.filter((pl.col("snp_a") == rsid) | (pl.col("snp_b") == rsid))
        
        if exists_in_ref.is_empty():
            return partners_r2, "NOT_IN_REFERENCE"
        # STAGE 3: Filter by R2 threshold
        filtered = exists_in_ref.filter(pl.col("r2") >= r2_threshold).with_columns(
            pl.when(pl.col("snp_a") == rsid)
            .then(pl.col("snp_b"))
            .otherwise(pl.col("snp_a"))
            .alias("partner_id")
        )
        if filtered.height > 1: # height > 1 because the self-match is often in the file
            new_partners = dict(zip(filtered["partner_id"], filtered["r2"]))
            partners_r2.update(new_partners)
            return partners_r2, "LD_CLUMPED"
        else:
            return partners_r2, "NO_LD_PARTNERS"
    except subprocess.CalledProcessError as e:
        return partners_r2, f"TABIX_ERROR: {e}"
    except Exception as e:
        return partners_r2, f"PARSING_ERROR: {e}"
    

def find_ind_sig_snps(chr_df, ld_path, lead_p_threshold, r2_clump_threshold, n_threads=1):
    # Filter for significant leads and sort by p-value (lowest first)
    remaining_leads = chr_df.filter(pl.col("pcol") <= lead_p_threshold).sort("pcol")
    available_pool = chr_df.clone()
    ind_sig_clumps = []
    if chr_df.is_empty():
        return pl.DataFrame()
    chrom = chr_df.select(pl.col("chrcol").first()).item()
    
    while not remaining_leads.is_empty():
        # 1. Pick the most significant remaining SNP as the Lead
        top_snp = remaining_leads.row(0, named=True)
        sid = top_snp['uniq_id']
        # 2. Get LD partners AND the status (Not in ref vs No LD)
        # Assumes get_ld_partners returns: (dict, str)
        partners_map, status = get_ld_partners(
            top_snp['chrcol'], top_snp['poscol'], sid, ld_path, r2_clump_threshold
        )
        candidate_ids = list(partners_map.keys())
        # 3. Identify which SNPs in our GWAS pool are in LD with this lead
        # Note: 'sid' is always in candidate_ids, so members will at least contain the lead
        members = available_pool.filter(pl.col("uniq_id").is_in(candidate_ids))
        if not members.is_empty():
            # Create a mapping DF for R2 values and the status
            r2_df = pl.DataFrame({
                "uniq_id": list(partners_map.keys()), 
                "r2_with_IndSig": list(partners_map.values())
            })
            # Join and annotate the clump members
            members = (
                members
                .join(r2_df, on="uniq_id", how="left")
                .with_columns([
                    pl.lit(sid).alias("ind_sig_SNP_id"),
                    pl.lit(len(partners_map)).alias("n_refsnps"),
                    pl.lit(status).alias("clump_status") # New column to distinguish
                ])
            )
            ind_sig_clumps.append(members)
            # 4. Remove assigned SNPs from both the pool and lead list
            assigned_ids = members["uniq_id"].to_list()
            available_pool = available_pool.filter( ~pl.col("uniq_id").is_in(assigned_ids) )
            remaining_leads = remaining_leads.filter( ~pl.col("uniq_id").is_in(assigned_ids) )
        else:
            # Fallback: This only triggers if something went wrong with the ID matching
            # We still remove the lead to avoid an infinite loop
            remaining_leads = remaining_leads.filter(pl.col("uniq_id") != sid)
        #print(f"[CHR {chrom}] Variants remaining: {remaining_leads.height}")
    print(f"[CHR {chrom}] Independent variant detection completed.")
    if not ind_sig_clumps:
        return pl.DataFrame()
    return pl.concat(ind_sig_clumps)


def process_chromosome(chrom, gwas_subset, ld_folder, pop, lead_p, r2_clump):
    """
    Worker function to process a single chromosome.
    Returns a Polars DataFrame if clumps are found, or an empty DataFrame with 
    the correct schema if no hits exist for that chromosome.
    """
    ld_path = os.path.join(ld_folder, f"{pop}_chr{chrom}.ld.gz")
    # 1. Handle missing LD reference files
    if not os.path.exists(ld_path):
        print(f"[!] Warning: LD file not found for Chromosome {chrom} at {ld_path}")
        return None
    # 2. Handle empty GWAS subset (no variants for this chrom)
    if gwas_subset.is_empty():
        return pl.DataFrame()
    try:
        # Step A: Find Independent Significant SNPs
        ind_sig = find_ind_sig_snps(gwas_subset, ld_path, lead_p, r2_clump)
        # 3. Ensure we return None or Empty consistently to the thread collector
        if ind_sig is None or ind_sig.is_empty():
            return pl.DataFrame()
        return ind_sig
    except Exception as e:
        print(f"[ERROR] Failed processing Chromosome {chrom}: {e}")
        return None


def run_parallel_clumping(
    gwas_df,
    main_gwas_id,
    ld_folder,
    pop,
    lead_p,
    r2_clump,
    outdir,
    snp_to_include=None
):
    """
    Orchestrates parallel LD clumping.
    Handles both Polars DataFrames and file paths as input.
    Fully Polars-version safe.
    """
    import polars as pl
    from packaging import version
    pl_ver = version.parse(pl.__version__)
    is_new_pl = pl_ver >= version.parse("0.19.0")
    csv_kwargs = {"separator": "\t"} if pl_ver >= version.parse("0.19.2") else {"sep": "\t"}
    # ------------------------------------------------------------------
    # 1. Load file if path provided
    # ------------------------------------------------------------------
    if isinstance(gwas_df, (str, os.PathLike)):
        print(f"[*] Input is a path. Loading GWAS data from {gwas_df}...")
        gwas_df = pl.read_csv(gwas_df, **csv_kwargs)
    # Ensure uniq_id exists
    if "uniq_id" not in gwas_df.columns:
        gwas_df = gwas_df.with_columns(
            pl.concat_str(
                [
                    pl.col("chrcol").cast(pl.Utf8),
                    pl.col("poscol").cast(pl.Utf8),
                    pl.col("neacol"),
                    pl.col("eacol"),
                ],
                separator="_",
            ).alias("uniq_id")
        )
    # ------------------------------------------------------------------
    # 2. Optional SNP filtering by reference include file
    # ------------------------------------------------------------------
    if snp_to_include is not None:
        snp_to_include_df = pl.read_csv(
            snp_to_include,
            has_header=False,
            new_columns=["uniq_id"],
            **csv_kwargs,
        )
        gwas_not_present = gwas_df.join(
            snp_to_include_df, on="uniq_id", how="anti"
        )
        gwas_not_present.write_csv(
            f"{outdir}/{main_gwas_id}_variants_not_present_in_reference.csv"
        )
        gwas_df = gwas_df.join(
            snp_to_include_df, on="uniq_id", how="semi"
        )
    # ------------------------------------------------------------------
    # 3. Balanced Chromosome Scheduling (FIXED)
    # ------------------------------------------------------------------
    if is_new_pl:
        chrom_stats = (
            gwas_df.group_by("chrcol")
            .len()
            .rename({"len": "count"})
        )
    else:
        chrom_stats = gwas_df.groupby("chrcol").count()
    chrom_stats = chrom_stats.sort("count", descending=True)
    balanced_chroms = []
    big_to_small = chrom_stats["chrcol"].to_list()
    while big_to_small:
        balanced_chroms.append(big_to_small.pop(0))
        if big_to_small:
            balanced_chroms.append(big_to_small.pop(-1))
    # ------------------------------------------------------------------
    # 4. Dynamic Worker Allocation
    # ------------------------------------------------------------------
    gwas_mem_gb = gwas_df.estimated_size() / (1024 ** 3)
    dynamic_gb_per_thread = max(1.0, gwas_mem_gb * 3)
    total_cpu = os.cpu_count() or 1
    free_mem_gb = psutil.virtual_memory().available / (1024 ** 3)
    n_workers = int(min(total_cpu, free_mem_gb // dynamic_gb_per_thread))
    n_workers = max(1, n_workers)
    print(
        f"[*] Resource Check: {free_mem_gb:.2f}GB free. "
        f"Using {n_workers} workers for {len(balanced_chroms)} chromosomes."
    )
    # ------------------------------------------------------------------
    # 5. Parallel Execution
    # ------------------------------------------------------------------
    res_list = []
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(
                process_chromosome,
                chrom,
                gwas_df.filter(pl.col("chrcol") == chrom),
                ld_folder,
                pop,
                lead_p,
                r2_clump,
            ): chrom
            for chrom in balanced_chroms
        }
        for future in as_completed(futures):
            chrom = futures[future]
            try:
                result = future.result()
                if result is not None and not result.is_empty():
                    res_list.append(result)
                    #print(f"[SUCCESS] Chromosome {chrom} completed.")
            except Exception as e:
                print(f"[!] Thread Error on chromosome {chrom}: {e}")
    # ------------------------------------------------------------------
    # 6. Final Concatenation
    # ------------------------------------------------------------------
    if not res_list:
        print("[!] Clumping resulted in zero variants.")
        return pl.DataFrame()
    return pl.concat(res_list)

    

def process_main_gwas(
    main_gwas_vcf_path,
    outdir,
    main_gwas_id,
    ld_folder,
    bcftools_bin="bcftools",
    n_threads=2,
    pop="EUR",
    lead_p=5e-8,
    r2_clump=0.05,
    snp_to_include=None,   # FIXED
):
    tsv_path = vcf_to_standered_ldclump(
        main_gwas_vcf_path,
        outdir,
        main_gwas_id,
        bcftools_bin,
        n_threads=n_threads,
    )
    ind_snp_df = run_parallel_clumping(
        gwas_df=tsv_path,
        main_gwas_id=main_gwas_id,
        ld_folder=ld_folder,
        pop=pop,
        lead_p=lead_p,
        r2_clump=r2_clump,
        outdir=outdir,
        snp_to_include=snp_to_include,
    )
    if ind_snp_df.is_empty():
        print("[!] No independent SNPs detected.")
        return None
    output_df = (
            ind_snp_df
            .filter(pl.col("uniq_id") == pl.col("ind_sig_SNP_id"))
            .select([
                pl.col("chrcol").alias("CHR"),
                pl.col("poscol").alias("POS"),
                pl.col("uniq_id"),
                pl.col("rsIDcol"),
                pl.col("neacol").alias("nea"),
                pl.col("eacol").alias("ea"),
                pl.col("pcol").alias("P_value"),
                pl.col("becol").alias("Beta"),
                pl.col("secol").alias("SE"),
                pl.col("eafcol").alias("eaf"),
            ])
            .sort(["CHR", "POS"])
        )
    output_path = f"{outdir}/{main_gwas_id}_independent_markers.tsv"
    output_df.write_csv(output_path, separator="\t")
    return {"main_sumstat_tsv": output_path}




# snp_to_include='/mnt/disks/sdd/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/EUR.chr1_22_Variants_to_include.txt'
# vcf_path='/mnt/disks/sdd/bnmf-clustering/filtered_vcf_files/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe_GRCh37_filtered.vcf.gz'
# sample_id='daner_PGC_SCZ_w3_90_0418b_ukbbdedupe'
# outdir=f'/mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/{sample_id}/'
# bcftools_bin='/usr/bin/bcftools'
# n_threads=10
# ld_folder='/mnt/disks/sdd/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/'
# pop="EUR"
# lead_p=5e-8
# r2_clump=0.05
# merge_dist=250000
                 

# process_main_gwas(
#     main_gwas_vcf_path=vcf_path,
#     outdir=outdir,
#     main_gwas_id=sample_id,
#     ld_folder=ld_folder,
#     bcftools_bin=bcftools_bin,
#     n_threads=10,
#     pop="EUR",
#     lead_p=5e-8,
#     r2_clump=0.05,
#     snp_to_include=snp_to_include,
# )
