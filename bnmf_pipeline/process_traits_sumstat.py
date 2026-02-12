import polars as pl
import textwrap
import subprocess
import os 
import time
import psutil
from io import StringIO
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Optional, List, Tuple
import tabix  # pip install pytabix
import os
import time
import polars as pl
from concurrent.futures import ProcessPoolExecutor
from typing import List, Optional



def map_variants_to_rsid_tabix(variant_ids,dbsnp_folder,genome_prefix="GRCh37_dbSNP157_chr"):
    results = {}
    for var in variant_ids:
        chrom,pos,ref,alt = var.split("_")
        region = f"{chrom}:{pos}-{pos}"
        dbsnp_vcf = Path(dbsnp_folder) / f"{genome_prefix}{chrom}.vcf.gz"
        if not dbsnp_vcf.exists():
            results[var] = None
            continue
        try:
            cmd = ["tabix",str(dbsnp_vcf),region]
            result = subprocess.run(cmd,capture_output=True,text=True,check=True)
        except subprocess.CalledProcessError:
            results[var] = None
            continue
        rsid_found = None
        for line in result.stdout.splitlines():
            parts = line.split("\t")
            chrom_db,pos_db,rsid,ref_db,alt_db = parts[:5]
            alt_alleles = alt_db.split(",")
            if ref_db == ref and alt in alt_alleles:
                rsid_found = rsid
                break
        results[var] = rsid_found
    return pl.DataFrame({"uniq_id":list(results.keys()),"rsid":list(results.values())})

# ==============================================================================
# 1. FILTERING BY MISSINGNESS
# ==============================================================================

def filter_by_missingness_polars(
    z_df: pl.DataFrame,
    ss_df: pl.DataFrame,
    variant_missing_cutoff: float = 0.5,
    sample_missing_cutoff: float = 0.5
) -> Tuple[pl.DataFrame, pl.DataFrame, List[str], List[str]]:
    
    trait_cols = [c for c in z_df.columns if c != "uniq_id"]
    n_traits = len(trait_cols)
    # 1. Variant-level filtering
    z_filtered = z_df.with_columns(
        missing_frac = pl.sum_horizontal(pl.all().exclude("uniq_id").is_null()) / n_traits
    ).filter(
        pl.col("missing_frac") < variant_missing_cutoff
    ).drop("missing_frac")
    # Anti join to find variants present in z_df but NOT in z_filtered
    removed_variants = (
        z_df.join(z_filtered.select("uniq_id"), on="uniq_id", how="anti")
        .select("uniq_id")
        .to_series()
        .to_list()
    )
    # 2. Trait-level filtering
    missing_rates = z_filtered.select([
        (pl.col(c).is_null().sum() / pl.len()).alias(c) 
        for c in trait_cols
    ])
    # Extract scalar rate using row 0
    keep_traits = ["uniq_id"] + [
        c for c in trait_cols if missing_rates[c][0] <= sample_missing_cutoff
    ]
    removed_traits = [c for c in trait_cols if c not in keep_traits]
    z_final = z_filtered.select(keep_traits)
    # 3. Align ss_df using a Semi Join
    ss_final = (
        ss_df.join(z_final.select("uniq_id"), on="uniq_id", how="semi")
        .select(keep_traits)
    )
    return z_final, ss_final, removed_variants, removed_traits


# ==============================================================================
# 2. SNP EXTRACTION (TRAIT VCFs)
# ==============================================================================

def _process_single_trait_polars(
    trait_id: str,
    vcf_path: str,
    out_dir: str,
    main_id: str,
    bcf_path: str,
    threads: int
) -> Tuple[str, pl.DataFrame, pl.DataFrame]:
    """Worker for VCF extraction returning Polars DataFrames."""
    output_file = f"{out_dir}/{trait_id}_selected_variants.tsv"
    log_file = f"{out_dir}/{main_id}_snp_extraction_{trait_id}.log"
    loci_file = f"{out_dir}/{main_id}_loci_id_ldvariants.txt"
    cmd = textwrap.dedent(f"""
        {{
            printf "uniq_id\\tSS\\tZ_score\\n"
            {bcf_path} annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' "{vcf_path}" | \\
            {bcf_path} view -i 'ID=@{loci_file}' | \\
            {bcf_path} view --threads {threads} --min-alleles 2 --max-alleles 2 | \\
            {bcf_path} query -f '%ID\\t[%NEF]\\t[%EZ]\\n' | \\
            sed 's|:|_|g'
        }} > "{output_file}"
    """)
    try:
        with open(log_file, "w") as lf:
            subprocess.run(cmd, shell=True, executable="/bin/bash", check=True, stderr=lf)
    except subprocess.CalledProcessError:
        raise RuntimeError(f"Extraction failed for {trait_id}. Check {log_file}")
    df = pl.read_csv(output_file, separator="\t", null_values=["NA", "."])
    z_df = df.select([pl.col("uniq_id"), pl.col("Z_score").alias(trait_id)])
    ss_df = df.select([pl.col("uniq_id"), pl.col("SS").alias(trait_id)])
    return trait_id, z_df, ss_df

def extract_trait_snps_parallel_polars(
    trait_file: str,
    output_folder: str,
    main_gwas_id: str,
    bcftools_path: str = "bcftools",
    bcftools_threads: int = 1,
    max_workers: int = 4
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """Parallel extraction using Polars joins for merging."""
    trait_meta = pl.read_csv(trait_file)
    os.makedirs(output_folder, exist_ok=True)
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _process_single_trait_polars, row["trait_sample_id"], 
                row["trait_input_path"], output_folder, main_gwas_id, 
                bcftools_path, bcftools_threads
            ): row["trait_sample_id"] for row in trait_meta.to_dicts()
        }
        for future in as_completed(futures):
            results.append(future.result())
    # Multi-way join in Polars
    z_final = results[0][1]
    ss_final = results[0][2]
    for _, z, ss in results[1:]:
        z_final = z_final.join(z, on="uniq_id", how="left")
        ss_final = ss_final.join(ss, on="uniq_id", how="left")
    return z_final, ss_final

# ==============================================================================
# 3. LD PROXY SEARCH
# ==============================================================================

def get_ld_partners_tabix(
    variant_id: str,
    ldref_folder: str,
    r2_threshold: float = 0.8,
    window: int = 1_000_000
) -> pl.DataFrame:
    schema = {"query_id": pl.Utf8, "partner_id": pl.Utf8, "r2": pl.Float64}
    try:
        chrom, pos = variant_id.split("_")[:2]
        pos = int(pos)
    except (ValueError, IndexError): return pl.DataFrame(schema=schema)
    ld_path = Path(ldref_folder) / f"EUR_chr{chrom}.ld.gz"
    if not ld_path.exists(): return pl.DataFrame(schema=schema)
    try:
        tb = tabix.open(str(ld_path))
        records = tb.query(chrom, max(1, pos - window), pos + window)
        data = []
        for r in records:
            r2 = float(r[6])
            if r2 >= r2_threshold:
                # LD files columns: chr_a, pos_a, snp_a, chr_b, pos_b, snp_b, r2
                if r[2] == variant_id: data.append((variant_id, r[5], r2))
                elif r[5] == variant_id: data.append((variant_id, r[2], r2))
        return pl.DataFrame(data, schema=schema, orient="row") if data else pl.DataFrame(schema=schema)
    except: return pl.DataFrame(schema=schema)


def _ld_worker(args):
    """
    Worker wrapper for parallel execution.
    Unpacks arguments and calls the core tabix search function.
    """
    # Unpack tuple arguments
    variant_id, ldref_folder, r2_threshold, window = args
    # Ensure this function returns a Polars DataFrame as defined previously
    return get_ld_partners_tabix(variant_id, ldref_folder, r2_threshold, window)


def run_ld_partner_search_parallel(
    missing_variants: List[str],
    ldref_folder: str,
    r2_threshold: float = 0.8,
    window: int = 1_000_000,
    max_workers: Optional[int] = None,
) -> pl.DataFrame:
    schema = {"query_id": pl.Utf8, "partner_id": pl.Utf8, "r2": pl.Float64}
    if not missing_variants:
        return pl.DataFrame(schema=schema)
    # Threads are much lighter than Processes
    if max_workers is None:
        max_workers = 10  # Tabix is fast, 10-20 threads is usually safe
    print(f"[LD] Querying {len(missing_variants)} variants using {max_workers} threads...")
    start_time = time.time()
    # We can call get_ld_partners_tabix directly without a wrapper in threads
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                get_ld_partners_tabix, v, ldref_folder, r2_threshold, window
            ): v for v in missing_variants
        }
        results = [future.result() for future in as_completed(futures)]
    valid_dfs = [df for df in results if df is not None and not df.is_empty()]
    print(f"[LD] Completed in {time.time() - start_time:.2f}s")
    if not valid_dfs:
        return pl.DataFrame(schema=schema)
    return pl.concat(valid_dfs).rechunk()

# ==============================================================================
# 4. IMPUTATION
# ==============================================================================

def impute_missing_with_ld_proxies_polars(
    target_z_df: pl.DataFrame,
    reference_z_df: pl.DataFrame,
    ld_proxy_df: pl.DataFrame
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """
    Vectorized imputation using Polars hash-joins.
    """
    trait_cols = [c for c in target_z_df.columns if c != "uniq_id"]
    # 1. Best proxy per variant
    best_proxies = (
        ld_proxy_df.sort("r2", descending=True)
        .unique(subset=["query_id"], keep="first")
    )
    # 2. Map target IDs to reference values via proxies
    impute_map = (
        target_z_df.select("uniq_id")
        .join(best_proxies, left_on="uniq_id", right_on="query_id", how="left")
        .join(
            reference_z_df.rename({c: f"{c}_imp" for c in trait_cols}),
            left_on="partner_id", right_on="uniq_id", how="left"
        )
    )
    # 3. Coalesce: Fill nulls only where we found a proxy value
    final_df = target_z_df.join(
        impute_map.select(["uniq_id"] + [f"{c}_imp" for c in trait_cols]),
        on="uniq_id", how="left"
    ).with_columns([
        pl.coalesce([pl.col(c), pl.col(f"{c}_imp")]).alias(c)
        for c in trait_cols
    ]).drop([f"{c}_imp" for c in trait_cols])
    return final_df, impute_map.filter(pl.col("partner_id").is_not_null())




def process_traits_gwas_polars(
    main_gwas_id, 
    main_gwas_loci, 
    ld_folder, 
    trait_input_file,
    dbsnp_folder, 
    output_folder,
    r2_threshold=0.8, 
    window=1_000_000, 
    bcftools_bin="bcftools"
):
    # ---------------------------------------------------------
    # 1. Resource Allocation & Setup
    # ---------------------------------------------------------
    os.makedirs(output_folder, exist_ok=True)
    # Dynamic worker calculation
    total_cpu = os.cpu_count() or 1
    free_mem_gb = psutil.virtual_memory().available / (1024**3)
    n_workers = int(min(total_cpu, free_mem_gb // 5))
    n_workers = max(1, n_workers)
    print(f"[*] Starting process with {n_workers} workers.")
    # ---------------------------------------------------------
    # 2. Load Main GWAS & Identify LD Partners
    # ---------------------------------------------------------
    # Load loci and calculate Z-score
    main_gwas_df = pl.read_csv(main_gwas_loci, separator="\t").select([
        "uniq_id", "rsIDcol", "CHR", "POS", "ea", "nea", "eaf", "P_value", "Beta", "SE"
    ]).with_columns(
        z_score = pl.col("Beta") / pl.col("SE")
    )
    main_gwas_variants = main_gwas_df.select("uniq_id").to_series().to_list()
    main_leads = main_gwas_df.select(pl.col("uniq_id").alias("variant_id"))
    # Find LD Proxies
    ld_partners_df = run_ld_partner_search_parallel(
        missing_variants=main_gwas_variants, 
        ldref_folder=ld_folder,
        r2_threshold=r2_threshold, 
        window=window
    )
    # Save the expanded list of variants (Lead + Proxies) for bcftools extraction
    extraction_loci = (
        pl.concat([
            ld_partners_df.select(pl.col("query_id").alias("variant_id")),
            ld_partners_df.select(pl.col("partner_id").alias("variant_id")),
            main_leads
        ])
        .unique()
        .sort("variant_id")
    )
    loci_file_path = f"{output_folder}/{main_gwas_id}_loci_id_ldvariants.txt"
    extraction_loci.write_csv(loci_file_path, include_header=False)
    # ---------------------------------------------------------
    # 3. Parallel Extraction of Traits
    # ---------------------------------------------------------
    zscore_ref_df, samplesize_ref_df = extract_trait_snps_parallel_polars(
        trait_file=trait_input_file,
        output_folder=output_folder,
        main_gwas_id=main_gwas_id,
        bcftools_path=bcftools_bin,
        bcftools_threads=1,
        max_workers=n_workers
    )
    # ---------------------------------------------------------
    # 4. Alignment & Filtering
    # ---------------------------------------------------------
    # Filter the matrices to only include lead variants initially
    zscore_lead_df = main_gwas_df.select("uniq_id").join(zscore_ref_df, on="uniq_id", how="left")
    ss_lead_df = main_gwas_df.select("uniq_id").join(samplesize_ref_df, on="uniq_id", how="left")
    # ---------------------------------------------------------
    # 5. Imputation via LD Proxies (Vectorized Join-Coalesce)
    # ---------------------------------------------------------
    z_full_df, z_log = impute_missing_with_ld_proxies_polars(
        target_z_df=zscore_lead_df, 
        reference_z_df=zscore_ref_df, 
        ld_proxy_df=ld_partners_df
    )
    ss_full_df, ss_log = impute_missing_with_ld_proxies_polars(
        target_z_df=ss_lead_df, 
        reference_z_df=samplesize_ref_df, 
        ld_proxy_df=ld_partners_df
    )
    z_log.write_csv(f"{output_folder}/{main_gwas_id}_variants_missing_zscore_imputed.csv")
    ss_log.write_csv(f"{output_folder}/{main_gwas_id}_variants_missing_samplesize_imputed.csv")
    # Apply Missingness Filter (Polars optimized)
    z_full_filt, ss_full_filt, rem_vars, rem_traits = filter_by_missingness_polars(
        z_df=z_full_df,
        ss_df=ss_full_df,
        variant_missing_cutoff=0.1,
        sample_missing_cutoff=0.1
    )
    print(f"number of variants removed {len(rem_vars)} : {rem_vars}")
    print("")
    print(f"number of traits removed {len(rem_traits)} : {rem_traits}")
    # ---------------------------------------------------------
    # 6. Final Median Imputation & Formatting
    # ---------------------------------------------------------
    trait_cols = [c for c in z_full_filt.columns if c != "uniq_id"]
    # Fill remaining NaNs with column medians using Polars
    z_final = z_full_filt.with_columns([
        pl.col(c).fill_null(pl.col(c).median()) for c in trait_cols
    ])
    ss_final = ss_full_filt.with_columns([
        pl.col(c).fill_null(pl.col(c).median()) for c in trait_cols
    ])
    # Create chr:pos and Allele columns
    z_final = z_final.with_columns(
        chr_pos = pl.col("uniq_id").str.split("_").list.get(0) + ":" + pl.col("uniq_id").str.split("_").list.get(1)
    )
    ss_final = ss_final.with_columns(
        chr_pos = pl.col("uniq_id").str.split("_").list.get(0) + ":" + pl.col("uniq_id").str.split("_").list.get(1)
    )
    # ---------------------------------------------------------
    # 7. Map rsIDs and Save Outputs
    # ---------------------------------------------------------
    # (Using the previous map function, but converting list to Polars)
    rsid_map_df = map_variants_to_rsid_tabix(
        variant_ids=z_final.select("uniq_id").to_series().to_list(),
        dbsnp_folder=dbsnp_folder
    )
    rsid_map_pl = pl.DataFrame(rsid_map_df.to_dict(as_series=False)).rename({"uniq_id": "VAR_ID", "rsid": "rsID"})
    rsid_map_pl.write_csv(f"{output_folder}/rsID_map.txt", separator=" ")
    # Save Indices
    zscore_file=f"{output_folder}/{main_gwas_id}_zscore_index.csv"
    sample_size_file=f"{output_folder}/{main_gwas_id}_sample_size_index.csv"
    z_final.write_csv(zscore_file)
    ss_final.write_csv(sample_size_file)
    # Alignment File (Risk Allele metadata)
    alignment_df = (
        z_final.select(["uniq_id", "chr_pos"])
        .with_columns([
            pl.col("uniq_id").str.split("_").list.get(2).alias("REF"),
            pl.col("uniq_id").str.split("_").list.get(3).alias("ALT")
        ])
        .with_columns([
            pl.col("ALT").alias("Risk_Allele_Orig"),
            pl.col("ALT").alias("Risk_Allele"),
            pl.col("REF").alias("Nonrisk_Allele")
        ])
        .join(main_gwas_df.select(["uniq_id", "P_value", "Beta", "SE", "z_score"]), on="uniq_id")
        .rename({"chr_pos": "SNP", "P_value": "P_VALUE", "Beta": "BETA"})
        .drop("uniq_id")
    )
    output_filename=f"{output_folder}/alignment_GWAS_summStats.csv"
    alignment_df.write_csv(output_filename)
    print(f"[SUCCESS] Pipeline completed. Results in {output_folder}")
    return {
        "zscore_file": zscore_file,
        "ss_final": ss_final
    }


# ld_folder='/mnt/disks/sdd/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/'
# main_gwas_id='daner_PGC_SCZ_w3_90_0418b_ukbbdedupe'
# main_gwas_loci='/mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe_independent_markers.tsv'

# output_folder=f'/mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/{main_gwas_id}/'
# input_loci_folder='/mnt/disks/sdd/bnmf-clustering/clustering_input/'
# os.system(f"mkdir -p {output_folder}")
# bcftools_path='/usr/bin/bcftools'
# n_threads=5
# r2_threshold=0.8
# variant_missing_cut_of=0.8
# sample_missing_cut_of=0.8
# dbsnp_folder="/mnt/disks/sdd/resourses/postgwas/gwas2vcf/GRCh37/dbSNP/vcf_files/"

# trait_input_file='/mnt/disks/sdd/bnmf-clustering/filtered_vcf_files/trait_gwas.csv'
# window=1_000_000, 


# process_traits_gwas_polars(
#     main_gwas_id, 
#     main_gwas_loci, 
#     ld_folder, 
#     trait_input_file,
#     dbsnp_folder, 
#     output_folder,
#     r2_threshold=0.8, 
#     window=1_000_000, 
#     bcftools_bin="bcftools"
# )