import pandas as pd
import polars as pl
import textwrap
import subprocess
import os,glob
from io import StringIO


sample_ids=['164IDP_5cluster_08Oct5_cluster_0_9_IDPs', '459IDP_6cluster_11Oct_luster_3_44_IDPs', 
                  '164IDP_5cluster_08Oct5_cluster_2_15_IDPs', '164IDP_5cluster_08Oct5_cluster_1_13_IDPs', 
                  '2315IDP_Cluster_6_129_IDPs', '164IDP_5cluster_08Oct5_cluster_3_14_IDPs', 
                  '459IDP_6cluster_11Oct_cluster_0_24_IDPs', '459IDP_6cluster_11Oct_cluster_5_31_IDPs', 
                  '2315IDP_Cluster_14_229_IDPs', '2315IDP_Cluster_4_121_IDPs', '2315IDP_Cluster_2_78_IDPs', 
                  '2315IDP_Cluster_7_186_IDPs',  '2315IDP_Cluster_1_176_IDPs', '2315IDP_Cluster_10_67_IDPs', 
                  '459IDP_6cluster_11Oct_cluster_2_31_IDPs', '2315IDP_Cluster_12_277_IDPs', 
                  '2315IDP_Cluster_5_143_IDPs', '164IDP_5cluster_08Oct5_cluster_4_11_IDPs', 
                  '459IDP_6cluster_11Oct_cluster_4_272_IDPs', '2315IDP_Cluster_13_47_IDPs', '2315IDP_Cluster_3_42_IDPs', 
                  '2315IDP_Cluster_0_68_IDPs', '2315IDP_Cluster_9_66_IDPs',  '2315IDP_Cluster_8_74_IDPs', 
                  '459IDP_6cluster_11Oct_cluster_1_41_IDPs','Lam_et_al_2021_CognitiveTaskPerformance',
                  'CRP_2024_PMC9033829_GCST90029070_buildGRCh37', 'MetS_2024_GCST90444487_PMC11549047',
                  'EA4_additive_excl_23andMe','PGC3_SCZ_wave3_EUR','daner_PGC_SCZ_w3_90_0418b_ukbbdedupe','pgc_bip2021_BDI']

output_folder='/mnt/disks/sdd/bnmf-clustering/clustering_input/traits_inpts/'
os.system(f"mkdir -p {output_folder}")
bcftools_path='/usr/bin/bcftools'
n_threads=5


for sample_id in sample_ids:
    vcf_path='/mnt/disks/sdd/bnmf-clustering/filtered_vcf_files/'
    vcf_path=f'{vcf_path}/{sample_id}_GRCh37_filtered.vcf.gz'
    raw_output_file=f'{output_folder}/{sample_id}_bnmf_clustering_input.txt'
    log_file = f'{output_folder}/{sample_id}_bnmf_clustering_input.log'
    cmd = textwrap.dedent(f"""
            set -o pipefail
            (
                printf "VAR_ID\\tEffect_Allele_PH\\tP_VALUE\\tID\\tBETA\\tSE\\tN_PH\\taf_meta\\n"
                {bcftools_path} view --min-alleles 2 --max-alleles 2 "{vcf_path}" | \
                bcftools query -f '%CHROM:%POS:%REF:%ALT\\t%ALT\\t[%LP]\\t%ID\\t[%ES]\\t[%SE]\\t[%NEF]\\t[%AF]\\n' | \
                sed 's|:|_|g'  | \
                    awk -F '\\t' 'BEGIN {{ OFS="\\t" }} {{
                        raw_p = exp(-$3 * log(10))
                        $3 = sprintf("%.6g", raw_p)
                        print
                    }}'
            ) > "{raw_output_file}"
        """)
    try:
        with open(log_file, "w") as lf:
            # Writing metadata and the exact command to the log
            lf.write(f"--- VCF Conversion Log ---\n")
            lf.write(f"Sample: {sample_id}\n")
            lf.write(f"Threads: {n_threads}\n")
            lf.write(f"Command used:\n{cmd}\n")
            lf.write(f"--- Execution Output ---\n")
            lf.flush() # Ensure the command is written before the process starts
            # Execute the command
            subprocess.run(cmd, shell=True, executable="/bin/bash", check=True, stdout=lf, stderr=lf)
            lf.write(f"\n[STAMP] Conversion completed successfully.\n")
    except subprocess.CalledProcessError:
        raise RuntimeError(f"VCF conversion failed for {sample_id}. Check log for details: {log_file}")






files = glob.glob(f'{output_folder}/*_bnmf_clustering_input.txt')
variant_rsid_df = pd.DataFrame(columns=["VAR_ID", "ID"])

for f in files:
    print(f"Reading {f}")
    df = pd.read_csv(f, sep="\t", usecols=["VAR_ID", "ID"])
    # Append
    variant_rsid_df = pd.concat(
        [variant_rsid_df, df],
        ignore_index=True
    )
    variant_rsid_df = variant_rsid_df.drop_duplicates()

# Remove duplicates
variant_rsid_df = variant_rsid_df.drop_duplicates()
variant_rsid_df.columns=['rsID','VAR_ID']

print("Final shape:", variant_rsid_df.shape)

variant_rsid_df.to_csv(f'{output_folder}/rsID_map_example.txt',sep="\t",index=None)






gsutil -m cp /mnt/disks/sdd/bnmf-clustering/clustering_input/traits_inpts/pgc_bip2021_BDI_bnmf_clustering_input.txt \
  "gs://gwas_summary_statistics/bnmf_clustering" 