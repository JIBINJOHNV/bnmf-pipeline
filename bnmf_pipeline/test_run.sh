

bnmf-pipeline \
    --main_gwas_id daner_PGC_SCZ_w3_90_0418b_ukbbdedupe \
    --ld_folder /mnt/disks/sdd/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/ \
    --bcftools_bin /usr/bin/bcftools \
    --n_threads 10 \
    --output_folder /mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe/ \
    --trait_input_file /mnt/disks/sdd/bnmf-clustering/filtered_vcf_files/trait_gwas.csv \
    --dbsnp_folder /mnt/disks/sdd/resourses/postgwas/gwas2vcf/GRCh37/dbSNP/vcf_files/ \
    --r2_threshold 0.8 \
    --window 1000000 \
    --main_gwas_vcf_path /mnt/disks/sdd/bnmf-clustering/filtered_vcf_files/daner_PGC_SCZ_w3_90_0418b_ukbbdedupe_GRCh37_filtered.vcf.gz \
    --snp_to_include /mnt/disks/sdd/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/EUR.chr1_22_Variants_to_include.txt \
    --pop EUR \
    --lead_p 5e-08 \
    --r2_clump 0.05 \
    --n_reps 100 \
    --tolerance 1e-06 \
    --corr_cutoff 0.8 \
    --maximum_k 30 




project_dir="/mnt/disks/sdd/bnmf-clustering/Cognitive_Phenotypes/BIPI/test_run/"
main_gwas_id="pgc_bip2021_BDI"
Rscript /mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/bnmf-pipeline/bnmf_pipeline/R/run_bnmf_cli.R \
    --project_dir ${project_dir} \
    --z_score_file ${project_dir}/pgc_bip2021_BDI_zscore_index.csv \
    --sample_size_file ${project_dir}/pgc_bip2021_BDI_sample_size_index.csv \
    --main_gwas_id ${main_gwas_id} \
    --n_reps 10 \
    --tolerance 1e-06 \
    --corr_cutoff 0.8 \
    --script_path /mnt/disks/sdd/bnmf-clustering/bnmf_cluster_analysis/bnmf-pipeline/bnmf_pipeline/R/ \
    --maximum_k 2
