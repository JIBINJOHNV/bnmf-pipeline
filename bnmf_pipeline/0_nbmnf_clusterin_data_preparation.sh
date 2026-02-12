## converting previously created gcpa_id cluster to tsv and then to vcf

for vcf_path in *build_GRCh37_withEAF.vcf.gz ; do

  n_threads=10
  out_file="${vcf_path%_build_GRCh37_withEAF.vcf.gz}.tsv"

  {
      printf "CHROM\tPOS\tREF\tALT\tID\tPVALUE\tBETA\tSE\tFREQ\tINFO\tN_CONTROLS\n"
      bcftools view --threads ${n_threads} --min-alleles 2 --max-alleles 2 "$vcf_path" | \
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t[%LP]\t[%ES]\t[%SE]\t[%AF]\t[%SI]\t[%SS]\n'
  } > "${out_file}"

done 





default_config='/mnt/disks/sdd/postgwas_analysis/harmonisation.yaml'
resource_folder='/mnt/disks/sdd/resourses/postgwas/gwas2vcf/'
ld_folder='/mnt/disks/sdd/resourses/postgwas/onekg_plinkfiles/GRCh37/LD_ref_EUR/'
base_folder='/mnt/disks/sdd/bnmf-clustering/'
output_folder='/mnt/disks/sdd/bnmf-clustering/final_sumstat_vcf/'
n_thread=10 



docker run --platform=linux/amd64 \
    -u $(id -u):$(id -g) \
    -v /mnt/disks/sdd/:/mnt/disks/sdd/ \
    --rm jibinjv/postgwas:1.3 python /opt/postgwas/src/postgwas/scripts/create_sumstat_map_pl.py \
    --input /mnt/disks/sdd/bnmf-clustering/sumstat/scz/ \
    --output-path ${base_folder}/final_sumstat_vcf/gpca_gwas2vcf_input4.tsv \
    --resource-folder ${resource_folder} \
    --harmonisation-output-path ${base_folder}/final_sumstat_vcf/
 


docker run --platform=linux/amd64 \
   -u $(id -u):$(id -g) \
    -v /mnt/disks/sdd/:/mnt/disks/sdd/ \
    -it jibinjv/postgwas:1.3 postgwas harmonisation \
        --nthreads ${n_thread} \
        --max-mem 50G \
        --config ${base_folder}/final_sumstat_vcf/gpca_gwas2vcf_input4.tsv \
        --defaults ${default_config}


for sample_id in ${output_folder}*/; do
    sample_id="${sample_id%/}"
    sample_id="${sample_id##*/}"

    # Skip if name contains "-1"
    if [[ "$sample_id" == *-1* ]]; then
        continue
    fi

    echo "$sample_id"

    docker run --platform=linux/amd64 \
        -u $(id -u):$(id -g) \
        -v /mnt/disks/sdd/:/mnt/disks/sdd/ \
        -it jibinjv/postgwas:1.3 postgwas sumstat_filter \
        --nthreads ${n_thread} \
        --vcf ${base_folder}/final_sumstat_vcf/${sample_id}/00_harmonised_sumstat/${sample_id}_GRCh37_merged.vcf.gz \
        --sample_id ${sample_id} \
        --outdir ${base_folder}/filtered_vcf_files \
        --maf-cutoff 0.01 \
        --external-af-name EUR \
        --allelefreq-diff-cutoff 0.2 \
        --info-cutoff 0.7 \
        --exclude-palindromic \
        --palindromic-af-lower 0.4 \
        --palindromic-af-upper 0.6 \
        --remove-mhc \
        --mhc-start 25000000 \
        --mhc-end 35000000
done


for sample_id in PGC3_SCZ_wave3_EUR daner_PGC_SCZ_w3_90_0418b_ukbbdedupe pgc_bip2021_BDI ; do 
    mkdir -p ${base_folder}/clustering_input/
    docker run --platform=linux/amd64 \
        -u $(id -u):$(id -g) \
        -v /mnt/disks/sdd/:/mnt/disks/sdd/ \
        -it jibinjv/postgwas:1.3 postgwas ld_clump \
        --ld_clump_population EUR \
        --nthreads ${n_thread} \
        --vcf ${base_folder}/filtered_vcf_files/${sample_id}_GRCh37_filtered.vcf.gz \
        --sample_id ${sample_id} \
        --outdir ${base_folder}/clustering_input/ \
        --ld-folder ${ld_folder} \
        --lead-p 0.00000005 \
        --r2-clump 0.6 \
        --r2-lead 0.1 \
        --merge-dist 250000 
done 













