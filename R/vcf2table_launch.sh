#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 12G
#SBATCH --time 04:00:00
#SBATCH --job-name=vcf2table
#SBATCH --output=../../data/log/filter_maf04/annotated/annotated_launch_chr%a_%A_%J.out
#SBATCH --error=../../data/log/filter_maf04/annotated/annotated_launch_chr%a_%A_%J.err
#SBATCH --array=1-25
# #SBATCH --array=1
# #SBATCH --array=21-22

# chr21 (.3G gz) test complete at 6-15 minutes, 120G
# chr21 (.3G gz) test fail     at -    minutes,  20G
# chr1 (1.3G gz) test complete at -    minutes, 120G
# chr1split (10MB gz) test complete at -   minutes, 8G

# test file: bcftools_gatk_decom_maf01_chrY.recode_vep_chrY_8661694_9661693.vcf.gz

set -e
echo "START AT $(date)"

source /cluster/work/phrtmma/dlawless/data/variables.sh
# overwrite data variable to work in the specific maf filter directory
DATA="/cluster/work/phrtmma/dlawless/data/filter_maf04"
PREANNOPROCESS="/cluster/work/phrtmma/dlawless/data/pre_annotation_output"
POSTANNOPROCESS="${DATA}/post_annotation_output"

vcf_pre_path=${PREANNOPROCESS}
vcf_out_path=${POSTANNOPROCESS}

# critical step to count number of samples - handle this better
vcf_sample="${vcf_pre_path}/bcftools_gatk_decom_maf04_chr21.recode_vep.vcf.gz"

module load bcftools
sample_count=$(bcftools query -l ${vcf_sample} | wc -l)

# Get specific file based on the array task ID
declare -a total_vcf_files
for i in ${vcf_out_path}/bcftools_gatk_decom_maf01_chr*.recode_vep.vcf.gz; do total_vcf_files+=($i); done
vcf_out=${total_vcf_files[$SLURM_ARRAY_TASK_ID]}

module load StdEnv/2020 r/4.3.1

version="v1"
output_dir="../../data/filter_maf04/annotated/"
file_suffix="phrtmma"

# vcf_file_to_process=${total_vcf_files[${SLURM_ARRAY_TASK_ID}-1]}
# Rscript "vcf2table_run.R" ${version} ${output_dir} ${file_suffix} ${sample_count} ${vcf_file_to_process}
Rscript "vcf2table_run.R" ${version} ${output_dir} ${file_suffix} ${sample_count} ${vcf_out_path} ${SLURM_ARRAY_TASK_ID}

echo "FINISH AT $(date)"
