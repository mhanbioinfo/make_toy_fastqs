#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./
#SBATCH --mem=16G
#SBATCH -J seqtk_sample_fastq
#SBATCH -p all
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o ./%j-%x.out
#SBATCH -e ./%j-%x.err

module load seqtk

echo "Job started at "$(date) 
time1=$(date +%s)

#IN_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/make_toy_fastqs_outputs/outputs/output_allchrs/toy_sample2"
#FQ1=${IN_DIR}/"toy_sample2_all_shuf222_F19K16_F24B22_BBBBBBBB.R1.fastq.gz"
#FQ2=${IN_DIR}/"toy_sample2_all_shuf222_F19K16_F24B22_BBBBBBBB.R2.fastq.gz"

IN_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/make_toy_fastqs_outputs/outputs/output_allchrs/toy_sample1"
FQ1=${IN_DIR}/"toy_sample1_all_shuf111_F19K16_F24B22_AAAAAAAA.R1.fastq.gz"
FQ2=${IN_DIR}/"toy_sample1_all_shuf111_F19K16_F24B22_AAAAAAAA.R2.fastq.gz"

FQ_ARR=($FQ1 $FQ2)

OUT_DIR="${IN_DIR}/seqtk_sampled"
mkdir -p ${OUT_DIR}
NUM_READS=50000

module load seqtk

#for FILE in ${IN_DIR}/*fastq.gz; do
for FILE in "${FQ_ARR[@]}"; do
    F_NAME="$(basename ${FILE} .fastq.gz)"
    echo "processing ${F_NAME}..."
    seqtk sample -s42 ${FILE} ${NUM_READS} | gzip -c - > ${OUT_DIR}/${F_NAME}.50k.fastq.gz
done


time2=$(date +%s)
echo "Job ended at "$(date)
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF
