#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./logs_slurm
#SBATCH --mem=32G
#SBATCH -J get_mock_fastqs_4_testing
#SBATCH -p himem
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -o ./%j-%x.out
#SBATCH -e ./%j-%x.err

echo "Job started at "$(date) 
time1=$(date +%s)

DATA_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/data"
OUT_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/get_mock_fastqs_4_testing/outputs"
TMP_DIR="${OUT_DIR}/TMP"
mkdir -p $TMP_DIR

## lines in BEDPE, so 400K*2reads*3samples=2.5m reads
NUM_LINES=400000


echo "Pivot fastq files."
## pivot wide FASTQs to get Phred
SAMPLE1_FASTQ_R1_FNAME="CMP-01-02-cfDNA-03_lib1_R1.fastq.gz"
SAMPLE1_FASTQ_R2_FNAME="CMP-01-02-cfDNA-03_lib1_R2.fastq.gz"
SAMPLE2_FASTQ_R1_FNAME="CMP-01-02-cfDNA-04_lib1_R1.fastq.gz"
SAMPLE2_FASTQ_R2_FNAME="CMP-01-02-cfDNA-04_lib1_R2.fastq.gz"
SAMPLE3_FASTQ_R1_FNAME="CMP-01-02-cfDNA-05_lib1_R1.fastq.gz"
SAMPLE3_FASTQ_R2_FNAME="CMP-01-02-cfDNA-05_lib1_R2.fastq.gz"

SAMPLES_FASTQ_ARR=( ${SAMPLE1_FASTQ_R1_FNAME} ${SAMPLE1_FASTQ_R2_FNAME}
                    ${SAMPLE2_FASTQ_R1_FNAME} ${SAMPLE2_FASTQ_R2_FNAME}
                    ${SAMPLE3_FASTQ_R1_FNAME} ${SAMPLE3_FASTQ_R2_FNAME} )

#for SAMPLE_FASTQ_R_FNAME in "${SAMPLES_FASTQ_ARR[@]:0:2}" ; do
for SAMPLE_FASTQ_R_FNAME in "${SAMPLES_FASTQ_ARR[@]}" ; do
    echo "Processing ${SAMPLE_FASTQ_R_FNAME}"
    
    SAMPLE_FASTQ_R_PIVOT_FNAME="${SAMPLE_FASTQ_R_FNAME%.fastq.gz}.fastq.pivot.sortd"
    zcat ${DATA_DIR}/${SAMPLE_FASTQ_R_FNAME} \
        | paste -d$'\t' - - - - \
        | sed 's/ /\t/' \
        | sort -k1,1 \
        > ${TMP_DIR}/${SAMPLE_FASTQ_R_PIVOT_FNAME}
done


echo "Grab Chr Start End FragName columns from BEDPE file."
## grab chr start end from bedpe.gz
SAMPLE1_BEDPE_FNAME="CMP-01-02-cfDNA-03.bedpe.gz"
SAMPLE2_BEDPE_FNAME="CMP-01-02-cfDNA-04.bedpe.gz"
SAMPLE3_BEDPE_FNAME="CMP-01-02-cfDNA-05.bedpe.gz"
SAMPLES_BEDPE_ARR=( ${SAMPLE1_BEDPE_FNAME} ${SAMPLE2_BEDPE_FNAME} ${SAMPLE3_BEDPE_FNAME} )

#for SAMPLE_FNAME in "${SAMPLES_BEDPE_ARR[@]:0:1}" ; do
for SAMPLE_FNAME in "${SAMPLES_BEDPE_ARR[@]}" ; do
    echo "Processing ${SAMPLE_FNAME}"
    
    SAMPLE_CHR21_SHUF_BED_FNAME="${SAMPLE_FNAME%.bedpe.gz*}.chr21.shuf.bedpe"
    zcat ${DATA_DIR}/${SAMPLE_FNAME} \
        | awk 'BEGIN {OFS="\t"} ($1 == "chr21" && $4 == "chr21") {print $1,$2,$3,$4,$5,$6,"@"$7}' \
        | shuf -n ${NUM_LINES} \
        > ${TMP_DIR}/${SAMPLE_CHR21_SHUF_BED_FNAME}

    SAMPLE_CHR21_SHUF_BED_R1_FNAME="${SAMPLE_CHR21_SHUF_BED_FNAME%.*}.R1.bed.sortd"
    SAMPLE_CHR21_SHUF_BED_R2_FNAME="${SAMPLE_CHR21_SHUF_BED_FNAME%.*}.R2.bed.sortd"
    cat ${TMP_DIR}/${SAMPLE_CHR21_SHUF_BED_FNAME} \
        | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$7}' \
        | sed 's/|/\t/' | sed 's/\./\t/' \
        | sort -k4,4 \
        > ${TMP_DIR}/${SAMPLE_CHR21_SHUF_BED_R1_FNAME}
    cat ${TMP_DIR}/${SAMPLE_CHR21_SHUF_BED_FNAME} \
        | awk 'BEGIN {OFS="\t"} {print $4,$5,$6,$7}' \
        | sed 's/|/\t/' | sed 's/\./\t/' \
        | sort -k4,4 \
        > ${TMP_DIR}/${SAMPLE_CHR21_SHUF_BED_R2_FNAME}
done


echo "Joining pivoted sorted FASTQ and sorted BED, then convert to fastqs."
## join pivoted FASTQ and chr star end
## and then
## Rscript to convert .joined to fastqs
module load R/4.1
## need packages: tidyverse, ggplot2, GenomicFeatures, BSgenome.Hsapiens.UCSC.hg38

RSCRIPT_PATH="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/get_mock_fastqs_4_testing/shuffled_joined_to_fastq_v1.R"
SAMPLE_INDEX="ACTGACTG:ACTGACTG"

#for SAMPLE_FNAME in "${SAMPLES_BEDPE_ARR[@]:0:1}" ; do
for SAMPLE_FNAME in "${SAMPLES_BEDPE_ARR[@]}" ; do
    echo "Processing ${SAMPLE_FNAME%.bedpe.gz*}"

    SAMPLE_CHR21_SHUF_BED_R1_FNAME="${SAMPLE_FNAME%.bedpe.gz*}.chr21.shuf.R1.bed.sortd"
    SAMPLE_CHR21_SHUF_BED_R2_FNAME="${SAMPLE_FNAME%.bedpe.gz*}.chr21.shuf.R2.bed.sortd"

    SAMPLE_FASTQ_R_PIVOT_R1_FNAME="${SAMPLE_FNAME%.bedpe.gz*}_lib1_R1.fastq.pivot.sortd"
    SAMPLE_FASTQ_R_PIVOT_R2_FNAME="${SAMPLE_FNAME%.bedpe.gz*}_lib1_R2.fastq.pivot.sortd"

    SAMPLE_R1_JOINED="${SAMPLE_FNAME%.bedpe.gz*}.chr21.shuf.R1.joined"
    SAMPLE_R2_JOINED="${SAMPLE_FNAME%.bedpe.gz*}.chr21.shuf.R2.joined"

    join -t $'\t' -1 1 -2 4 \
        ${TMP_DIR}/${SAMPLE_FASTQ_R_PIVOT_R1_FNAME} \
        ${TMP_DIR}/${SAMPLE_CHR21_SHUF_BED_R1_FNAME} \
        > ${TMP_DIR}/${SAMPLE_R1_JOINED}
    join -t $'\t' -1 1 -2 4 \
        ${TMP_DIR}/${SAMPLE_FASTQ_R_PIVOT_R2_FNAME} \
        ${TMP_DIR}/${SAMPLE_CHR21_SHUF_BED_R2_FNAME} \
        > ${TMP_DIR}/${SAMPLE_R2_JOINED}

    Rscript ${RSCRIPT_PATH} \
        --R1_fpath ${TMP_DIR}/${SAMPLE_R1_JOINED} \
        --R2_fpath ${TMP_DIR}/${SAMPLE_R2_JOINED} \
        --output_dir ${TMP_DIR} \
        --sample_index ${SAMPLE_INDEX}
done


echo "Concatenating fastqs from all samples subsampled, for R1, then R2."
## cat all fastqs for R1, then R2

R1_FASTQ_FNAME="test_sample_hg38_${SAMPLE_INDEX}.R1.fastq.gz"
R2_FASTQ_FNAME="test_sample_hg38_${SAMPLE_INDEX}.R2.fastq.gz"

find ~+ -type f \
    -name "*R1*${SAMPLE_INDEX}.fastq" \
    | sort | tr '\n' ' ' | xargs cat \
    | gzip -c \
    > ${OUT_DIR}/${R1_FASTQ_FNAME}

find ~+ -type f \
    -name "*R2*${SAMPLE_INDEX}.fastq" \
    | sort | tr '\n' ' ' | xargs cat \
    | gzip -c \
    > ${OUT_DIR}/${R2_FASTQ_FNAME} 


time2=$(date +%s)
echo "Job ended at "$(date)
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF
