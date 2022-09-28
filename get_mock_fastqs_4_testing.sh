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

## set params #####################################

DATA_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/data"
OUT_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/get_mock_fastqs_4_testing/outputs"
TMP_DIR="${OUT_DIR}/TMP"
mkdir -p $TMP_DIR

FASTQ_FLIST_PATH="/cluster/home/t110409uhn/git/make_toy_fastqs/backups/fastq_flist.bak"
BEDPE_FLIST_PATH="/cluster/home/t110409uhn/git/make_toy_fastqs/backups/bedpe_flist.bak"

## chromosome to sample
CHR="chr21"

## lines in BEDPE, so 400K*2reads*3samples=2.5m reads
NUM_LINES=400000

## seed for repeatable shuffle
SHUF_SEED=42


###################################################

echo "Pivot fastq files."
## pivot wide FASTQs to get Phred

FASTQ_FLIST=$(cat ${FASTQ_FLIST_PATH} | tr "\n" " ")
SAMPLES_FASTQ_ARR=($FASTQ_FLIST)

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


echo "Grab Chr Start End FragName Strand columns from BEDPE file."
## grab chr start end from bedpe.gz

BEDPE_FLIST=$(cat ${BEDPE_FLIST_PATH} | tr "\n" " ")
SAMPLES_BEDPE_ARR=($BEDPE_FLIST)

#for SAMPLE_FNAME in "${SAMPLES_BEDPE_ARR[@]:0:1}" ; do
for SAMPLE_FNAME in "${SAMPLES_BEDPE_ARR[@]}" ; do
    echo "Processing ${SAMPLE_FNAME}"
    
    SAMPLE_CHR_SHUF_BED_FNAME="${SAMPLE_FNAME%.bedpe.gz*}.${CHR}.shuf${SHUF_SEED}.bedpe"
    zcat ${DATA_DIR}/${SAMPLE_FNAME} \
        | awk 'BEGIN {OFS="\t"} ($1 == "chr21" && $4 == "chr21") {print $1,$2,$3,$4,$5,$6,"@"$7,$10,$11}' \
        | shuf -n ${NUM_LINES} --random-source=<(yes ${SHUF_SEED}) \
        > ${TMP_DIR}/${SAMPLE_CHR_SHUF_BED_FNAME}

    SAMPLE_CHR_SHUF_BED_R1_FNAME="${SAMPLE_CHR_SHUF_BED_FNAME%.*}.R1.bed.sortd"
    SAMPLE_CHR_SHUF_BED_R2_FNAME="${SAMPLE_CHR_SHUF_BED_FNAME%.*}.R2.bed.sortd"
    cat ${TMP_DIR}/${SAMPLE_CHR_SHUF_BED_FNAME} \
        | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$7,$8}' \
        | sed 's/|/\t/' | sed 's/\./\t/' \
        | sort -k4,4 \
        > ${TMP_DIR}/${SAMPLE_CHR_SHUF_BED_R1_FNAME}
    cat ${TMP_DIR}/${SAMPLE_CHR_SHUF_BED_FNAME} \
        | awk 'BEGIN {OFS="\t"} {print $4,$5,$6,$7,$9}' \
        | sed 's/|/\t/' | sed 's/\./\t/' \
        | sort -k4,4 \
        > ${TMP_DIR}/${SAMPLE_CHR_SHUF_BED_R2_FNAME}
done


echo "Joining pivoted sorted FASTQ and sorted BED, then convert to fastqs."
## join pivoted FASTQ and chr star end
## and then
## Rscript to convert .joined to fastqs
module load R/4.1
## need packages: tidyverse, ggplot2, GenomicFeatures, BSgenome.Hsapiens.UCSC.hg38

RSCRIPT_PATH="/cluster/home/t110409uhn/git/make_toy_fastqs/shuffled_joined_to_fastq_v3.R"
SAMPLE_INDEX="ACTGACTG:ACTGACTG"

#for SAMPLE_FNAME in "${SAMPLES_BEDPE_ARR[@]:0:1}" ; do
for SAMPLE_FNAME in "${SAMPLES_BEDPE_ARR[@]}" ; do
    echo "Processing ${SAMPLE_FNAME%.bedpe.gz*}"

    SAMPLE_CHR_SHUF_BED_R1_FNAME="${SAMPLE_FNAME%.bedpe.gz*}.${CHR}.shuf${SHUF_SEED}.R1.bed.sortd"
    SAMPLE_CHR_SHUF_BED_R2_FNAME="${SAMPLE_FNAME%.bedpe.gz*}.${CHR}.shuf${SHUF_SEED}.R2.bed.sortd"

    SAMPLE_FASTQ_R_PIVOT_R1_FNAME="${SAMPLE_FNAME%.bedpe.gz*}_lib1_R1.fastq.pivot.sortd"
    SAMPLE_FASTQ_R_PIVOT_R2_FNAME="${SAMPLE_FNAME%.bedpe.gz*}_lib1_R2.fastq.pivot.sortd"

    SAMPLE_R1_JOINED="${SAMPLE_FNAME%.bedpe.gz*}.${CHR}.shuf${SHUF_SEED}.R1.joined"
    SAMPLE_R2_JOINED="${SAMPLE_FNAME%.bedpe.gz*}.${CHR}.shuf${SHUF_SEED}.R2.joined"

    join -t $'\t' -1 1 -2 4 \
        ${TMP_DIR}/${SAMPLE_FASTQ_R_PIVOT_R1_FNAME} \
        ${TMP_DIR}/${SAMPLE_CHR_SHUF_BED_R1_FNAME} \
        > ${TMP_DIR}/${SAMPLE_R1_JOINED}
    join -t $'\t' -1 1 -2 4 \
        ${TMP_DIR}/${SAMPLE_FASTQ_R_PIVOT_R2_FNAME} \
        ${TMP_DIR}/${SAMPLE_CHR_SHUF_BED_R2_FNAME} \
        > ${TMP_DIR}/${SAMPLE_R2_JOINED}

    Rscript ${RSCRIPT_PATH} \
        --R1_fpath ${TMP_DIR}/${SAMPLE_R1_JOINED} \
        --R2_fpath ${TMP_DIR}/${SAMPLE_R2_JOINED} \
        --output_dir ${TMP_DIR} \
        --sample_index ${SAMPLE_INDEX}
done


echo "Concatenating fastqs from all samples subsampled, for R1, then R2."
## cat all fastqs for R1, then R2

SAMPLE_INDEX2="${SAMPLE_INDEX//:/}"
R1_FASTQ_FNAME="test_sample_shuf${SHUF_SEED}_hg38_${SAMPLE_INDEX2}.R1.fastq.gz"
R2_FASTQ_FNAME="test_sample_shuf${SHUF_SEED}_hg38_${SAMPLE_INDEX2}.R2.fastq.gz"

find ${TMP_DIR} ~+ -type f \
    -name "*R1*${SAMPLE_INDEX2}.fastq" \
    | sort | tr '\n' ' ' | xargs cat \
    | gzip -c \
    > ${OUT_DIR}/${R1_FASTQ_FNAME}

find ${TMP_DIR} ~+ -type f \
    -name "*R2*${SAMPLE_INDEX2}.fastq" \
    | sort | tr '\n' ' ' | xargs cat \
    | gzip -c \
    > ${OUT_DIR}/${R2_FASTQ_FNAME} 


time2=$(date +%s)
echo "Job ended at "$(date)
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF
