#!/bin/bash
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH --account=cgptip
#SBATCH -t 1-00:00:00
#SBATCH --mem=32G
#SBATCH -J get_mock_fastqs_4_testing
#SBATCH -p himem
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -o ./logs_slurm/%j-%x.out
#SBATCH -e ./logs_slurm/%j-%x.err

echo "Job started at "$(date) 
time1=$(date +%s)

module load R/4.1
## need packages: tidyverse, ggplot2, GenomicFeatures, BSgenome.Hsapiens.UCSC.hg38


# set params for debuggin ###################################

#DATA_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/data"
#OUT_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/get_mock_fastqs_4_testing/output_allchrs"
#FASTQ_FLIST_PATH="./backups/fastq_flist.bak"
#BEDPE_FLIST_PATH="./backups/bedpe_flist.bak"
##CHR="all"
#CHR="chr21"
#NUM_LINES=400000
#SHUF_SEED=42
#SAMPLE_INDEX="ACTGACTG:ACTGACTG"

# getopts ###################################################
usage(){
    echo 
    echo "Usage: bash/sbatch get_mock_fastqs_4_testing.sh -d [DATA_DIR] -f [FASTQ_FLIST] -b [BEDPE_FLIST] -m [SAMPLE_NAME] -i [SAMPLE_INDEX] -n [NUM_LINES] -s [SHUFFLE_SEED] -c [CHROMOSOME] -o [OUT_DIR] [-t [KEEP_TMP]]"
    echo 
}
no_args="true"
KEEP_TMP="false"

## Help 
Help()
{
    # Display Help
    echo 
    echo "Makes toy fastqs subsampled from several samples and with hg38 sequences to avoid privacy issues."
    echo
    echo "Usage: bash/sbatch get_mock_fastqs_4_testing.sh -d [DATA_DIR] -f [FASTQ_FLIST] -b [BEDPE_FLIST] -m [SAMPLE_NAME] -i [SAMPLE_INDEX] -n [NUM_LINES] -s [SHUFFLE_SEED] -c [CHROMOSOME] -o [OUT_DIR] [-t [KEEP_TMP]]"
    echo "options:"
    echo "-h   [HELP]      Print help"
    echo "-d   [REQUIRED]  Data directory with original fastq and bedpe.gz files." 
    echo "-f   [REQUIRED]  FASTQ file list, containing names of fastqs to be subsampled."
    echo "-b   [REQUIRED]  BEDPE file list, containing names of bedpe.gz to be subsampled."
    echo "-m   [REQUIRED]  Output sample file name."
    echo "-i   [REQUIRED]  Sample index for new toy fastq (e.g. ACTGACTG:ACTGACTG)."
    echo "-n   [REQUIRED]  Number of lines in bedpe.gz to subsample (e.g. 400000)"
    echo "-s   [REQUIRED]  Seed for consistent subsampling (e.g. 42)"
    echo "-c   [REQUIRED]  Specific chromosome to subsample, or all to subsample from all chromosomes [e.g. chr21, or all]"
    echo "-o   [REQUIRED]  Output directory (full path)"
    echo "-t   [OPTIONAL]  Keep TMP directory?"
    echo
}

## Get the options
while getopts ":hd:f:b:m:i:n:s:c:o:t" option; do
    case "${option}" in
        h) Help
           exit;;
        d) DATA_DIR=${OPTARG};;
        f) FASTQ_FLIST_PATH=${OPTARG};;
        b) BEDPE_FLIST_PATH=${OPTARG};;
        m) OUT_FNAME=${OPTARG};;
        i) SAMPLE_INDEX=${OPTARG};;
        n) NUM_LINES=${OPTARG};;
        s) SHUF_SEED=${OPTARG};;
        c) CHR=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
        t) KEEP_TMP="true";;
       \?) echo "Error: Invalid option"
           exit;;
    esac
    no_args="false"
done

[[ "$no_args" == "true" ]] && { usage; exit 1; }

echo "Running script... "
echo "Number of lines to subsample from each bedpe.gz: $NUM_LINES"
echo "Chromosome to subsample from:                    $CHR"
echo "Seed set for subsampling:                        $SHUF_SEED"
echo "New toy sample index is:                         $SAMPLE_INDEX"
echo "Output path:                                     $OUT_DIR"
echo "Output filename:                                 $OUT_FNAME"

CHR_LIST_ALL="./chr_list_wBAC.txt"
CHR_LIST_BACs="./chr_list_BACs.txt"
RSCRIPT_PATH="./shuffled_joined_to_fastq_wBAC.R"

TMP_DIR="${OUT_DIR}/TMP"
mkdir -p $TMP_DIR

# Main program ###############################################

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

    if [ ${CHR} == "all" ]; then
        echo "Subsampling ${CHR} chromosomes..."
        SAMPLE_CHR_SHUF_BED_FNAME="${SAMPLE_FNAME%.bedpe.gz*}.${CHR}.shuf${SHUF_SEED}.bedpe"
        zcat ${DATA_DIR}/${SAMPLE_FNAME} \
            | awk 'BEGIN {OFS="\t"} ($1 == $4) {print $1,$2,$3,$4,$5,$6,"@"$7,$10,$11}' \
            | awk 'FNR==NR {a[$1]; next}; $1 in a' ${CHR_LIST_ALL} - \
            | shuf -n ${NUM_LINES} --random-source=<(yes ${SHUF_SEED}) \
            > ${TMP_DIR}/${SAMPLE_CHR_SHUF_BED_FNAME}
    elif [ ${CHR} == "BACs" ]; then
        echo "Subsampling ${CHR} chromosomes..."
        SAMPLE_CHR_SHUF_BED_FNAME="${SAMPLE_FNAME%.bedpe.gz*}.${CHR}.shuf${SHUF_SEED}.bedpe"
        zcat ${DATA_DIR}/${SAMPLE_FNAME} \
            | awk 'BEGIN {OFS="\t"} ($1 == $4) {print $1,$2,$3,$4,$5,$6,"@"$7,$10,$11}' \
            | awk 'FNR==NR {a[$1]; next}; $1 in a' ${CHR_LIST_BACs} - \
            > ${TMP_DIR}/${SAMPLE_CHR_SHUF_BED_FNAME}
            #| shuf -n ${NUM_LINES} --random-source=<(yes ${SHUF_SEED}) \
    else
        echo "Subsampling ${CHR} only..."
        SAMPLE_CHR_SHUF_BED_FNAME="${SAMPLE_FNAME%.bedpe.gz*}.${CHR}.shuf${SHUF_SEED}.bedpe"
        zcat ${DATA_DIR}/${SAMPLE_FNAME} \
            | awk -v CHR="${CHR}" 'BEGIN {OFS="\t"} ($1 == CHR && $4 == CHR) {print $1,$2,$3,$4,$5,$6,"@"$7,$10,$11}' \
            | shuf -n ${NUM_LINES} --random-source=<(yes ${SHUF_SEED}) \
            > ${TMP_DIR}/${SAMPLE_CHR_SHUF_BED_FNAME}
    fi
    
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
R1_FASTQ_FNAME="${OUT_FNAME}_${CHR}_shuf${SHUF_SEED}_F19K16_F24B22_${SAMPLE_INDEX2}.R1.fastq.gz"
R2_FASTQ_FNAME="${OUT_FNAME}_${CHR}_shuf${SHUF_SEED}_F19K16_F24B22_${SAMPLE_INDEX2}.R2.fastq.gz"

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


if "${KEEP_TMP}"; then
    echo "Keeping TMP directory."
else
    rm -r ${TMP_DIR}
    echo "Removing TMP directory."
fi


time2=$(date +%s)
echo "Job ended at "$(date)
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF
