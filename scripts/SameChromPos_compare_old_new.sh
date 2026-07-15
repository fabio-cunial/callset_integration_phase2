#!/bin/bash
#
REMOTE_INDIR_OLD="gs://fc-7f861a33-ddb4-4b2f-8d10-5679c9df6108/v3/15x/workpackage_3_second_attempt"
REMOTE_INDIR_NEW="gs://fc-7f861a33-ddb4-4b2f-8d10-5679c9df6108/v3/15x/testing_bcftools_annotate_fix_workpackage_3"
SAMPLE_IDS="HG00097 HG00099 HG00126 HG00128 HG00133 HG00140 HG00146 HG002 HG00232 HG00235"

set -euxo pipefail

N_THREADS="8"



# --------------------------------- Functions ----------------------------------

cat << 'END' > download_chunks.sh
#!/bin/bash
REMOTE_FILE=$1
LOCAL_FILE=$2
gcloud storage cp ${REMOTE_FILE} ${LOCAL_FILE}
END
chmod +x download_chunks.sh


function Download() {
    local SAMPLE_ID=$1
    local REMOTE_INDIR=$2
    local SUFFIX=$3

    # Recomposing chunks
    rm -f list.wsv list.txt
    for CHUNK_ID in $(seq 0 115); do
        echo "${REMOTE_INDIR}/chunk_${CHUNK_ID}/${SAMPLE_ID}.bcf ${SAMPLE_ID}_${CHUNK_ID}.bcf" >> list.wsv
        echo "${SAMPLE_ID}_${CHUNK_ID}.bcf" >> list.txt
    done
    xargs --arg-file=list.wsv --max-lines=1 --max-procs=${N_THREADS} ./download_chunks.sh
    bcftools concat --naive --file-list list.txt --output-type b --output ${SAMPLE_ID}.bcf
    rm -f ${SAMPLE_ID}_*.bcf*

    # Preparing for comparisons
    bcftools view --include 'ABS(SVLEN)<50' --output-type b ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}_${SUFFIX}_lt_50.bcf
    bcftools view --include 'ABS(SVLEN)>=50' --output-type b ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}_${SUFFIX}_ge_50.bcf
    rm -f ${SAMPLE_ID}.bcf*
    bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\n' ${SAMPLE_ID}_${SUFFIX}_lt_50.bcf | sort > ${SAMPLE_ID}_${SUFFIX}_lt_50.txt
    bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\n' ${SAMPLE_ID}_${SUFFIX}_ge_50.bcf | sort > ${SAMPLE_ID}_${SUFFIX}_ge_50.txt
}


function Compare() {
    local SAMPLE_ID=$1
    local MODE=$2
    local OLD_TXT=$3
    local NEW_TXT=$4
    local OUTPUT_CSV=$5

    comm -2 -3 ${OLD_TXT} ${NEW_TXT} > ${SAMPLE_ID}_${MODE}_old_only.txt
    comm -1 -3 ${OLD_TXT} ${NEW_TXT} > ${SAMPLE_ID}_${MODE}_new_only.txt
    comm -1 -2 ${OLD_TXT} ${NEW_TXT} > ${SAMPLE_ID}_${MODE}_common.txt
    OLD_ONLY=$(wc -l < ${SAMPLE_ID}_${MODE}_old_only.txt)
    NEW_ONLY=$(wc -l < ${SAMPLE_ID}_${MODE}_new_only.txt)
    COMMON=$(wc -l < ${SAMPLE_ID}_${MODE}_common.txt)
    echo "${OLD_ONLY},${NEW_ONLY},${COMMON}" >> ${OUTPUT_CSV}
}




# ------------------------------- Main program ---------------------------------

rm -f results_lt_50.csv results_ge_50.csv
for SAMPLE_ID in ${SAMPLE_IDS}; do
    Download ${SAMPLE_ID} ${REMOTE_INDIR_OLD} old
    Download ${SAMPLE_ID} ${REMOTE_INDIR_NEW} new
    Compare ${SAMPLE_ID} lt_50 ${SAMPLE_ID}_old_lt_50.txt ${SAMPLE_ID}_new_lt_50.txt results_lt_50.csv
    Compare ${SAMPLE_ID} ge_50 ${SAMPLE_ID}_old_ge_50.txt ${SAMPLE_ID}_new_ge_50.txt results_ge_50.csv
done
