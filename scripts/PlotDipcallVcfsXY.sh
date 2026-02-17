#!/bin/bash
#
INPUT_TSV=$1
MIN_SV_LENGTH=$2
MAX_SV_LENGTH=$3
CHROMOSOME=$4

N_THREADS="8"

set -euxo pipefail


# Downloading all dipcall VCFs and BEDs and assigning them standard names
rm -f id_list.txt
while read LINE; do
    SAMPLE_ID=$(echo ${LINE} | tr ' ' '\t' | cut -f 1)
    echo ${SAMPLE_ID} >> id_list.txt
    URI=$(echo ${LINE} | tr ' ' '\t' | cut -f 2)
    gcloud storage cp ${URI} ./${SAMPLE_ID}.vcf.gz
    URI=$(echo ${LINE} | tr ' ' '\t' | cut -f 3)
    gcloud storage cp ${URI} ./${SAMPLE_ID}.bed
done < ${INPUT_TSV}

# Counting in parallel
rm -f *.csv
xargs --arg-file=id_list.txt --max-lines=1 --max-procs=${N_THREADS} ./PlotDipcallVcfsXY_Impl.sh ${MIN_SV_LENGTH} ${MAX_SV_LENGTH} ${CHROMOSOME}
cat *.csv > matrix.csv
