#!/bin/bash
#
INPUT_TSV=$1
CHRX_LENGTH="156040895"
CHRY_LENGTH="57227415"
QUANTUM="1000"

set -euxo pipefail

rm -f uri_list.txt file_list.txt
while read LINE; do
    URI=$(echo ${LINE} | tr ' ' '\t' | cut -f 2)
    echo ${URI} >> uri_list.txt
    echo $(basename ${URI}) >> file_list.txt
done < ${INPUT_TSV}
cat uri_list.txt | gcloud storage cp -I .
N_FILES=$(wc -l < file_list.txt)
java PlotDipcallBeds ${CHRX_LENGTH} ${CHRY_LENGTH} ${QUANTUM} file_list.txt ${N_FILES}
