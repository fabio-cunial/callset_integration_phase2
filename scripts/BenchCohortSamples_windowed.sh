#!/bin/bash
#
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortSamples_windowed"
MAX_WINDOW="98";
SAMPLE="HG002"

set -euxo pipefail

for STEP in kanpig 07 09 cohort_merged_07 cohort_merged_09 cohort_regenotyped_07 cohort_regenotyped_09; do
    rm -f ${STEP}.csv
    for i in $(seq 0 ${MAX_WINDOW}); do
        P=$(grep precision ${INPUT_DIR}/${SAMPLE}_${STEP}_${i}.txt | cut -w -f 3)
        R=$(grep recall ${INPUT_DIR}/${SAMPLE}_${STEP}_${i}.txt | cut -w -f 3)
        echo "${P}${R}" >> ${STEP}.csv
    done
done
