#!/bin/bash
#
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortSamples_vcfdist"
SAMPLES="HG002"
#SAMPLES="HG002 HG00438 HG005 HG00621 HG00673 HG00733 HG00735 HG00741 HG01071 HG01106 HG01109 HG01123 HG01175 HG01243 HG01258 HG01358 HG01361 HG01891 HG01928 HG01952 HG01978 HG02055 HG02080 HG02109 HG02145 HG02148 HG02257 HG02486 HG02559 HG02572 HG02622 HG02630 HG02717 HG02723 HG02818 HG02886 HG03098 HG03453 HG03486 HG03492 HG03516 HG03540 HG03579 NA18906 NA19240 NA20129 NA21309"

set -euxo pipefail

rm -f *.csv
for SAMPLE in ${SAMPLES}; do
    for STEP in kanpig 07 09 cohort_merged_07 cohort_merged_09 cohort_regenotyped_07 cohort_regenotyped_09; do
        for REGION in all tr not_tr; do
            P=$(grep ^SV ${INPUT_DIR}/${SAMPLE}_${STEP}_${REGION}.txt | grep BEST | cut -f 8)
            R=$(grep ^SV ${INPUT_DIR}/${SAMPLE}_${STEP}_${REGION}.txt | grep BEST | cut -f 9)
            echo "${P},${R}" >> ${STEP}_${REGION}.csv
        done
    done
done
