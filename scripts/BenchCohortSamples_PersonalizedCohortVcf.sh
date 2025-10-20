#!/bin/bash
#
SAMPLES="HG03579 NA18906 NA19240 NA20129 NA21309"
#SAMPLES="HG002 HG00438 HG005 HG00621 HG00673 HG00733 HG00735 HG00741 HG01071 HG01106 HG01109 HG01123 HG01175 HG01243 HG01258 HG01358 HG01361 HG01891 HG01928 HG01952 HG01978 HG02055 HG02080 HG02109 HG02145 HG02148 HG02257 HG02486 HG02559 HG02572 HG02622 HG02630 HG02717 HG02723 HG02818 HG02886 HG03098 HG03453 HG03486 HG03492 HG03516 HG03540 HG03579 NA18906 NA19240 NA20129 NA21309"
MIN_N_SAMPLES="2 3 4 8 16 32 64 128 256 512 1024 2048"

set -euxo pipefail


INPUT_DIR="/Users/fcunial/Downloads/BenchCohortSamples_PersonalizedCohortVcf/kanpig_1_1_0"
for REGION in all tr not_tr; do
    rm -f ${INPUT_DIR}/personalized_${REGION}.csv
    for SAMPLE in ${SAMPLES}; do
        P=$(grep precision ${INPUT_DIR}/${SAMPLE}_v1_07_${REGION}.txt | cut -w -f 3)
        R=$(grep recall ${INPUT_DIR}/${SAMPLE}_v1_07_${REGION}.txt | cut -w -f 3)
        ROW="${P},${R}"
        for MIN in ${MIN_N_SAMPLES}; do
            P=$(grep precision ${INPUT_DIR}/${SAMPLE}_${MIN}_${REGION}.txt | cut -w -f 3)
            R=$(grep recall ${INPUT_DIR}/${SAMPLE}_${MIN}_${REGION}.txt | cut -w -f 3)
            ROW="${ROW},${P},${R}"
        done
        echo ${ROW} >> ${INPUT_DIR}/personalized_${REGION}.csv
    done
done

INPUT_DIR="/Users/fcunial/Downloads/BenchCohortSamples_PersonalizedCohortVcf/kanpig_1_0_2"
for REGION in all tr not_tr; do
    rm -f ${INPUT_DIR}/personalized_${REGION}.csv
    for SAMPLE in ${SAMPLES}; do
        P=$(grep precision ${INPUT_DIR}/${SAMPLE}_v1_07_${REGION}.txt | cut -w -f 3)
        R=$(grep recall ${INPUT_DIR}/${SAMPLE}_v1_07_${REGION}.txt | cut -w -f 3)
        ROW="${P},${R}"
        for MIN in ${MIN_N_SAMPLES}; do
            P=$(grep precision ${INPUT_DIR}/${SAMPLE}_${MIN}_${REGION}.txt | cut -w -f 3)
            R=$(grep recall ${INPUT_DIR}/${SAMPLE}_${MIN}_${REGION}.txt | cut -w -f 3)
            ROW="${ROW},${P},${R}"
        done
        echo ${ROW} >> ${INPUT_DIR}/personalized_${REGION}.csv
    done
done
