#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios_only_5.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_squish"
INPUT_DIR_PHAB_CALL_B50="/Users/fcunial/Downloads/phab_trios/construct_all_bench_50"
INPUT_DIR_PHAB_C50_B50="/Users/fcunial/Downloads/phab_trios/construct_50_bench_50"

MAX_AD="50"

set -euxo pipefail


rm -f triophab_all.csv triophab_tr.csv triophab_not_tr.csv
cut -f 2 ${TRIOS_TSV} > children.txt
while read CHILD_ID; do
    for SUFFIX in all tr not_tr; do
        # 1. Mendelian error
        # V1
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        ROW="${N_GOOD_ALT},${N_MERR}"
        # phab
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_PHAB_CALL_B50}/${CHILD_ID}_phab_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_PHAB_CALL_B50}/${CHILD_ID}_phab_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_PHAB_C50_B50}/${CHILD_ID}_phab_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_PHAB_C50_B50}/${CHILD_ID}_phab_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        
        # 2. De novo rate
        # V1
        DENOVO=$(java CheckDeNovoOld ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}_gtmatrix.txt)
        ROW="${ROW},${DENOVO}"
        # phab
        DENOVO=$(java CheckDeNovoOldPhab ${INPUT_DIR_PHAB_CALL_B50}/${CHILD_ID}_phab_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        DENOVO=$(java CheckDeNovoOldPhab ${INPUT_DIR_PHAB_C50_B50}/${CHILD_ID}_phab_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        echo "${ROW}" >> triophab_${SUFFIX}.csv
    done
done < children.txt
rm -f children.txt
