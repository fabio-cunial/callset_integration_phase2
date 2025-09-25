#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios_only_5.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_squish"
INPUT_DIR_DIPCALL="/Users/fcunial/Downloads/BenchCohortTrios_stratifications/dipcall"
INPUT_DIR_NIST="/Users/fcunial/Downloads/BenchCohortTrios_stratifications/nist"

MAX_AD="50"

set -euxo pipefail


rm -f stratifications_all.csv
cut -f 2 ${TRIOS_TSV} > children.txt
while read CHILD_ID; do
    # 1. Mendelian error
    # V1
    N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_all.txt | cut -f 2)
    N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_all.txt | cut -f 2)
    ROW="${N_GOOD_ALT},${N_MERR}"
    # DIPCALL
    N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_DIPCALL}/${CHILD_ID}_v1_tr.txt | cut -f 2)
    N_MERR=$(grep ^nmerr ${INPUT_DIR_DIPCALL}/${CHILD_ID}_v1_tr.txt | cut -f 2)
    ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
    # NIST
    N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_NIST}/${CHILD_ID}_v1_tr.txt | cut -f 2)
    N_MERR=$(grep ^nmerr ${INPUT_DIR_NIST}/${CHILD_ID}_v1_tr.txt | cut -f 2)
    ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
    
    # 2. De novo rate
    # V1
    DENOVO=$(java CheckDeNovoOld ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_all_gtmatrix.txt)
    ROW="${ROW},${DENOVO}"
    # DIPCALL
    DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_DIPCALL}/${CHILD_ID}_v1_tr_gtmatrix.txt ${MAX_AD})
    ROW="${ROW},${DENOVO}"
    DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_NIST}/${CHILD_ID}_v1_tr_gtmatrix.txt ${MAX_AD})
    ROW="${ROW},${DENOVO}"
    
    echo "${ROW}" >> stratifications_all.csv
done < children.txt
rm -f children.txt
