#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios_only_5.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_filtered_truvari"

MAX_AD="50"

set -euxo pipefail


rm -f filtered_truvari_all.csv filtered_truvari_tr.csv filtered_truvari_not_tr.csv
cut -f 2 ${TRIOS_TSV} > children.txt
while read CHILD_ID; do
    for SUFFIX in all tr not_tr; do
        # 1. Mendelian error
        ## V1 original
        #N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        #N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        #ROW="${N_GOOD_ALT},${N_MERR}"
        
        # V1 latest run
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_v1_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_v1_${SUFFIX}.txt | cut -f 2)
        ROW="${N_GOOD_ALT},${N_MERR}"
        
        # 3_25_3
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_3_25_3_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_3_25_3_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        # 10_25_3
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_10_25_3_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_10_25_3_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        
        # 2. De novo rate
        ## V1 original
        #DENOVO=$(java CheckDeNovoOld ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}_gtmatrix.txt)
        #ROW="${ROW},${DENOVO}"
        # V1 latest run
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR}/${CHILD_ID}_v1_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # 3_25_3
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR}/${CHILD_ID}_3_25_3_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # 10_25_3
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR}/${CHILD_ID}_10_25_3_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"

        echo "${ROW}" >> filtered_truvari_${SUFFIX}.csv
    done
done < children.txt
rm -f children.txt
