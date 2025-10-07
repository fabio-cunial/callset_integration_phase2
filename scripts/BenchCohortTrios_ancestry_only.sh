#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/BenchCohortTrios_ancestry_only/trios_5_afr.tsv"
TRIOS_TSV_V1="/Users/fcunial/Downloads/trios_only_5.tsv"

INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_ancestry_only"
INPUT_DIR_V1="/Users/fcunial/Downloads/BenchCohortTrios_filtered_truvari2"

MAX_AD="50"

set -euxo pipefail


rm -f ancestry_only_all.csv ancestry_only_tr.csv ancestry_only_not_tr.csv
cut -f 2 ${TRIOS_TSV_V1} > children_v1.txt
cut -f 2 ${TRIOS_TSV} > children_ancestry_only.txt
paste -d , children_v1.txt children_ancestry_only.txt > children.tsv
while read CHILDREN_ROW; do
    for SUFFIX in all tr not_tr; do
        # 1. Mendelian error
        ## V1 original
        #N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        #N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        #ROW="${N_GOOD_ALT},${N_MERR}"
        
        CHILD_ID_V1=$(echo ${CHILDREN_ROW} | cut -d , -f 1)
        CHILD_ID=$(echo ${CHILDREN_ROW} | cut -d , -f 2)
        
        # V1 latest run
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_V1}/${CHILD_ID_V1}_v1_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_V1}/${CHILD_ID_V1}_v1_${SUFFIX}.txt | cut -f 2)
        ROW="${N_GOOD_ALT},${N_MERR}"
        
        # ancestry_only
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_afr_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_afr_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        
        # 2. De novo rate
        ## V1 original
        #DENOVO=$(java CheckDeNovoOld ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}_gtmatrix.txt)
        #ROW="${ROW},${DENOVO}"
        # V1 latest run
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_V1}/${CHILD_ID_V1}_v1_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # ancestry_only
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR}/${CHILD_ID}_afr_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"

        echo "${ROW}" >> ancestry_only_${SUFFIX}.csv
    done
done < children.tsv
rm -f children.tsv
