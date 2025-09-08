#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios_only_5.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_squish"
INPUT_DIR_NEW="/Users/fcunial/Downloads/BenchCohortTrios_squish/fixed_gts_2_prime"

IDS="fixed_gts_2_prime"

set -euxo pipefail


rm -f squish_trios_all.csv squish_trios_tr.csv squish_trios_not_tr.csv
cut -f 2 ${TRIOS_TSV} > children.txt
while read CHILD_ID; do
    for SUFFIX in all tr not_tr; do
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        ROW="${N_GOOD_ALT},${N_MERR}"
        for ID in ${IDS}; do
            N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_NEW}/${CHILD_ID}_${ID}_${SUFFIX}.txt | cut -f 2)
            N_MERR=$(grep ^nmerr ${INPUT_DIR_NEW}/${CHILD_ID}_${ID}_${SUFFIX}.txt | cut -f 2)
            ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        done
        DENOVO=$(java CheckDeNovoOld ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}_gtmatrix.txt)
        ROW="${ROW},${DENOVO}"
        for ID in ${IDS}; do
            DENOVO=$(java CheckDeNovo ${INPUT_DIR_NEW}/${CHILD_ID}_${ID}_${SUFFIX}_gtmatrix.txt 1)
            ROW="${ROW},${DENOVO}"
        done
        echo "${ROW}" >> squish_trios_${SUFFIX}.csv
    done
done < children.txt
rm -f children.txt
