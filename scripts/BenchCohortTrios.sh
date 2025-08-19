#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios"

set -euxo pipefail


rm -f trios_all.csv trios_tr.csv trios_not_tr.csv
cut -f 2 ${TRIOS_TSV} > children.txt
while read CHILD_ID; do
    # All
    KANPIG_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_kanpig_all.txt | cut -f 2)
    KANPIG_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_kanpig_all.txt | cut -f 2)
    
    F07_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_07_all.txt | cut -f 2)
    F07_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_07_all.txt | cut -f 2)
    
    F09_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_09_all.txt | cut -f 2)
    F09_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_09_all.txt | cut -f 2)
    
    COHORT_MERGED_07_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_merged_07_all.txt | cut -f 2)
    COHORT_MERGED_07_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_merged_07_all.txt | cut -f 2)
    
    COHORT_MERGED_09_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_merged_09_all.txt | cut -f 2)
    COHORT_MERGED_09_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_merged_09_all.txt | cut -f 2)
    
    COHORT_REGENOTYPED_07_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_all.txt | cut -f 2)
    COHORT_REGENOTYPED_07_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_all.txt | cut -f 2)
    
    COHORT_REGENOTYPED_09_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_09_all.txt | cut -f 2)
    COHORT_REGENOTYPED_09_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_09_all.txt | cut -f 2)
    
    echo "${KANPIG_N_GOOD_ALT},${KANPIG_N_MERR},${F07_N_GOOD_ALT},${F07_N_MERR},${F09_N_GOOD_ALT},${F09_N_MERR},${COHORT_MERGED_07_N_GOOD_ALT},${COHORT_MERGED_07_N_MERR},${COHORT_MERGED_09_N_GOOD_ALT},${COHORT_MERGED_09_N_MERR},${COHORT_REGENOTYPED_07_N_GOOD_ALT},${COHORT_REGENOTYPED_07_N_MERR},${COHORT_REGENOTYPED_09_N_GOOD_ALT},${COHORT_REGENOTYPED_09_N_MERR}" >> trios_all.csv
    
    # Inside TRs
    KANPIG_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_kanpig_tr.txt | cut -f 2)
    KANPIG_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_kanpig_tr.txt | cut -f 2)
    
    F07_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_07_tr.txt | cut -f 2)
    F07_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_07_tr.txt | cut -f 2)
    
    F09_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_09_tr.txt | cut -f 2)
    F09_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_09_tr.txt | cut -f 2)
    
    COHORT_MERGED_07_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_merged_07_tr.txt | cut -f 2)
    COHORT_MERGED_07_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_merged_07_tr.txt | cut -f 2)
    
    COHORT_MERGED_09_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_merged_09_tr.txt | cut -f 2)
    COHORT_MERGED_09_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_merged_09_tr.txt | cut -f 2)
    
    COHORT_REGENOTYPED_07_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_tr.txt | cut -f 2)
    COHORT_REGENOTYPED_07_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_tr.txt | cut -f 2)
    
    COHORT_REGENOTYPED_09_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_09_tr.txt | cut -f 2)
    COHORT_REGENOTYPED_09_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_09_tr.txt | cut -f 2)
    
    echo "${KANPIG_N_GOOD_ALT},${KANPIG_N_MERR},${F07_N_GOOD_ALT},${F07_N_MERR},${F09_N_GOOD_ALT},${F09_N_MERR},${COHORT_MERGED_07_N_GOOD_ALT},${COHORT_MERGED_07_N_MERR},${COHORT_MERGED_09_N_GOOD_ALT},${COHORT_MERGED_09_N_MERR},${COHORT_REGENOTYPED_07_N_GOOD_ALT},${COHORT_REGENOTYPED_07_N_MERR},${COHORT_REGENOTYPED_09_N_GOOD_ALT},${COHORT_REGENOTYPED_09_N_MERR}" >> trios_tr.csv
    
    # Outside TRs
    KANPIG_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_kanpig_not_tr.txt | cut -f 2)
    KANPIG_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_kanpig_not_tr.txt | cut -f 2)
    
    F07_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_07_not_tr.txt | cut -f 2)
    F07_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_07_not_tr.txt | cut -f 2)
    
    F09_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_09_not_tr.txt | cut -f 2)
    F09_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_09_not_tr.txt | cut -f 2)
    
    COHORT_MERGED_07_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_merged_07_not_tr.txt | cut -f 2)
    COHORT_MERGED_07_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_merged_07_not_tr.txt | cut -f 2)
    
    COHORT_MERGED_09_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_merged_09_not_tr.txt | cut -f 2)
    COHORT_MERGED_09_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_merged_09_not_tr.txt | cut -f 2)
    
    COHORT_REGENOTYPED_07_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_not_tr.txt | cut -f 2)
    COHORT_REGENOTYPED_07_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_not_tr.txt | cut -f 2)
    
    COHORT_REGENOTYPED_09_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_09_not_tr.txt | cut -f 2)
    COHORT_REGENOTYPED_09_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_09_not_tr.txt | cut -f 2)
    
    echo "${KANPIG_N_GOOD_ALT},${KANPIG_N_MERR},${F07_N_GOOD_ALT},${F07_N_MERR},${F09_N_GOOD_ALT},${F09_N_MERR},${COHORT_MERGED_07_N_GOOD_ALT},${COHORT_MERGED_07_N_MERR},${COHORT_MERGED_09_N_GOOD_ALT},${COHORT_MERGED_09_N_MERR},${COHORT_REGENOTYPED_07_N_GOOD_ALT},${COHORT_REGENOTYPED_07_N_MERR},${COHORT_REGENOTYPED_09_N_GOOD_ALT},${COHORT_REGENOTYPED_09_N_MERR}" >> trios_not_tr.csv    
done < children.txt
rm -f children.txt
