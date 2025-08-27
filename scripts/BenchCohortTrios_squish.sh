#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios_only_5.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_squish"

set -euxo pipefail


rm -f squish_trios_all.csv squish_trios_tr.csv squish_trios_not_tr.csv
cut -f 2 ${TRIOS_TSV} > children.txt
while read CHILD_ID; do
    for SUFFIX in all tr not_tr; do
        COHORT_REGENOTYPED_07_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        COHORT_REGENOTYPED_07_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
    
        SQUISH_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_squish_${SUFFIX}.txt | cut -f 2)
        SQUISH_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_squish_${SUFFIX}.txt | cut -f 2)
    
        AB_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_ab_${SUFFIX}.txt | cut -f 2)
        AB_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_ab_${SUFFIX}.txt | cut -f 2)
    
        MAXHOM_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_maxhom_${SUFFIX}.txt | cut -f 2)
        MAXHOM_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_maxhom_${SUFFIX}.txt | cut -f 2)
    
        FPENALTY1_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_fpenalty1_${SUFFIX}.txt | cut -f 2)
        FPENALTY1_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_fpenalty1_${SUFFIX}.txt | cut -f 2)
        FPENALTY2_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_fpenalty2_${SUFFIX}.txt | cut -f 2)
        FPENALTY2_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_fpenalty2_${SUFFIX}.txt | cut -f 2)
    
        GPENALTY_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_gpenalty_${SUFFIX}.txt | cut -f 2)
        GPENALTY_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_gpenalty_${SUFFIX}.txt | cut -f 2)
    
        FNMAX1_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_fnmax1_${SUFFIX}.txt | cut -f 2)
        FNMAX1_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_fnmax1_${SUFFIX}.txt | cut -f 2)
        FNMAX2_N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_fnmax2_${SUFFIX}.txt | cut -f 2)
        FNMAX2_N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_fnmax2_${SUFFIX}.txt | cut -f 2)
    
        COHORT_REGENOTYPED_07_DENOVO=$(java CheckDeNovo ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}_gtmatrix.txt)
        SQUISH_DENOVO=$(java CheckDeNovo ${INPUT_DIR}/${CHILD_ID}_squish_${SUFFIX}_gtmatrix.txt)
        AB_DENOVO=$(java CheckDeNovo ${INPUT_DIR}/${CHILD_ID}_ab_${SUFFIX}_gtmatrix.txt)
        MAXHOM_DENOVO=$(java CheckDeNovo ${INPUT_DIR}/${CHILD_ID}_maxhom_${SUFFIX}_gtmatrix.txt)
        FPENALTY1_DENOVO=$(java CheckDeNovo ${INPUT_DIR}/${CHILD_ID}_fpenalty1_${SUFFIX}_gtmatrix.txt)
        FPENALTY2_DENOVO=$(java CheckDeNovo ${INPUT_DIR}/${CHILD_ID}_fpenalty2_${SUFFIX}_gtmatrix.txt)
        GPENALTY_DENOVO=$(java CheckDeNovo ${INPUT_DIR}/${CHILD_ID}_gpenalty_${SUFFIX}_gtmatrix.txt)
        FNMAX1_DENOVO=$(java CheckDeNovo ${INPUT_DIR}/${CHILD_ID}_fnmax1_${SUFFIX}_gtmatrix.txt)
        FNMAX2_DENOVO=$(java CheckDeNovo ${INPUT_DIR}/${CHILD_ID}_fnmax2_${SUFFIX}_gtmatrix.txt)
    
        echo "${COHORT_REGENOTYPED_07_N_GOOD_ALT},${COHORT_REGENOTYPED_07_N_MERR},${SQUISH_N_GOOD_ALT},${SQUISH_N_MERR},${AB_N_GOOD_ALT},${AB_N_MERR},${MAXHOM_N_GOOD_ALT},${MAXHOM_N_MERR},${FPENALTY1_N_GOOD_ALT},${FPENALTY1_N_MERR},${FPENALTY2_N_GOOD_ALT},${FPENALTY2_N_MERR},${GPENALTY_N_GOOD_ALT},${GPENALTY_N_MERR},${FNMAX1_N_GOOD_ALT},${FNMAX1_N_MERR},${FNMAX2_N_GOOD_ALT},${FNMAX2_N_MERR},${COHORT_REGENOTYPED_07_DENOVO},${SQUISH_DENOVO},${AB_DENOVO},${MAXHOM_DENOVO},${FPENALTY1_DENOVO},${FPENALTY2_DENOVO},${GPENALTY_DENOVO},${FNMAX1_DENOVO},${FNMAX2_DENOVO}" >> squish_trios_${SUFFIX}.csv
    done
done < children.txt
rm -f children.txt
