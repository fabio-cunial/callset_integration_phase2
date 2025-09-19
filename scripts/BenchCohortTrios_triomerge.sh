#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios_only_5.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_squish"
INPUT_DIR_50BP="/Users/fcunial/Downloads/BenchCohortTrios_squish/ad_denovo_plot/gt_matrix_only50bp"
INPUT_DIR_TRIOMERGE="/Users/fcunial/Downloads/BenchCohortTrios_squish/ad_denovo_plot/gt_matrix_triomerge"
INPUT_DIR_SAM_PARAMS_ALL="/Users/fcunial/Downloads/BenchCohortTrios_squish/ad_denovo_plot/gt_matrix_sam_all_calls"
INPUT_DIR_SAM_PARAMS_50BP="/Users/fcunial/Downloads/BenchCohortTrios_squish/ad_denovo_plot/gt_matrix_sam_50bp"

MAX_AD="50"

set -euxo pipefail


rm -f triomerge_all.csv triomerge_tr.csv triomerge_not_tr.csv
cut -f 2 ${TRIOS_TSV} > children.txt
while read CHILD_ID; do
    for SUFFIX in all tr not_tr; do
        # 1. Mendelian error
        # V1
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        ROW="${N_GOOD_ALT},${N_MERR}"
        # >=50bp cohort merge, kanpig.
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_50BP}/${CHILD_ID}_only_50bp_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_50BP}/${CHILD_ID}_only_50bp_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        # Trio merge, kanpig.
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_TRIOMERGE}/${CHILD_ID}_kanpig_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_TRIOMERGE}/${CHILD_ID}_kanpig_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        # Sam's params, all calls.
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_SAM_PARAMS_ALL}/${CHILD_ID}_sam_params_all_calls_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_SAM_PARAMS_ALL}/${CHILD_ID}_sam_params_all_calls_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        # Sam's params, >=50bp.
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_SAM_PARAMS_50BP}/${CHILD_ID}_only_50bp_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_SAM_PARAMS_50BP}/${CHILD_ID}_only_50bp_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        
        # 2. De novo rate
        # V1
        DENOVO=$(java CheckDeNovoOld ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}_gtmatrix.txt)
        ROW="${ROW},${DENOVO}"
        # >=50bp cohort merge, kanpig.
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_50BP}/${CHILD_ID}_only_50bp_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # Trio merge, kanpig.
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_TRIOMERGE}/${CHILD_ID}_kanpig_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # Sam's params, all calls.
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_SAM_PARAMS_ALL}/${CHILD_ID}_sam_params_all_calls_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # Sam's params, >=50bp.
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_SAM_PARAMS_50BP}/${CHILD_ID}_only_50bp_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        echo "${ROW}" >> triomerge_${SUFFIX}.csv
    done
done < children.txt
rm -f children.txt
