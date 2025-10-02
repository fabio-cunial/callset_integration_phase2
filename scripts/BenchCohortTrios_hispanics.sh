#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios_only_5.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_squish"
INPUT_DIR_HISPANICS="/Users/fcunial/Downloads/BenchCohortTrios_hardfilters/fraction_10"

MAX_AD="50"

set -euxo pipefail


rm -f triohispanics_all.csv triohispanics_tr.csv triohispanics_not_tr.csv
cut -f 2 ${TRIOS_TSV} > children.txt
while read CHILD_ID; do
    for SUFFIX in all tr not_tr; do
        # 1. Mendelian error
        ## V1 original
        #N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        #N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}.txt | cut -f 2)
        #ROW="${N_GOOD_ALT},${N_MERR}"
        
        # V1 latest run
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_HISPANICS}/${CHILD_ID}_v1_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_HISPANICS}/${CHILD_ID}_v1_${SUFFIX}.txt | cut -f 2)
        ROW="${N_GOOD_ALT},${N_MERR}"
        
        # minkfreq_4
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_HISPANICS}/${CHILD_ID}_minkfreq_4_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_HISPANICS}/${CHILD_ID}_minkfreq_4_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        # minkfreq_6
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_HISPANICS}/${CHILD_ID}_minkfreq_6_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_HISPANICS}/${CHILD_ID}_minkfreq_6_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        # pileupmax_20
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_HISPANICS}/${CHILD_ID}_pileupmax_20_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_HISPANICS}/${CHILD_ID}_pileupmax_20_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        # hispanics
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_HISPANICS}/${CHILD_ID}_hispanics_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_HISPANICS}/${CHILD_ID}_hispanics_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        # hispanics_original
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_HISPANICS}/${CHILD_ID}_hispanics_original_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_HISPANICS}/${CHILD_ID}_hispanics_original_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        # k_8
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_HISPANICS}/${CHILD_ID}_k_8_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR_HISPANICS}/${CHILD_ID}_k_8_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
        
        # 2. De novo rate
        ## V1 original
        #DENOVO=$(java CheckDeNovoOld ${INPUT_DIR}/${CHILD_ID}_cohort_regenotyped_07_${SUFFIX}_gtmatrix.txt)
        #ROW="${ROW},${DENOVO}"
        # V1 latest run
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_HISPANICS}/${CHILD_ID}_v1_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # minkfreq_4
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_HISPANICS}/${CHILD_ID}_minkfreq_4_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # minkfreq_6
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_HISPANICS}/${CHILD_ID}_minkfreq_6_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # pileupmax_20
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_HISPANICS}/${CHILD_ID}_pileupmax_20_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # hispanics
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_HISPANICS}/${CHILD_ID}_hispanics_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # hispanics_original
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_HISPANICS}/${CHILD_ID}_hispanics_original_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
        # k_8
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_HISPANICS}/${CHILD_ID}_k_8_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"

        echo "${ROW}" >> triohispanics_${SUFFIX}.csv
    done
done < children.txt
rm -f children.txt
