#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios_only_5.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_tr_stratifications2"
INPUT_DIR_V1="/Users/fcunial/Downloads/BenchCohortTrios_filtered_truvari2"

#IDS="AllHomopolymers    AllTandemRepeats AllTandemRepeats_le50bp AllTandemRepeats_51to200bp AllTandemRepeats_201to10000bp AllTandemRepeats_ge10001bp AllTandemRepeats_ge101bp    satellites microsatellite SimpleRepeat_diTR_10to49 SimpleRepeat_diTR_50to149 SimpleRepeat_diTR_ge150 SimpleRepeat_homopolymer_4to6 SimpleRepeat_homopolymer_7to11 SimpleRepeat_homopolymer_ge12 SimpleRepeat_homopolymer_ge21 SimpleRepeat_imperfecthomopolge11 SimpleRepeat_imperfecthomopolge21 SimpleRepeat_quadTR_19to49 SimpleRepeat_quadTR_50to149 SimpleRepeat_quadTR_ge150 SimpleRepeat_triTR_14to49 SimpleRepeat_triTR_50to149 SimpleRepeat_triTR_ge150    segdups segdups_gt10kb    numts    gc15 gc15to20 gc20to25 gc25to30 gc30to55 gc55to60 gc60to65 gc65to70 gc70to75 gc75to80 gc80to85 gc85    TRCompDB"
IDS="motifLength_1 motifLength_2 motifLength_4 motifLength_8 motifLength_16 motifLength_32 motifLength_64 motifLength_128 nRepeatsInRef_1 nRepeatsInRef_2 nRepeatsInRef_4 nRepeatsInRef_8 nRepeatsInRef_16 nRepeatsInRef_32 nRepeatsInRef_64 nRepeatsInRef_128 purity_05 purity_06 purity_07 purity_08 purity_09 purity_095 purity_096 purity_097 purity_098 purity_099 purity_1 lpsLengthStdev_0 lpsLengthStdev_1 lpsLengthStdev_2 lpsLengthStdev_4 lpsLengthStdev_8 lpsLengthStdev_16 lpsLengthStdev_32 lpsMotifFrac_01 lpsMotifFrac_02 lpsMotifFrac_03 lpsMotifFrac_04 lpsMotifFrac_05 lpsMotifFrac_06 lpsMotifFrac_07 lpsMotifFrac_08 lpsMotifFrac_09 lpsMotifFrac_095 lpsMotifFrac_096 lpsMotifFrac_097 lpsMotifFrac_098 lpsMotifFrac_099 lpsMotifFrac_1 stdevIllumina_0 stdevIllumina_1 stdevIllumina_2 stdevIllumina_4 stdevIllumina_8 stdevIllumina_16 stdevIllumina_32 stdevT2T_0 stdevT2T_1 stdevT2T_2 stdevT2T_4 stdevT2T_8 stdevT2T_16 stdevT2T_32        intervalLength_2 intervalLength_4 intervalLength_6 intervalLength_8 intervalLength_10 intervalLength_12 intervalLength_14 intervalLength_16 intervalLength_18 intervalLength_20 intervalLength_40 intervalLength_80 intervalLength_160 intervalLength_320 intervalLength_640 intervalLength_1280 intervalLength_2560"

MAX_AD="50"

set -euxo pipefail


# Counts in V1 07
rm -f tr_substrat_counts.txt
echo 0 >> tr_substrat_counts.txt  # Placeholder for V1
for PREFIX in ${IDS}; do
    if [ -e ${INPUT_DIR}/${PREFIX}_count.txt ]; then
        cat ${INPUT_DIR}/${PREFIX}_count.txt >> tr_substrat_counts.txt
    else
        echo 0 >> tr_substrat_counts.txt
    fi
done


# Errors
rm -f tr_substrat.csv
cut -f 2 ${TRIOS_TSV} > children.txt
while read CHILD_ID; do
    # 1. Mendelian error
    N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR_V1}/${CHILD_ID}_v1_tr.txt | cut -f 2)
    N_MERR=$(grep ^nmerr ${INPUT_DIR_V1}/${CHILD_ID}_v1_tr.txt | cut -f 2)
    ROW="${N_GOOD_ALT},${N_MERR}"
    for SUFFIX in ${IDS}; do
        N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_${SUFFIX}.txt | cut -f 2)
        N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_${SUFFIX}.txt | cut -f 2)
        ROW="${ROW},${N_GOOD_ALT},${N_MERR}"
    done
    
    # 2. De novo rate
    DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR_V1}/${CHILD_ID}_v1_tr_gtmatrix.txt ${MAX_AD})
    ROW="${ROW},${DENOVO}"
    for SUFFIX in ${IDS}; do
        DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR}/${CHILD_ID}_${SUFFIX}_gtmatrix.txt ${MAX_AD})
        ROW="${ROW},${DENOVO}"
    done
    
    echo "${ROW}" >> tr_substrat.csv
done < children.txt
rm -f children.txt
