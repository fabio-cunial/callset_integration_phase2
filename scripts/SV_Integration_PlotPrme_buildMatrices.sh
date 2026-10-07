#!/bin/bash
#
PRECISION_RECALL_SAMPLES_15x="HG00097 HG00272 HG00408 HG00733 HG01928 HG02572 HG02602 HG02922 HG03098 HG03688 HG03742 HG03804 HG06807 NA18505 NA18565 NA18940 NA19036 NA19240 NA19682 NA19700 NA20503 NA20850 NA21309 HG00512 HG00513 HG00514 HG00731 HG00732 HG01890 HG03683 NA18534 NA18939 NA19238 NA19239 NA19317 NA19650 NA19705 NA20509 NA20847 NA21487"
PRECISION_RECALL_SAMPLES_30x="HG01123 HG01530 HG01884 HG02015 HG02155 HG01573 HG02587 HG02953 HG03009 NA12329"

MENDELIAN_ERROR_SAMPLES_15x_CONTROL="HG00514 HG00733 NA19240"  # Children only
MENDELIAN_ERROR_SAMPLES_15x_AOU=" "
MENDELIAN_ERROR_SAMPLES_30x_AOU=" "


set -euxo pipefail


# ---------------------------- Precision/recall --------------------------------

for REGION in all tr not_tr; do
    rm -f precision_recall_20bp_49bp_${REGION}.csv
    for SAMPLE in ${PRECISION_RECALL_SAMPLES_15x}; do
        FILE="./precision_recall/${SAMPLE}_20bp_49bp_${REGION}.txt"
        P=$(grep precision ${FILE} | cut -w -f 3)
        R=$(grep recall ${FILE} | cut -w -f 3)
        F=$(grep f1 ${FILE} | cut -w -f 3)
        C=$(grep gt_concordance ${FILE} | cut -w -f 3)
        ROW="${P}${R}${F}${C}"
        echo ${ROW} >> precision_recall_20bp_49bp_${REGION}.csv
    done
    for SAMPLE in ${PRECISION_RECALL_SAMPLES_30x}; do
        FILE="./precision_recall/${SAMPLE}_20bp_49bp_${REGION}.txt"
        P=$(grep precision ${FILE} | cut -w -f 3)
        R=$(grep recall ${FILE} | cut -w -f 3)
        F=$(grep f1 ${FILE} | cut -w -f 3)
        C=$(grep gt_concordance ${FILE} | cut -w -f 3)
        ROW="${P}${R}${F}${C}"
        echo ${ROW} >> precision_recall_20bp_49bp_${REGION}.csv
    done
    
    rm -f precision_recall_50bp_10000bp_${REGION}.csv
    for SAMPLE in ${PRECISION_RECALL_SAMPLES_15x}; do
        FILE="./precision_recall/${SAMPLE}_50bp_10000bp_${REGION}.txt"
        P=$(grep precision ${FILE} | cut -w -f 3)
        R=$(grep recall ${FILE} | cut -w -f 3)
        F=$(grep f1 ${FILE} | cut -w -f 3)
        C=$(grep gt_concordance ${FILE} | cut -w -f 3)
        ROW="${P}${R}${F}${C}"
        echo ${ROW} >> precision_recall_50bp_10000bp_${REGION}.csv
    done
    for SAMPLE in ${PRECISION_RECALL_SAMPLES_30x}; do
        FILE="./precision_recall/${SAMPLE}_50bp_10000bp_${REGION}.txt"
        P=$(grep precision ${FILE} | cut -w -f 3)
        R=$(grep recall ${FILE} | cut -w -f 3)
        F=$(grep f1 ${FILE} | cut -w -f 3)
        C=$(grep gt_concordance ${FILE} | cut -w -f 3)
        ROW="${P}${R}${F}${C}"
        echo ${ROW} >> precision_recall_50bp_10000bp_${REGION}.csv
    done
done


# ---------------------------- Mendelian error ---------------------------------

for REGION in all tr not_tr; do
    for SUFFIX in ${REGION} ${REGION}_no_missing; do
        rm -f mendelian_error_20bp_49bp_${SUFFIX}.csv
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_CONTROL}; do
            FILE="./mendelian/${SAMPLE}_mendelian_20bp_49bp_${SUFFIX}.txt"
            N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
            N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
            ROW="${N_GOOD_ALT},${N_MERR}"
            echo ${ROW} >> mendelian_error_20bp_49bp_${SUFFIX}.csv
        done
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_AOU}; do
            FILE="./mendelian/${SAMPLE}_mendelian_20bp_49bp_${SUFFIX}.txt"
            N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
            N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
            ROW="${N_GOOD_ALT},${N_MERR}"
            echo ${ROW} >> mendelian_error_20bp_49bp_${SUFFIX}.csv
        done
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_30x_AOU}; do
            FILE="./mendelian/${SAMPLE}_mendelian_20bp_49bp_${SUFFIX}.txt"
            N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
            N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
            ROW="${N_GOOD_ALT},${N_MERR}"
            echo ${ROW} >> mendelian_error_20bp_49bp_${SUFFIX}.csv
        done
        
        rm -f mendelian_error_50bp_10000bp_${SUFFIX}.csv
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_CONTROL}; do
            FILE="./mendelian/${SAMPLE}_mendelian_50bp_10000bp_${SUFFIX}.txt"
            N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
            N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
            ROW="${N_GOOD_ALT},${N_MERR}"
            echo ${ROW} >> mendelian_error_50bp_10000bp_${SUFFIX}.csv
        done
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_AOU}; do
            FILE="./mendelian/${SAMPLE}_mendelian_50bp_10000bp_${SUFFIX}.txt"
            N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
            N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
            ROW="${N_GOOD_ALT},${N_MERR}"
            echo ${ROW} >> mendelian_error_50bp_10000bp_${SUFFIX}.csv
        done
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_30x_AOU}; do
            FILE="./mendelian/${SAMPLE}_mendelian_50bp_10000bp_${SUFFIX}.txt"
            N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
            N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
            ROW="${N_GOOD_ALT},${N_MERR}"
            echo ${ROW} >> mendelian_error_50bp_10000bp_${SUFFIX}.csv
        done
    done
done


# ------------------------------ De novo rate ----------------------------------

for REGION in all tr not_tr; do
    for SUFFIX in ${REGION} ${REGION}_no_missing; do
        rm -f denovo_20bp_49bp_${SUFFIX}.csv
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_CONTROL}; do
            FILE="./mendelian/${SAMPLE}_dnm1_20bp_49bp_${SUFFIX}.txt"
            ROW=$(cat ${FILE})
            FILE="./mendelian/${SAMPLE}_dnm2_20bp_49bp_${SUFFIX}.txt"
            ROW="${ROW},$(cat ${FILE})"
            echo ${ROW} >> denovo_20bp_49bp_${SUFFIX}.csv
        done
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_AOU}; do
            FILE="./mendelian/${SAMPLE}_dnm1_20bp_49bp_${SUFFIX}.txt"
            ROW=$(cat ${FILE})
            FILE="./mendelian/${SAMPLE}_dnm2_20bp_49bp_${SUFFIX}.txt"
            ROW="${ROW},$(cat ${FILE})"
            echo ${ROW} >> denovo_20bp_49bp_${SUFFIX}.csv
        done
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_30x_AOU}; do
            FILE="./mendelian/${SAMPLE}_dnm1_20bp_49bp_${SUFFIX}.txt"
            ROW=$(cat ${FILE})
            FILE="./mendelian/${SAMPLE}_dnm2_20bp_49bp_${SUFFIX}.txt"
            ROW="${ROW},$(cat ${FILE})"
            echo ${ROW} >> denovo_20bp_49bp_${SUFFIX}.csv
        done
        
        rm -f denovo_50bp_10000bp_${SUFFIX}.csv
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_CONTROL}; do
            FILE="./mendelian/${SAMPLE}_dnm1_50bp_10000bp_${SUFFIX}.txt"
            ROW=$(cat ${FILE})
            FILE="./mendelian/${SAMPLE}_dnm2_50bp_10000bp_${SUFFIX}.txt"
            ROW="${ROW},$(cat ${FILE})"
            echo ${ROW} >> denovo_50bp_10000bp_${SUFFIX}.csv
        done
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_AOU}; do
            FILE="./mendelian/${SAMPLE}_dnm1_50bp_10000bp_${SUFFIX}.txt"
            ROW=$(cat ${FILE})
            FILE="./mendelian/${SAMPLE}_dnm2_50bp_10000bp_${SUFFIX}.txt"
            ROW="${ROW},$(cat ${FILE})"
            echo ${ROW} >> denovo_50bp_10000bp_${SUFFIX}.csv
        done
        for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_30x_AOU}; do
            FILE="./mendelian/${SAMPLE}_dnm1_50bp_10000bp_${SUFFIX}.txt"
            ROW=$(cat ${FILE})
            FILE="./mendelian/${SAMPLE}_dnm2_50bp_10000bp_${SUFFIX}.txt"
            ROW="${ROW},$(cat ${FILE})"
            echo ${ROW} >> denovo_50bp_10000bp_${SUFFIX}.csv
        done
    done
done
