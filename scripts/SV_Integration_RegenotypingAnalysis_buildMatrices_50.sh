#!/bin/bash
#
N_SAMPLES="1 2 4 8 16 32 64 128 256 512 1024 2048"

PRECISION_RECALL_SAMPLES_15x="HG00097 HG01890 HG00272 NA18534"
PRECISION_RECALL_SAMPLES_30x="HG01123 NA18943 HG02015 HG02984"

MENDELIAN_ERROR_SAMPLES_15x_CONTROL="HG00514 HG00733 NA19240"
MENDELIAN_ERROR_SAMPLES_15x_AOU="3498199 1665275 2611029 2656909"
MENDELIAN_ERROR_SAMPLES_30x_AOU="1981828 1113963 1458572 1752535"


set -euxo pipefail

# ---------------------------- Precision/recall --------------------------------
for REGION in all tr not_tr; do
    rm -f precision_recall_15x_50bp_${REGION}.csv
    for SAMPLE in ${PRECISION_RECALL_SAMPLES_15x}; do
        FILE="truvari/precision_recall/${SAMPLE}_truvari_50bp_${REGION}.txt"
        P=$(grep precision ${FILE} | cut -w -f 3)
        R=$(grep recall ${FILE} | cut -w -f 3)
        F=$(grep f1 ${FILE} | cut -w -f 3)
        C=$(grep gt_concordance ${FILE} | cut -w -f 3)
        ROW="${P}${R}${F}${C}"
        for MIN in ${N_SAMPLES}; do
            FILE="${MIN}_samples/precision_recall/${SAMPLE}_kanpig_50bp_${REGION}.txt"
            P=$(grep precision ${FILE} | cut -w -f 3)
            R=$(grep recall ${FILE} | cut -w -f 3)
            F=$(grep f1 ${FILE} | cut -w -f 3)
            C=$(grep gt_concordance ${FILE} | cut -w -f 3)
            ROW="${ROW}${P}${R}${F}${C}"
        done
        echo ${ROW} >> precision_recall_15x_50bp_${REGION}.csv
    done
    
    rm -f precision_recall_30x_50bp_${REGION}.csv
    for SAMPLE in ${PRECISION_RECALL_SAMPLES_30x}; do
        FILE="truvari/precision_recall/${SAMPLE}_truvari_50bp_${REGION}.txt"
        P=$(grep precision ${FILE} | cut -w -f 3)
        R=$(grep recall ${FILE} | cut -w -f 3)
        F=$(grep f1 ${FILE} | cut -w -f 3)
        C=$(grep gt_concordance ${FILE} | cut -w -f 3)
        ROW="${P}${R}${F}${C}"
        for MIN in ${N_SAMPLES}; do
            FILE="${MIN}_samples/precision_recall/${SAMPLE}_kanpig_50bp_${REGION}.txt"
            P=$(grep precision ${FILE} | cut -w -f 3)
            R=$(grep recall ${FILE} | cut -w -f 3)
            F=$(grep f1 ${FILE} | cut -w -f 3)
            C=$(grep gt_concordance ${FILE} | cut -w -f 3)
            ROW="${ROW}${P}${R}${F}${C}"
        done
        echo ${ROW} >> precision_recall_30x_50bp_${REGION}.csv
    done
done

# ---------------------------- Mendelian error ---------------------------------
for REGION in all tr not_tr; do
    rm -f mendelian_error_15x_50bp_control_${REGION}.csv
    for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_CONTROL}; do
        FILE="truvari/mendelian/${SAMPLE}_mendelian_50bp_${REGION}.txt"
        N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
        N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
        ROW="${N_GOOD_ALT},${N_MERR},"
        for MIN in ${N_SAMPLES}; do
            FILE="${MIN}_samples/mendelian/${SAMPLE}_mendelian_50bp_${REGION}.txt"
            N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
            N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
            ROW="${ROW}${N_GOOD_ALT},${N_MERR},"
        done
        echo ${ROW} >> mendelian_error_15x_50bp_control_${REGION}.csv
    done
    
    rm -f mendelian_error_15x_50bp_aou_${REGION}.csv
    for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_AOU}; do
        FILE="truvari/mendelian/${SAMPLE}_mendelian_50bp_${REGION}.txt"
        N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
        N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
        ROW="${N_GOOD_ALT},${N_MERR},"
        for MIN in ${N_SAMPLES}; do
            FILE="${MIN}_samples/mendelian/${SAMPLE}_mendelian_50bp_${REGION}.txt"
            N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
            N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
            ROW="${ROW}${N_GOOD_ALT},${N_MERR},"
        done
        echo ${ROW} >> mendelian_error_15x_50bp_aou_${REGION}.csv
    done
    
    rm -f mendelian_error_30x_50bp_aou_${REGION}.csv
    for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_30x_AOU}; do
        FILE="truvari/mendelian/${SAMPLE}_mendelian_50bp_${REGION}.txt"
        N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
        N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
        ROW="${N_GOOD_ALT},${N_MERR},"
        for MIN in ${N_SAMPLES}; do
            FILE="${MIN}_samples/mendelian/${SAMPLE}_mendelian_50bp_${REGION}.txt"
            N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
            N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
            ROW="${ROW}${N_GOOD_ALT},${N_MERR},"
        done
        echo ${ROW} >> mendelian_error_30x_50bp_aou_${REGION}.csv
    done
done

# ------------------------------ De novo rate ----------------------------------
for REGION in all tr not_tr; do
    rm -f denovo_15x_50bp_control_${REGION}.csv
    for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_CONTROL}; do
        FILE="truvari/mendelian/${SAMPLE}_dnm2_50bp_${REGION}.txt"
        ROW=$(cat ${FILE})
        for MIN in ${N_SAMPLES}; do
            FILE="${MIN}_samples/mendelian/${SAMPLE}_dnm2_50bp_${REGION}.txt"
            ROW="${ROW},$(cat ${FILE})"
        done
        echo ${ROW} >> denovo_15x_50bp_control_${REGION}.csv
    done
    
    rm -f denovo_15x_50bp_aou_${REGION}.csv
    for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_15x_AOU}; do
        FILE="truvari/mendelian/${SAMPLE}_dnm2_50bp_${REGION}.txt"
        ROW=$(cat ${FILE})
        for MIN in ${N_SAMPLES}; do
            FILE="${MIN}_samples/mendelian/${SAMPLE}_dnm2_50bp_${REGION}.txt"
            ROW="${ROW},$(cat ${FILE})"
        done
        echo ${ROW} >> denovo_15x_50bp_aou_${REGION}.csv
    done
    
    rm -f denovo_30x_50bp_aou_${REGION}.csv
    for SAMPLE in ${MENDELIAN_ERROR_SAMPLES_30x_AOU}; do
        FILE="truvari/mendelian/${SAMPLE}_dnm2_50bp_${REGION}.txt"
        ROW=$(cat ${FILE})
        for MIN in ${N_SAMPLES}; do
            FILE="${MIN}_samples/mendelian/${SAMPLE}_dnm2_50bp_${REGION}.txt"
            ROW="${ROW},$(cat ${FILE})"
        done
        echo ${ROW} >> denovo_30x_50bp_aou_${REGION}.csv
    done
done
