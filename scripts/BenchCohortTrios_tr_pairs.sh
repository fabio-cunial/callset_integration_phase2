#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios_only_5.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_tr_pairs"

ML_IDS="1 2 4 8 16 32 64 128"
NR_IDS="1 2 4 8 16 32 64"
IL_IDS="2 4 6 8 10 12 14 16 18 20 40 80 160 320"

set -euxo pipefail




# ------------------------------ Count matrices --------------------------------
rm -f ml_nr_counts.csv
for IDI in ${ML_IDS}; do
    ROW=""
    for IDJ in ${NR_IDS}; do
        if [ -e ${INPUT_DIR}/ml_nr/ml_${IDI}_nr_${IDJ}_count.txt ]; then
            COUNT=$(cat ${INPUT_DIR}/ml_nr/ml_${IDI}_nr_${IDJ}_count.txt)
            ROW="${ROW}${COUNT},"
        else
            ROW="${ROW}0,"
        fi
    done
    echo ${ROW} >> ml_nr_counts.csv
done

rm -f ml_il_counts.csv
for IDI in ${ML_IDS}; do
    ROW=""
    for IDJ in ${IL_IDS}; do
        if [ -e ${INPUT_DIR}/ml_il/ml_${IDI}_il_${IDJ}_count.txt ]; then
            COUNT=$(cat ${INPUT_DIR}/ml_il/ml_${IDI}_il_${IDJ}_count.txt)
            ROW="${ROW}${COUNT},"
        else
            ROW="${ROW}0,"
        fi
    done
    echo ${ROW} >> ml_il_counts.csv
done

rm -f nr_il_counts.csv
for IDI in ${NR_IDS}; do
    ROW=""
    for IDJ in ${IL_IDS}; do
        if [ -e ${INPUT_DIR}/nr_il/nr_${IDI}_il_${IDJ}_count.txt ]; then
            COUNT=$(cat ${INPUT_DIR}/nr_il/nr_${IDI}_il_${IDJ}_count.txt)
            ROW="${ROW}${COUNT},"
        else
            ROW="${ROW}0,"
        fi
    done
    echo ${ROW} >> nr_il_counts.csv
done




# ----------------------------- Error matrices ---------------------------------
cut -f 2 ${TRIOS_TSV} > children.txt

rm -f ml_nr_mendelian.csv
for IDI in ${ML_IDS}; do
    ROW=""
    for IDJ in ${NR_IDS}; do
        SUM="scale=8; 0"; N_CHILDREN="0";
        while read CHILD_ID; do
            FILE="${INPUT_DIR}/ml_nr/${CHILD_ID}_ml_${IDI}_nr_${IDJ}.txt"
            if [ -e ${FILE} ]; then
                N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
                N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
                if [ ${N_MERR} -eq 0 ]; then
                    VALUE="0"
                else
                    VALUE=$(echo "scale=8; ${N_MERR}/(${N_MERR}+${N_GOOD_ALT})" | bc)
                fi
                SUM="${SUM} + ${VALUE}"
                N_CHILDREN=$(( ${N_CHILDREN} + 1 ))
            else
                SUM="${SUM} + 0"
            fi
        done < children.txt
        NUMERATOR=$(echo ${SUM} | bc)
        if [ ${NUMERATOR} -eq 0 ]; then
            AVG="0"
        else
            AVG="scale=8; (${NUMERATOR} / ${N_CHILDREN})"
            AVG=$(echo ${AVG} | bc)
        fi
        ROW="${ROW}${AVG},"
    done
    echo ${ROW} >> ml_nr_mendelian.csv
done

rm -f ml_il_mendelian.csv
for IDI in ${ML_IDS}; do
    ROW=""
    for IDJ in ${IL_IDS}; do
        SUM="scale=8; 0"; N_CHILDREN="0";
        while read CHILD_ID; do
            FILE="${INPUT_DIR}/ml_il/${CHILD_ID}_ml_${IDI}_il_${IDJ}.txt"
            if [ -e ${FILE} ]; then
                N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
                N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
                if [ ${N_MERR} -eq 0 ]; then
                    VALUE="0"
                else
                    VALUE=$(echo "scale=8; ${N_MERR}/(${N_MERR}+${N_GOOD_ALT})" | bc)
                fi
                SUM="${SUM} + ${VALUE}"
                N_CHILDREN=$(( ${N_CHILDREN} + 1 ))
            else
                SUM="${SUM} + 0"
            fi
        done < children.txt
        NUMERATOR=$(echo ${SUM} | bc)
        if [ ${NUMERATOR} -eq 0 ]; then
            AVG="0"
        else
            AVG="scale=8; (${NUMERATOR} / ${N_CHILDREN})"
            AVG=$(echo ${AVG} | bc)
        fi
        ROW="${ROW}${AVG},"
    done
    echo ${ROW} >> ml_il_mendelian.csv
done

rm -f nr_il_mendelian.csv
for IDI in ${NR_IDS}; do
    ROW=""
    for IDJ in ${IL_IDS}; do
        SUM="scale=8; 0"; N_CHILDREN="0";
        while read CHILD_ID; do
            FILE="${INPUT_DIR}/nr_il/${CHILD_ID}_nr_${IDI}_il_${IDJ}.txt"
            if [ -e ${FILE} ]; then
                N_GOOD_ALT=$(grep ^ngood_alt ${FILE} | cut -f 2)
                N_MERR=$(grep ^nmerr ${FILE} | cut -f 2)
                if [ ${N_MERR} -eq 0 ]; then
                    VALUE="0"
                else
                    VALUE=$(echo "scale=8; ${N_MERR}/(${N_MERR}+${N_GOOD_ALT})" | bc)
                fi
                SUM="${SUM} + ${VALUE}"
                N_CHILDREN=$(( ${N_CHILDREN} + 1 ))
            else
                SUM="${SUM} + 0"
            fi
        done < children.txt
        NUMERATOR=$(echo ${SUM} | bc)
        if [ ${NUMERATOR} -eq 0 ]; then
            AVG="0"
        else
            AVG="scale=8; (${NUMERATOR} / ${N_CHILDREN})"
            AVG=$(echo ${AVG} | bc)
        fi
        ROW="${ROW}${AVG},"
    done
    echo ${ROW} >> nr_il_mendelian.csv
done

rm -f children.txt
