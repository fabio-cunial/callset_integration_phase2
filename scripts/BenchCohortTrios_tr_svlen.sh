#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios_only_5.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_tr_svlen"
SV_LENGTHS="50 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 1000000"

set -euxo pipefail


cut -f 2 ${TRIOS_TSV} > children.txt

# Mendelian error
for SUFFIX in del ins; do
    rm -f tr_svlen_${SUFFIX}.csv
    while read CHILD_ID; do
        ROW=""
        for SVLEN in ${SV_LENGTHS}; do
            N_GOOD_ALT=$(grep ^ngood_alt ${INPUT_DIR}/${CHILD_ID}_${SVLEN}_${SUFFIX}.txt | cut -f 2)
            N_MERR=$(grep ^nmerr ${INPUT_DIR}/${CHILD_ID}_${SVLEN}_${SUFFIX}.txt | cut -f 2)
            ROW="${ROW}${N_GOOD_ALT},${N_MERR},"
        done
        echo "${ROW}" >> tr_svlen_${SUFFIX}.csv
    done < children.txt
done

# Count
for SUFFIX in del ins; do
    rm -f tr_svlen_${SUFFIX}_count.csv
    while read CHILD_ID; do
        ROW=""
        for SVLEN in ${SV_LENGTHS}; do
            COUNT=$(cat ${INPUT_DIR}/${CHILD_ID}_${SVLEN}_${SUFFIX}_count.txt)
            ROW="${ROW}${COUNT},"
        done
        echo "${ROW}" >> tr_svlen_${SUFFIX}_count.csv
    done < children.txt
done

rm -f children.txt
