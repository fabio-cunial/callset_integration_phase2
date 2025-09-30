#!/bin/bash
#
TRIOS_TSV="/Users/fcunial/Downloads/trios_only_5.tsv"
INPUT_DIR="/Users/fcunial/Downloads/BenchCohortTrios_hispanics"

MAX_AD="50"
IDS="pileupmax_20"   #"only_50bp"  #"v1 beta_binomial"


set -euxo pipefail


rm -f ad_denovo_plot_all.csv ad_denovo_plot_tr.csv ad_denovo_plot_not_tr.csv
cut -f 2 ${TRIOS_TSV} > children.txt
while read CHILD_ID; do
    ROW=""
    for SUFFIX in all tr not_tr; do
        for ID in ${IDS}; do
            DENOVO=$(java CheckDeNovoNumneigh ${INPUT_DIR}/${CHILD_ID}_${ID}_${SUFFIX}_gtmatrix.txt ${MAX_AD})
            ROW="${ROW},${DENOVO}"
        done
        echo "${ROW}" >> ad_denovo_plot_${SUFFIX}.csv
    done
done < children.txt
rm -f children.txt
