#!/bin/bash
#
INPUT_DIR="/Users/fcunial/Downloads/denovo_by_region"
SAMPLES="1004931 1402381 1665275 1806012 3498199"

set -euxo pipefail


for SAMPLE in ${SAMPLES}; do
    for CHR in $(seq 1 22) X Y M; do
        cat ${INPUT_DIR}/${SAMPLE}_all.bed | awk -v chr="chr${CHR}" 'BEGIN { } { if ($1==chr) print $0 }' | cut -f 2- > ${INPUT_DIR}/${SAMPLE}_all_chr${CHR}.tsv
        cat ${INPUT_DIR}/${SAMPLE}_tr.bed | awk -v chr="chr${CHR}" 'BEGIN { } { if ($1==chr) print $0 }' | cut -f 2- > ${INPUT_DIR}/${SAMPLE}_tr_chr${CHR}.tsv
        cat ${INPUT_DIR}/${SAMPLE}_not_tr.bed | awk -v chr="chr${CHR}" 'BEGIN { } { if ($1==chr) print $0 }' | cut -f 2- > ${INPUT_DIR}/${SAMPLE}_not_tr_chr${CHR}.tsv
    done
done
