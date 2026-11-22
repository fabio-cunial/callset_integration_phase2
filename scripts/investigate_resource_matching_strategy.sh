#!/bin/bash
#
SAMPLE_IDS=" "
START_POS_DIR="./sp"
START_POS_AND_GIVEN_REPR_DIR="./spagr"

set -euxo pipefail


for SAMPLE_ID in ${SAMPLE_IDS}; do
    bcftools query --format '%CHROM\t%POS\t%SVTYPE\t%SCORE\n' ${START_POS_DIR}/${SAMPLE_ID}_scored.bcf > ${START_POS_DIR}/${SAMPLE_ID}_scored.tsv
    bcftools query --format '%CHROM\t%POS\t%SVTYPE\t%SCORE\n' ${START_POS_AND_GIVEN_REPR_DIR}/${SAMPLE_ID}_scored.bcf > ${START_POS_AND_GIVEN_REPR_DIR}/${SAMPLE_ID}_scored.tsv
    diff ${START_POS_DIR}/${SAMPLE_ID}_scored.tsv ${START_POS_AND_GIVEN_REPR_DIR}/${SAMPLE_ID}_scored.tsv
done
