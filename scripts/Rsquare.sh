#!/bin/bash
#
# Collects and plots the output of `Rsquare.wdl`.
#
RUN_ID="all_50_49_1000000_0.001"
REMOTE_BUCKET="..."
REMOTE_INDIR="${REMOTE_BUCKET}/scratch/cunial_intersample_vcf/v3/rsquare_analysis"
REMOTE_TR_BED="${REMOTE_BUCKET}/scratch/cunial_intersample_vcf/v2/plot_hwe/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed"
MIN_OVERLAP_FRACTION="0.9"
N_SAMPLES="12554"

set -euxo pipefail


function buildInput() {
    mkdir ./chromosomes/
    gcloud storage cp ${REMOTE_INDIR}/'chr*_'${RUN_ID}'.bed' ./chromosomes/
    cat $(ls ./chromosomes/*.bed | sort -V) > all.bed
    rm -rf ./chromosomes/
    bedtools coverage -a all.bed -b tr.bed > coverage.tsv
    rm -f all.bed
    awk -F'\t' -v OFS='\t' -v min_overlap_fraction=${MIN_OVERLAP_FRACTION} '/^#/ { print; next } $NF >= min_overlap_fraction { NF -= 4; print }' coverage.tsv > inside_tr.bed
    awk -F'\t' -v OFS='\t' -v min_overlap_fraction=${MIN_OVERLAP_FRACTION} '/^#/ { print; next } $NF < min_overlap_fraction { NF -= 4; print }' coverage.tsv > outside_tr.bed
    rm -f coverage.tsv
}


# ----------------------------- Main program -----------------------------------

gcloud storage cp ${REMOTE_TR_BED} ./tr.bed
buildInput
for SVTYPE in DEL INS; do
    python Rsquare.py inside_tr.bed outside_tr.bed ${N_SAMPLES} ${SVTYPE} --out output_${SVTYPE}.png --x-bins 64 --y-bins 40 --title-prefix "${SVTYPE},"
done
