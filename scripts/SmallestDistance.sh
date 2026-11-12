#!/bin/bash
#
# Collects and plots the output of `SmallestDistance.wdl`.
#
REMOTE_BUCKET="..."
REMOTE_TR_BED="${REMOTE_BUCKET}/scratch/cunial_intersample_vcf/v2/plot_hwe/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed"
MIN_OVERLAP_FRACTION="0.9"

set -euxo pipefail


function buildInput() {
    local MIN_LENGTH=$1

    REMOTE_INDIR="${REMOTE_BUCKET}/scratch/cunial_intersample_vcf/v3/smallest_distance_analysis_minaf_0_minlength_${MIN_LENGTH}"
    mkdir ./chromosomes/
    gcloud storage cp ${REMOTE_INDIR}/'*.bed' ./chromosomes/
    cat $(ls ./chromosomes/*.bed | sort -V) > all.bed
    rm -rf ./chromosomes/
    bedtools coverage -a all.bed -b tr.bed > coverage.tsv
    rm -f all.bed
    awk -F'\t' -v OFS='\t' -v min_overlap_fraction=${MIN_OVERLAP_FRACTION} '/^#/ { print; next } $NF >= min_overlap_fraction { NF -= 4; print }' coverage.tsv > ${MIN_LENGTH}_inside_tr.bed
    awk -F'\t' -v OFS='\t' -v min_overlap_fraction=${MIN_OVERLAP_FRACTION} '/^#/ { print; next } $NF < min_overlap_fraction { NF -= 4; print }' coverage.tsv > ${MIN_LENGTH}_outside_tr.bed
    rm -f coverage.tsv
}


# ----------------------------- Main program -----------------------------------

gcloud storage cp ${REMOTE_TR_BED} ./tr.bed
buildInput 0
buildInput 50
python SmallestDistance.py output.png " " 0_inside_tr.bed 0_outside_tr.bed 50_inside_tr.bed 50_outside_tr.bed
