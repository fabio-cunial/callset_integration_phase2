#!/bin/bash
#
# Collects and plots the output of `NumNeighbors.wdl`.
#
REMOTE_BUCKET="..."
REMOTE_TR_BED="${REMOTE_BUCKET}/scratch/cunial_intersample_vcf/v2/plot_hwe/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed"
MIN_OVERLAP_FRACTION="0.9"

set -euxo pipefail


function buildInput() {
    local MIN_LENGTH=$1
    local SVTYPE=$2
    local TR_SUFFIX=$3
    local OUT_FILENAME=$4

    REMOTE_INDIR="${REMOTE_BUCKET}/scratch/cunial_intersample_vcf/v3/numneighbors_analysis_minlength_${MIN_LENGTH}"
    mkdir ./chromosomes/
    gcloud storage cp ${REMOTE_INDIR}/'*_'${SVTYPE}'_'${TR_SUFFIX}'.tsv' ./chromosomes/
    cat $(ls ./chromosomes/*.tsv | sort -V) > ${OUT_FILENAME}
    rm -rf ./chromosomes/
}


# ----------------------------- Main program -----------------------------------

gcloud storage cp ${REMOTE_TR_BED} ./tr.bed
for SVTYPE in INS DEL; do
    for TR_SUFFIX in tr not_tr; do
        for MIN_LENGTH in 20 50; do
            if [ ${MIN_LENGTH} -eq 20 ]; then
                LENGTH="0"
            else
                LENGTH=${MIN_LENGTH}
            fi
            OUT_FILENAME="${SVTYPE}_${TR_SUFFIX}_${MIN_LENGTH}.tsv"
            OUT_FILENAME="${OUT_FILENAME,,}"
            buildInput ${LENGTH} ${SVTYPE} ${TR_SUFFIX} ${OUT_FILENAME}
        done
    done
done
python NumNeighbors.py output.png " " .
