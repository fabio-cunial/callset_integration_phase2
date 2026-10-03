#!/bin/bash
#
REMOTE_DIR="gs://fc-secure-95bbd6eb-6d63-49aa-a980-47f3c1342b1e/scratch/cunial_intersample_vcf/v3/manuscript/regenotyping_analysis"
N_SAMPLES="1 2 4 8 16 32 64 128 256 512 1024 2048"

set -euxo pipefail

mkdir -p truvari/precision_recall truvari/mendelian
gcloud storage cp ${REMOTE_DIR}/truvari/precision_recall/'*_truvari_*' ./truvari/precision_recall/
gcloud storage cp ${REMOTE_DIR}/truvari/mendelian/'*_mendelian_*' ./truvari/mendelian/
gcloud storage cp ${REMOTE_DIR}/truvari/mendelian/'*_dnm*' ./truvari/mendelian/
for N in ${N_SAMPLES}; do
    mkdir -p ${N}_samples/precision_recall ${N}_samples/mendelian
    gcloud storage cp ${REMOTE_DIR}/${N}_samples/precision_recall/'*_kanpig_*' ./${N}_samples/precision_recall/
    gcloud storage cp ${REMOTE_DIR}/${N}_samples/mendelian/'*_mendelian_*' ./${N}_samples/mendelian/
    gcloud storage cp ${REMOTE_DIR}/${N}_samples/mendelian/'*_dnm*' ./${N}_samples/mendelian/
done
