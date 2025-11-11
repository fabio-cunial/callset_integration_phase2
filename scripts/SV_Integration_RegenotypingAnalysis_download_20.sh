#!/bin/bash
#
REMOTE_DIR="${BUCKET_ROOT}/scratch/cunial_intersample_vcf/v2/regenotyping_analysis_v2_20bp"
N_SAMPLES="1 2 4 8 16 32 64 128 256 512 1024 2048"

set -euxo pipefail

mkdir -p truvari/precision_recall truvari/mendelian
gsutil cp ${REMOTE_DIR}/truvari/precision_recall/'*_truvari_*' ./truvari/precision_recall/
gsutil cp ${REMOTE_DIR}/truvari/mendelian/'*_mendelian_*' ./truvari/mendelian/
gsutil cp ${REMOTE_DIR}/truvari/mendelian/'*_dnm*' ./truvari/mendelian/
for N in ${N_SAMPLES}; do
    mkdir -p ${N}_samples/precision_recall ${N}_samples/mendelian
    gsutil cp ${REMOTE_DIR}/${N}_samples/precision_recall/'*_kanpig_*' ./${N}_samples/precision_recall/
    gsutil cp ${REMOTE_DIR}/${N}_samples/mendelian/'*_mendelian_*' ./${N}_samples/mendelian/
    gsutil cp ${REMOTE_DIR}/${N}_samples/mendelian/'*_dnm*' ./${N}_samples/mendelian/
done
