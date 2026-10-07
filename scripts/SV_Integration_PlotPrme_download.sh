#!/bin/bash
#
REMOTE_DIR="gs://_____________________/scratch/cunial_intersample_vcf/v3/manuscript/prme_analysis_truvari_collapse_vcf"  # prme_analysis_final_vcf    prme_analysis_truvari_collapse_vcf

set -euxo pipefail

mkdir -p ./precision_recall ./mendelian
gcloud storage cp ${REMOTE_DIR}/precision_recall/'*' ./precision_recall/
gcloud storage cp ${REMOTE_DIR}/mendelian/'*_mendelian_*' ./mendelian/
gcloud storage cp ${REMOTE_DIR}/mendelian/'*_dnm*' ./mendelian/
