#!/bin/bash
#
REMOTE_VCF_GZ_TBI=$1
REMOTE_COUNTS_DIR=$2  # Without final slash
PREFIX=$3

set -euxo pipefail

rm -f ${PREFIX}_counts.txt
gcloud storage cp ${REMOTE_VCF_GZ_TBI} ./all.vcf.gz.tbi
echo $(bcftools index --nrecords ./all.vcf.gz.tbi) >> ${PREFIX}_counts.txt
gcloud storage cp ${REMOTE_COUNTS_DIR}/'start_*.bed' ${REMOTE_COUNTS_DIR}/'end_*.bed' ${REMOTE_COUNTS_DIR}/'interval_*.bed' .
for FILE in start_tr.bed end_tr.bed interval_tr.bed  start_segdup.bed end_segdup.bed interval_segdup.bed  start_segdup_gt10kb.bed end_segdup_gt10kb.bed interval_segdup_gt10kb.bed ; do
    echo $(wc -l < ${FILE}) >> ${PREFIX}_counts.txt
done
rm -f *.bed *.tbi
