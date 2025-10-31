#!/bin/bash
#
# Total samples in output: 12680
#
REMOTE_OUTPUT_DIR="${BI_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_4_all_centers"

REMOTE_INPUT_DIR_BI="${BI_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_4"
REMOTE_INPUT_DIR_HA="${HA_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_4"
REMOTE_INPUT_DIR_BCM="${BCM_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_4"
REMOTE_INPUT_DIR_UW="${UW_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_4"
REMOTE_INPUT_DIR_CONTROLS_15X="${HPRC_BUCKET}/v2/15x/workpackage_4"
REMOTE_INPUT_DIR_CONTROLS_30X="${HPRC_BUCKET}/v2/30x/workpackage_4"

BI_SAMPLES_TO_COPY="1841126 1846203"
CONTROL_30X_SAMPLES_TO_COPY="HG01123 HG01530 HG01884 HG02015 HG02155 NA12329 HG01573 HG02587 HG02953 HG03009 HG00128 HG00344 HG00438 HG00658 HG02018 HG02717 HG02984 HG03874 HG04184 NA18943"


set -euxo pipefail

gsutil -m cp ${REMOTE_INPUT_DIR_BI}/'*.vcf.gz*' ${REMOTE_OUTPUT_DIR}/
gsutil -m cp ${REMOTE_INPUT_DIR_HA}/'*.vcf.gz*' ${REMOTE_OUTPUT_DIR}/
for SAMPLE_ID in ${BI_SAMPLES_TO_COPY}; do
    gsutil -m cp ${REMOTE_INPUT_DIR_BI}/${SAMPLE_ID}_'*.vcf.gz*' ${REMOTE_OUTPUT_DIR}/
done
gsutil -m cp ${REMOTE_INPUT_DIR_UW}/'*.vcf.gz*' ${REMOTE_OUTPUT_DIR}/
gsutil -m cp ${REMOTE_INPUT_DIR_BCM}/'*.vcf.gz*' ${REMOTE_OUTPUT_DIR}/
gsutil -m cp ${REMOTE_INPUT_DIR_CONTROLS_15X}/'*.vcf.gz*' ${REMOTE_OUTPUT_DIR}/
for SAMPLE_ID in ${CONTROL_30X_SAMPLES_TO_COPY}; do
    gsutil -m cp ${REMOTE_INPUT_DIR_CONTROLS_30X}/${SAMPLE_ID}_'*.vcf.gz*' ${REMOTE_OUTPUT_DIR}/
done
