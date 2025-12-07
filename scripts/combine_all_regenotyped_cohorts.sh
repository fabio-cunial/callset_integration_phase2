#!/bin/bash
#
# Total samples in output: 12680
#
BI_BUCKET="gs://"
HA_BUCKET="gs://"
BCM_BUCKET="gs://"
UW_BUCKET="gs://"

REMOTE_OUTPUT_DIR="${BI_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_9_all_centers"

REMOTE_INPUT_DIR_BI="${BI_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_9"
REMOTE_INPUT_DIR_HA="${HA_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_9"
REMOTE_INPUT_DIR_BCM="${BCM_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_9"
REMOTE_INPUT_DIR_UW="${UW_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_9"
REMOTE_INPUT_DIR_CONTROLS_15X="${BI_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_9_hprc_hgsvc_15x"
REMOTE_INPUT_DIR_CONTROLS_30X="${BI_BUCKET}/scratch/cunial_intersample_vcf/v2/workpackage_9_hprc_hgsvc_30x"

BI_SAMPLES_TO_COPY="1841126 1846203"
CONTROL_30X_SAMPLES_TO_COPY="HG01123 HG01530 HG01884 HG02015 HG02155 NA12329 HG01573 HG02587 HG02953 HG03009 HG00128 HG00344 HG00438 HG00658 HG02018 HG02717 HG02984 HG03874 HG04184 NA18943"
N_SAMPLES_BI="10023"
N_SAMPLES_HA="1234"
N_SAMPLES_BCM="756"
N_SAMPLES_UW="402"
N_CONTROLS_30X="20"
N_CONTROLS_15X=$(( 292 - ${N_CONTROLS_30X} ))


set -euxo pipefail

# Checking that the number of samples is correct
N_FILES=$(gsutil ls ${REMOTE_INPUT_DIR_BI}/'*_chunk_0.bcf' | wc -l)
if [ ${N_FILES} -ne ${N_SAMPLES_BI} ]; then
    echo "ERROR: BI has ${N_FILES} files != ${N_SAMPLES_BI}"
    exit
fi
N_FILES=$(gsutil ls ${REMOTE_INPUT_DIR_HA}/'*_chunk_0.bcf' | wc -l)
if [ ${N_FILES} -ne ${N_SAMPLES_HA} ]; then
    echo "ERROR: HA has ${N_FILES} files != ${N_SAMPLES_HA}"
    exit
fi
N_FILES=$(gsutil ls ${REMOTE_INPUT_DIR_BCM}/'*_chunk_0.bcf' | wc -l)
if [ ${N_FILES} -ne ${N_SAMPLES_BCM} ]; then
    echo "ERROR: BCM has ${N_FILES} files != ${N_SAMPLES_BCM}"
    exit
fi
N_FILES=$(gsutil ls ${REMOTE_INPUT_DIR_UW}/'*_chunk_0.bcf' | wc -l)
if [ ${N_FILES} -ne ${N_SAMPLES_UW} ]; then
    echo "ERROR: UW has ${N_FILES} files != ${N_SAMPLES_UW}"
    exit
fi
N_FILES=$(gsutil ls ${REMOTE_INPUT_DIR_CONTROLS_15X}/'*_chunk_0.bcf' | wc -l)
if [ ${N_FILES} -ne ${N_CONTROLS_15X} ]; then
    echo "ERROR: CONTROLS15X has ${N_FILES} files != ${N_CONTROLS_15X}"
    exit
fi
N_FILES=$(gsutil ls ${REMOTE_INPUT_DIR_CONTROLS_30X}/'*_chunk_0.bcf' | wc -l)
if [ ${N_FILES} -ne ${N_CONTROLS_30X} ]; then
    echo "ERROR: CONTROLS30X has ${N_FILES} files != ${N_CONTROLS_30X}"
    exit
fi

# Copying all files to a single directory
gsutil -m cp ${REMOTE_INPUT_DIR_BI}/'*.bcf*' ${REMOTE_OUTPUT_DIR}/
gsutil -m cp ${REMOTE_INPUT_DIR_HA}/'*.bcf*' ${REMOTE_OUTPUT_DIR}/
for SAMPLE_ID in ${BI_SAMPLES_TO_COPY}; do
    gsutil -m cp ${REMOTE_INPUT_DIR_BI}/${SAMPLE_ID}_'*.bcf*' ${REMOTE_OUTPUT_DIR}/
done
gsutil -m cp ${REMOTE_INPUT_DIR_UW}/'*.bcf*' ${REMOTE_OUTPUT_DIR}/
gsutil -m cp ${REMOTE_INPUT_DIR_BCM}/'*.bcf*' ${REMOTE_OUTPUT_DIR}/
gsutil -m cp ${REMOTE_INPUT_DIR_CONTROLS_15X}/'*.bcf*' ${REMOTE_OUTPUT_DIR}/
for SAMPLE_ID in ${CONTROL_30X_SAMPLES_TO_COPY}; do
    gsutil -m cp ${REMOTE_INPUT_DIR_CONTROLS_30X}/${SAMPLE_ID}_'*.bcf*' ${REMOTE_OUTPUT_DIR}/
done

# It is now safe to remove the input directories
