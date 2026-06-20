#!/bin/bash
#
# Copies some 30x files to a subdirectory to make integration easier in 
# `SV_Integration_Workpackage12.wdl`.
#
# Remark: the original files are kept, for safety.
#
BUCKET="gs://fc-7f861a33-ddb4-4b2f-8d10-5679c9df6108"
SAMPLES="HG01123 HG01530 HG01884 HG02015 HG02155 NA12329 HG01573 HG02587 HG02953 HG03009 HG00128 HG00344 HG00438 HG00658 HG02018 HG02717 HG02984 HG03874 HG04184 NA18943"

set -euxo pipefail


for SOURCE in all lenient stringent ; do
    REMOTE_INDIR="${BUCKET}/v3/30x/workpackage_3_ultralong_${SOURCE}"
    REMOTE_OUTDIR="${BUCKET}/v3/30x/workpackage_3_ultralong_${SOURCE}/for_integration"
    for SAMPLE_ID in ${SAMPLES} ; do
        gcloud storage cp ${REMOTE_INDIR}/${SAMPLE_ID}_${SOURCE}.'bcf*' ${REMOTE_OUTDIR}/
    done
done
