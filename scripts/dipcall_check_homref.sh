#!/bin/bash
#
# Given a set of dipcall VCFs on HPRC, checks how many calls have a non-ALT GT.
#
INPUT_DIR=$1

set -euo pipefail

for FILE in $(ls ${INPUT_DIR}/*.vcf.gz); do
    N_RECORDS=$(bcftools index --nrecords ${FILE}.tbi)
    N_ALT=$(bcftools filter --include 'GT="alt"' ${FILE} | grep -v ^# | wc -l | xargs)
    if [ ${N_ALT} -ne ${N_RECORDS} ]; then
        echo "ERROR: ${N_RECORDS},${N_ALT},${FILE}"
    fi
done
