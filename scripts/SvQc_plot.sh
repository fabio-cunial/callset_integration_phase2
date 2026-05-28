#!/bin/bash
#
COHORT_ID="UW_ONT"
ROOT_DIR="/Users/fcunial/Downloads/svqc/plot/${COHORT_ID}"

set -euxo pipefail

for SUFFIX in all female male revio sequel; do
    for CALLER in 0 1 2; do
        for SVTYPE in DEL INS DUP INV; do
            python SvQc.py --input-dir ${ROOT_DIR} --cohort ${COHORT_ID} --suffix ${SUFFIX} --caller ${CALLER} --svtype ${SVTYPE}
        done
    done
done
