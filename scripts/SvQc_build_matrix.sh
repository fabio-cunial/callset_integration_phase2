#!/bin/bash
#
REMOTE_INDIR_BUCKET="???"
REMOTE_SAMPLES_TSV_BUCKET="???"

REMOTE_ROOT_DIR="${REMOTE_INDIR_BUCKET}/scratch/cunial_svqc"
REMOTE_INDIR="${REMOTE_ROOT_DIR}/counts"
REMOTE_PLATFORM_DIR="${REMOTE_ROOT_DIR}"
REMOTE_OUTDIR="${REMOTE_ROOT_DIR}"
REMOTE_SAMPLES_TSV="${REMOTE_SAMPLES_TSV_BUCKET}/for_Fabio_Kiran_v9_QC_reports/cohort_sample_list/v9.public_facing.research_ids.tsv"

set -euxo pipefail

# Localizing input files
gcloud storage cp ${REMOTE_SAMPLES_TSV} ./tmp.tsv
tail -n +2 tmp.tsv | cut -f 1 > samples.tsv  # Discarding header
rm -f tmp.tsv
gcloud storage cp ${REMOTE_PLATFORM_DIR}/revio.txt .
gcloud storage cp ${REMOTE_PLATFORM_DIR}/sequel.txt .
gcloud storage cp ${REMOTE_PLATFORM_DIR}/male.txt .
gcloud storage cp ${REMOTE_PLATFORM_DIR}/female.txt .
gcloud storage cp ${REMOTE_INDIR}/'*.counts' .

# Building counts tables
set +x
FILES=""; FILES_REVIO=""; FILES_SEQUEL=""; FILES_MALE=""; FILES_FEMALE=""
while read LINE; do
    SAMPLE_ID=$(echo ${LINE} | cut -f 1)
    if [ -e ${SAMPLE_ID}.counts ]; then
        FILES="${FILES} ${SAMPLE_ID}.counts"
        if grep -q ${SAMPLE_ID} revio.txt; then
            FILES_REVIO="${FILES_REVIO} ${SAMPLE_ID}.counts"
        fi
        if grep -q ${SAMPLE_ID} sequel.txt; then
            FILES_SEQUEL="${FILES_SEQUEL} ${SAMPLE_ID}.counts"
        fi
        if grep -q ${SAMPLE_ID} male.txt; then
            FILES_MALE="${FILES_MALE} ${SAMPLE_ID}.counts"
        fi
        if grep -q ${SAMPLE_ID} female.txt; then
            FILES_FEMALE="${FILES_FEMALE} ${SAMPLE_ID}.counts"
        fi
    fi
done < samples.tsv
set -x
if [ -n "${FILES}" ]; then cat ${FILES} > matrix_all.csv; else touch matrix_all.csv; fi
if [ -n "${FILES}" ]; then cat ${FILES} | cut -d , -f 2- > counts_all.csv; else touch counts_all.csv; fi
if [ -n "${FILES_REVIO}" ]; then cat ${FILES_REVIO} | cut -d , -f 2- > counts_revio.csv; else touch counts_revio.csv; fi
if [ -n "${FILES_SEQUEL}" ]; then cat ${FILES_SEQUEL} | cut -d , -f 2- > counts_sequel.csv; else touch counts_sequel.csv; fi
if [ -n "${FILES_MALE}" ]; then cat ${FILES_MALE} | cut -d , -f 2- > counts_male.csv; else touch counts_male.csv; fi
if [ -n "${FILES_FEMALE}" ]; then cat ${FILES_FEMALE} | cut -d , -f 2- > counts_female.csv; else touch counts_female.csv; fi

gcloud storage cp 'counts_*.csv' ${REMOTE_OUTDIR}/
