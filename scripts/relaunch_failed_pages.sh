#!/bin/bash
#
# 1. Click on a Terra submission with failed workflows, select "Completion 
#    status = Failed" and "Download TSV". Do this for every page with failed
#    workflows, and save all the TSVs in the same input directory.
# 2. Run this script, providing that input directory.
# 3. Import every generated TSV into the `table_name_set` table in the 
#    workspace (Data > Import data > Upload TSV).
# 4. Relaunch every set in Terra.
#
INPUT_DIR=$1
OUTPUT_DIR=$2

TABLE_NAME="sv_integration_hg38_workpackage1"

set -euxo pipefail

for FILE in ${INPUT_DIR}/*.tsv; do
    OUTPUT_FILE=${OUTPUT_DIR}/$(basename ${FILE})
    PAGE_ID=$(basename ${FILE} .tsv)
    PAGE_ID=${PAGE_ID#page}
    echo -e "membership:${TABLE_NAME}_set_id\t${TABLE_NAME}" > ${OUTPUT_FILE}
    cut -f 8 ${FILE} | tail -n +2 | sed "s/{\"entityType\":\"${TABLE_NAME}\",\"entityName\":\"//g" | sed 's/"}//g' | sed "s/^/failed_page_${PAGE_ID}$(printf '\t')/g" >> ${OUTPUT_FILE}
done
