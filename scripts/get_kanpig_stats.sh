#!/bin/bash
#
REMOTE_DIR=$1

set -euxo pipefail


gcloud storage cp ${REMOTE_DIR}/'*_kanpig.csv' .
rm -f table1.csv table2.csv
ls *.csv > file_list.txt
while read FILE; do
    head -n 1 ${FILE} | cut -d , -f 1,2,3  >> table1.csv
    tail -n 1 ${FILE} | cut -d , -f 1,2,3  >> table2.csv
done < file_list.txt
