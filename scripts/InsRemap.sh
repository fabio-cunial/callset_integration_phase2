#!/bin/bash
#
INPUT_DIR=$1

set -euxo pipefail


function CountClassifications() {
    MIN_SVLEN=$1
    MAX_SVLEN=$2

    rm -f remap_classification_${MIN_SVLEN}_${MAX_SVLEN}.csv
    for FILE in ${INPUT_DIR}/*.csv ; do
        awk -v MIN_SVLEN=${MIN_SVLEN} -v MAX_SVLEN=${MAX_SVLEN} 'BEGIN { FS=","; OFS=","; COUNT_ALL=0; COUNT_NOVEL=0; COUNT_TANDEM=0; COUNT_TANDEM_INVERTED=0; COUNT_TANDEM_COMPLEX=0; COUNT_INTERSPERSED=0; COUNT_PARTIAL=0; } { \
            if ($4<MIN_SVLEN || $4>=MAX_SVLEN) { next } \
            COUNT_ALL++ ; \
            if ($5=="novel") { COUNT_NOVEL++ } \
            else if ($5=="tandem") { COUNT_TANDEM++ } \
            else if ($5=="tandem_inverted") { COUNT_TANDEM_INVERTED++ } \
            else if ($5=="tandem_complex") { COUNT_TANDEM_COMPLEX++ } \
            else if ($5=="interspersed") { COUNT_INTERSPERSED++ } \
            else if ($5=="partial") { COUNT_PARTIAL++ } \
        } END { print COUNT_ALL, COUNT_NOVEL, COUNT_TANDEM, COUNT_TANDEM_INVERTED, COUNT_TANDEM_COMPLEX, COUNT_INTERSPERSED, COUNT_PARTIAL }' ${FILE} >> remap_classification_${MIN_SVLEN}_${MAX_SVLEN}.csv
    done
}


function CountTandems() {
    rm -f remap_tandems.csv
    for FILE in ${INPUT_DIR}/*.csv ; do
        awk 'BEGIN { FS=","; OFS=","; COUNT_ALL=0; COUNT_SNIFFLES=0; COUNT_PBSV=0; COUNT_PAV=0; } { \
            if ($5!="tandem") { next } \
            COUNT_ALL++ ; \
            if ($1==1) { COUNT_SNIFFLES++ } \
            if ($2==1) { COUNT_PBSV++ } \
            if ($3==1) { COUNT_PAV++ } \
        } END { print COUNT_ALL, COUNT_SNIFFLES, COUNT_PBSV, COUNT_PAV }' ${FILE} >> remap_tandems.csv
    done
}


function CountCaller() {
    CALLER=$1

    rm -f remap_${CALLER}.csv
    for FILE in ${INPUT_DIR}/*.csv ; do
        awk -v CALLER=${CALLER} 'BEGIN { FS=","; OFS=","; COUNT_ALL=0; COUNT_NOVEL=0; COUNT_TANDEM=0; COUNT_TANDEM_INVERTED=0; COUNT_TANDEM_COMPLEX=0; COUNT_INTERSPERSED=0; COUNT_PARTIAL=0; } { \
            if ( (CALLER=="sniffles" && $1==1) || (CALLER=="pbsv" && $2==1) || (CALLER=="pav" && $3==1) ) { \
                COUNT_ALL++ ; \
                if ($5=="novel") { COUNT_NOVEL++ } \
                else if ($5=="tandem") { COUNT_TANDEM++ } \
                else if ($5=="tandem_inverted") { COUNT_TANDEM_INVERTED++ } \
                else if ($5=="tandem_complex") { COUNT_TANDEM_COMPLEX++ } \
                else if ($5=="interspersed") { COUNT_INTERSPERSED++ } \
                else if ($5=="partial") { COUNT_PARTIAL++ } \
            } \
        } END { print COUNT_ALL, COUNT_NOVEL, COUNT_TANDEM, COUNT_TANDEM_INVERTED, COUNT_TANDEM_COMPLEX, COUNT_INTERSPERSED, COUNT_PARTIAL }' ${FILE} >> remap_${CALLER}.csv
    done
}


function CountTandemLengths() {
    rm -f remap_tandem_lengths.csv
    for FILE in ${INPUT_DIR}/*.csv ; do
        awk 'BEGIN { FS=","; OFS=","; COUNT_ALL=0; COUNT_1=0; COUNT_2=0; COUNT_3=0; COUNT_4=0; } { \
            if ($5!="tandem") { next } \
            COUNT_ALL++; \
            if ($4>=10000 && $4<20000) { COUNT_1++ } \
            else if ($4>=20000 && $4<50000) { COUNT_2++ } \
            else if ($4>=50000 && $4<100000) { COUNT_3++ } \
            else if ($4>=100000) { COUNT_4++ } \
        } END { print COUNT_ALL, COUNT_1, COUNT_2, COUNT_3, COUNT_4 }' ${FILE} >> remap_tandem_lengths.csv
    done
}


function CountTandemByLength() {
    rm -f remap_tandem_by_length.csv
    for FILE in ${INPUT_DIR}/*.csv ; do
        awk 'BEGIN { FS=","; OFS=","; COUNT_1_ALL=0; COUNT_1=0; COUNT_2_ALL=0; COUNT_2=0; COUNT_3_ALL=0; COUNT_3=0; COUNT_4_ALL=0; COUNT_4=0; } { \
            if ($4>=10000 && $4<20000) { \
                COUNT_1_ALL++; \
                if ($5=="tandem") { COUNT_1++ } \
            } else if ($4>=20000 && $4<50000) { \
                COUNT_2_ALL++; \
                if ($5=="tandem") { COUNT_2++ } \
            } else if ($4>=50000 && $4<100000) { \
                COUNT_3_ALL++; \
                if ($5=="tandem") { COUNT_3++ } \
            } else if ($4>=100000) { \
                COUNT_4_ALL++; \
                if ($5=="tandem") { COUNT_4++ } \
            } \
        } END { print COUNT_1_ALL, COUNT_1, COUNT_2_ALL, COUNT_2, COUNT_3_ALL, COUNT_3, COUNT_4_ALL, COUNT_4 }' ${FILE} >> remap_tandem_by_length.csv
    done
}


CountClassifications 0 10000000
CountClassifications 10000 20000
CountClassifications 20000 50000
CountClassifications 50000 100000
CountTandems
CountCaller sniffles
CountCaller pbsv
CountCaller pav
CountTandemLengths
CountTandemByLength