#!/bin/bash
#
INPUT_TSV="hprc_y2_males.tsv"

set -euxo pipefail

while read ROW; do
    ID=$(echo "${ROW}" | cut -f 1)
    HAP1_BAM=$(echo "${ROW}" | cut -f 2)
    HAP2_BAM=$(echo "${ROW}" | cut -f 3)
    echo ${HAP1_BAM} > list.txt
    echo ${HAP1_BAM}.bai >> list.txt
    echo ${HAP2_BAM} >> list.txt
    echo ${HAP2_BAM}.bai >> list.txt
    cat list.txt | gcloud storage cp -I ./
    samtools coverage --region chrY $(basename ${HAP1_BAM}) > ${ID}_y_hap1.txt &
    samtools coverage --region chrY $(basename ${HAP2_BAM}) > ${ID}_y_hap2.txt &
    samtools coverage --region chrX $(basename ${HAP1_BAM}) > ${ID}_x_hap1.txt &
    samtools coverage --region chrX $(basename ${HAP2_BAM}) > ${ID}_x_hap2.txt &
    wait
    rm -rf *.bam *.bai
done < ${INPUT_TSV}
