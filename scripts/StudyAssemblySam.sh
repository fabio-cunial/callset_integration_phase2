#!/bin/bash
#
DISTANCE_THRESHOLD="1000"

N_THREADS="8"
RAM_MB="29000"
RAM_PER_THREAD_MB=$(( ${RAM_MB} / ${N_THREADS} ))

set -euxo pipefail

rm -f lengths.txt
while read -r ROW; do
    SAMPLE_ID=$(echo ${ROW} | cut -d , -f 1)
    HAP1_BAM=$(echo ${ROW} | cut -d , -f 2)
    HAP2_BAM=$(echo ${ROW} | cut -d , -f 3)

    gcloud storage cp ${HAP1_BAM} ./input.bam
    samtools sort -@ ${N_THREADS} -m ${RAM_PER_THREAD_MB}M -n -O sam input.bam > sorted.sam
    rm -f input.bam
    java -Xmx${RAM_MB}M StudyAssemblySam sorted.sam lengths.tsv ${SAMPLE_ID}_hap1_matrix.tsv ${DISTANCE_THRESHOLD}
    rm -f sorted.sam

    gcloud storage cp ${HAP2_BAM} ./input.bam
    samtools sort -@ ${N_THREADS} -m ${RAM_PER_THREAD_MB}M -n -O sam input.bam > sorted.sam
    rm -f input.bam
    java -Xmx${RAM_MB}M StudyAssemblySam sorted.sam lengths.tsv ${SAMPLE_ID}_hap2_matrix.tsv ${DISTANCE_THRESHOLD}
    rm -f sorted.sam
done < table.csv
