#!/bin/bash
#
INPUT_DIR="gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HGSVC3/Assemblies/fasta/hifiasm"
OUTPUT_DIR="gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HGSVC3/Assemblies/by_hap/"
SEQKIT_COMMAND="./seqkit"

set -euxo pipefail

gsutil ls ${INPUT_DIR} > assemblies.txt
while read FILE; do
    LOCAL_FILE=$(basename ${FILE})
    gsutil -m cp ${FILE} .
    zcat ${LOCAL_FILE} | grep '^>' > headers.txt
    grep '.h1tg' headers.txt > tmp.txt
    cat tmp.txt
    cut -d ' ' -f 1 tmp.txt | tr -d '>' > contigs1.txt
    rm -f tmp.txt
    grep '.h2tg' headers.txt > tmp.txt
    cat tmp.txt
    cut -d ' ' -f 1 tmp.txt | tr -d '>' > contigs2.txt
    rm -f tmp.txt
    rm -f headers.txt
    ID=${LOCAL_FILE%.fna.gz}
    ID=${LOCAL_FILE#*_}
    ID=${ID#*_}
    ID=${ID%%.*}
    ${SEQKIT_COMMAND} grep --pattern-file contigs1.txt ${LOCAL_FILE} | bgzip --compress-level 1 > ${ID}_hap1.fna.gz &
    ${SEQKIT_COMMAND} grep --pattern-file contigs2.txt ${LOCAL_FILE} | bgzip --compress-level 2 > ${ID}_hap1.fna.gz &
    wait
    gsutil -m mv ${ID}_hap'*'.fna.gz ${OUTPUT_DIR}
    rm -f ${LOCAL_FILE}.fna.gz contigs*.txt
done < assemblies.txt
