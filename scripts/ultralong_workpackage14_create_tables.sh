#!/bin/bash
#
# Creates the input Terra tables for `SV_Integration_Workpackage14.wdl` based on
# the output of `SV_Integration_Workpackage13.wdl`.
#
BUCKET="???"
REMOTE_WORKPACKAGE_13_PREFIX="${BUCKET}/scratch/cunial_intersample_vcf/v3/workpackage_13_ultralong_filtered"
REMOTE_OUTPREFIX="${BUCKET}/scratch/cunial_intersample_vcf/v3/workpackages/workpackage14_ultralong_filtered"
LOCAL_OUTDIR="/Users/fcunial/Downloads/ultralong_workpackage_14"
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

set -euxo pipefail


for TYPE in all lenient stringent ; do
    rm -rf ${LOCAL_OUTDIR}/${TYPE}/ ; mkdir -p ${LOCAL_OUTDIR}/${TYPE}/
    TABLE_TSV="${LOCAL_OUTDIR}/${TYPE}/sv_integration_hg38_workpackage14_ultralong_filtered_${TYPE}.tsv"
    rm -f ${TABLE_TSV}
    echo -e "entity:sv_integration_hg38_workpackage14_ultralong_filtered_${TYPE}_id\tchromosome_id\tchunk_ids_file" >> ${TABLE_TSV}
    id="0"
    for CHROMOSOME in ${CHROMOSOMES}; do
        gcloud storage ls ${REMOTE_WORKPACKAGE_13_PREFIX}_${TYPE}/${CHROMOSOME}/'chunk_*.bcf' > chunks.txt
        touch ${LOCAL_OUTDIR}/${TYPE}/workpackage_14_${CHROMOSOME}
        while read CHUNK_URI; do
            CHUNK_ID=$( basename ${CHUNK_URI} .bcf | cut -d '_' -f 2 )
            echo ${CHUNK_ID} >> ${LOCAL_OUTDIR}/${TYPE}/workpackage_14_${CHROMOSOME}
        done < chunks.txt
        echo -e "${id}\t${CHROMOSOME}\t${REMOTE_OUTPREFIX}_${TYPE}/workpackage_14_${CHROMOSOME}" >> ${TABLE_TSV}
        id=$(( ${id} + 1 ))
    done
done
