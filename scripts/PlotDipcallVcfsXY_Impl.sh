#!/bin/bash
#
MIN_SV_LENGTH=$1
MAX_SV_LENGTH=$2
CHROMOSOME=$3
SAMPLE_ID=$4


set -euxo pipefail


# Similar to the same function in
# `SV_Integration_RegenotypingAnalysis2.wdl`.
#
function CanonizeDipcallVcf() {
    local SAMPLE_ID=$1
    local INPUT_VCF_GZ=$2
    local MIN_SV_LENGTH=$3
    local MAX_SV_LENGTH=$4
    local CHROMOSOME=$5
    
    # Limiting to the desired chromosome
    ${TIME_COMMAND} bcftools view --output-type z ${INPUT_VCF_GZ} ${CHROMOSOME} > ${SAMPLE_ID}_out.vcf.gz
    mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    # Splitting multiallelic records into biallelic records
    ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    # Removing SNVs, records that are not marked as present, and records with
    # unresolved REF/ALT.
    ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || REF="*" || ALT="*"' --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    # Making sure SVLEN and SVTYPE are consistently annotated
    truvari anno svinfo --minsize 1 ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    # Keeping only records in the given length range
    ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_canonized.vcf.gz
    mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_canonized.vcf.gz.tbi
}


function CountDipcallVcf() {
    local SAMPLE_ID=$1
    local INPUT_VCF_GZ=$2
    local INPUT_BED=$3
    
    # Whole chromosome
    N_RECORDS_TOTAL=$(bcftools index --nrecords ${INPUT_VCF_GZ}.tbi)
    ${TIME_COMMAND} bcftools filter --exclude 'FILTER!="PASS" && FILTER!="."' --output-type z ${INPUT_VCF_GZ} > ${SAMPLE_ID}_out.vcf.gz
    tabix -f ${SAMPLE_ID}_out.vcf.gz
    N_RECORDS_NO_FILTER=$(bcftools index --nrecords ${SAMPLE_ID}_out.vcf.gz.tbi)
    rm -f ${SAMPLE_ID}_out.vcf.gz*
    
    # Dipcall BED
    ${TIME_COMMAND} bcftools filter --regions-file ${INPUT_BED} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > ${SAMPLE_ID}_bed.vcf.gz
    tabix -f ${SAMPLE_ID}_bed.vcf.gz
    N_RECORDS_BED_TOTAL=$(bcftools index --nrecords ${SAMPLE_ID}_bed.vcf.gz.tbi)
    ${TIME_COMMAND} bcftools filter --exclude 'FILTER!="PASS" && FILTER!="."' --output-type z ${SAMPLE_ID}_bed.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
    tabix -f ${SAMPLE_ID}_out.vcf.gz
    N_RECORDS_BED_NO_FILTER=$(bcftools index --nrecords ${SAMPLE_ID}_out.vcf.gz.tbi)
    rm -f ${SAMPLE_ID}_out.vcf.gz* ${SAMPLE_ID}_bed.vcf.gz*
    
    echo "${N_RECORDS_TOTAL},${N_RECORDS_NO_FILTER},${N_RECORDS_BED_TOTAL},${N_RECORDS_BED_NO_FILTER}" > ${SAMPLE_ID}_dipcall.csv
}




# ------------------------------- Main program ---------------------------------

TIME_COMMAND=" "
tabix -f ${SAMPLE_ID}.vcf.gz
CanonizeDipcallVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz ${MIN_SV_LENGTH} ${MAX_SV_LENGTH} ${CHROMOSOME}
CountDipcallVcf ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}.bed
