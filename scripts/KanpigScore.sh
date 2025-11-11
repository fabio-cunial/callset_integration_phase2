#!/bin/bash
#
MIN_SV_LENGTH="20"
MAX_SV_LENGTH="10000"




function Download() {
    # gsutil -m cp gs://fc-secure-95bbd6eb-6d63-49aa-a980-47f3c1342b1e/scratch/cunial_intersample_vcf/v2/regenotyping_analysis_v2_50bp/128_samples/kanpig/HG00097_kanpig.vcf.'gz*' .
#     gsutil cp gs://fc-7f861a33-ddb4-4b2f-8d10-5679c9df6108/submissions/0a354214-9225-4416-97b5-73a86057ad06/runDipcall/31377972-e48a-4785-83f1-6635edf5f615/call-dipcall/attempt-2/glob-5b09b846642a9bb3e1843523e167602a/HG00097_hap2_hprc_r2_v1.0.1.dipcall.vcf.gz ./HG00097_dipcall.vcf.gz
#     tabix HG00097_dipcall.vcf.gz
#     gsutil cp gs://fc-7f861a33-ddb4-4b2f-8d10-5679c9df6108/submissions/0a354214-9225-4416-97b5-73a86057ad06/runDipcall/31377972-e48a-4785-83f1-6635edf5f615/call-dipcall/attempt-2/glob-6a71e0052edde14f09e7720a93686506/HG00097_hap2_hprc_r2_v1.0.1.dipcall.bed ./HG00097_dipcall.bed
#     gsutil cp gs://fc-secure-95bbd6eb-6d63-49aa-a980-47f3c1342b1e/submissions/28a86f40-dd3b-407c-b3c3-45ea0d785756/SV_Integration_RegenotypingAnalysis/ec034476-86b6-4fe9-ad3b-5334f2945c3a/call-ComplementBed/cacheCopy/sorted.bed ./tr.bed
#     gsutil cp gs://fc-secure-95bbd6eb-6d63-49aa-a980-47f3c1342b1e/submissions/28a86f40-dd3b-407c-b3c3-45ea0d785756/SV_Integration_RegenotypingAnalysis/ec034476-86b6-4fe9-ad3b-5334f2945c3a/call-ComplementBed/cacheCopy/complement.bed ./not_tr.bed
    
    SAMPLE_ID="HG00097"
    SAMPLE_VCF_GZ="HG00097_kanpig.vcf.gz"
    SAMPLE_TBI="HG00097_kanpig.vcf.gz.tbi"
    DIPCALL_VCF_GZ="HG00097_dipcall.vcf.gz"
    DIPCALL_BED="HG00097_dipcall.bed"
    TANDEM_BED="tr.bed"
    NOT_TANDEM_BED="not_tr.bed"
}


# Puts in canonical form a raw VCF from dipcall. This is identical to
# `SV_Integration_BuildTrainingResource.wdl`.
#
function CanonizeDipcallVcf() {
    local SAMPLE_ID=$1
    local INPUT_VCF_GZ=$2
    local INPUT_TBI=$3
    local INPUT_BED=$4
    
    cp ${INPUT_VCF_GZ} ${SAMPLE_ID}_in.vcf.gz
    cp ${INPUT_TBI} ${SAMPLE_ID}_in.vcf.gz.tbi
    
    # Splitting multiallelic records into biallelic records
    ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    # Removing SNVs, records that are not marked as present, records
    # with a FILTER, and records with unresolved REF/ALT.
    ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || COUNT(GT="alt")=0 || (FILTER!="PASS" && FILTER!=".") || REF="*" || ALT="*"' --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    # Keeping only records in the dipcall BED
    ${TIME_COMMAND} bcftools filter --regions-file ${INPUT_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    # Making sure SVLEN and SVTYPE are consistently annotated
    truvari anno svinfo --minsize 1 ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    # Keeping only records in the given length range
    ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_truth.vcf.gz
    mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_truth.vcf.gz.tbi
}


# Puts in canonical form an inter-sample-regenotyped kanpig VCF.
#
function CanonizeKanpigVcf() {
    local SAMPLE_ID=$1
    local INPUT_VCF_GZ=$2
    local INPUT_TBI=$3
    local DIPCALL_BED=$4
    
    cp ${INPUT_VCF_GZ} ${SAMPLE_ID}_in.vcf.gz
    cp ${INPUT_TBI} ${SAMPLE_ID}_in.vcf.gz.tbi
    
    # Keeping only records that are genotyped as present. This is important,
    # since truvari bench does not consider GTs when matching.
    ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    # Keeping only records in the dipcall BED
    ${TIME_COMMAND} bcftools filter --regions-file ${DIPCALL_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    # Keeping only records in the given length range
    ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
    rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
    
    mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_kanpig_canonized.vcf.gz
    mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_kanpig_canonized.vcf.gz.tbi
}


# Adds IS_POSITIVE and IS_TR to the INFO column of a canonized kanpig VCF.
#
function Annotate() {
    local SAMPLE_ID=$1
    local TANDEM_BED=$2
    local NOT_TANDEM_BED=$3
    
    rm -rf ./${SAMPLE_ID}_truvari_*
    ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_truth.vcf.gz -c ${SAMPLE_ID}_kanpig_canonized.vcf.gz --sizemin 1 --sizefilt 1 --sizemax ${INFINITY} -o ./${SAMPLE_ID}_truvari/
    
    # Annotating positives
    ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t1\n' ./${SAMPLE_ID}_truvari/tp-comp.vcf.gz | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
    tabix -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
    echo '##INFO=<ID=IS_POSITIVE,Number=1,Type=Integer,Description="According to truvari bench">' > ${SAMPLE_ID}_header.txt
    ${TIME_COMMAND} bcftools annotate --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns CHROM,POS,REF,ALT,INFO/IS_POSITIVE --output-type z ${SAMPLE_ID}_kanpig_canonized.vcf.gz > ${SAMPLE_ID}_kanpig_tmp1.vcf.gz
    tabix -f ${SAMPLE_ID}_kanpig_tmp1.vcf.gz
    rm -f ${SAMPLE_ID}_annotations.tsv.gz
    
    # Annotating negatives
    ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t0\n' ./${SAMPLE_ID}_truvari/fp.vcf.gz | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
    tabix -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
    ${TIME_COMMAND} bcftools annotate --annotations ${SAMPLE_ID}_annotations.tsv.gz --columns CHROM,POS,REF,ALT,INFO/IS_POSITIVE --output-type z ${SAMPLE_ID}_kanpig_tmp1.vcf.gz > ${SAMPLE_ID}_kanpig_tmp2.vcf.gz
    tabix -f ${SAMPLE_ID}_kanpig_tmp2.vcf.gz
    rm -f ${SAMPLE_ID}_annotations.tsv.gz
    
    # Annotating TRs
    ${TIME_COMMAND} bcftools view --regions-file ${TANDEM_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_kanpig_tmp2.vcf.gz > ${SAMPLE_ID}_kanpig_tr.vcf.gz
    ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t1\n' ${SAMPLE_ID}_kanpig_tr.vcf.gz | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
    tabix -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
    rm -f ${SAMPLE_ID}_kanpig_tr.vcf.gz
    echo '##INFO=<ID=IS_TR,Number=1,Type=Integer,Description="In TR BED or not">' > ${SAMPLE_ID}_header.txt
    ${TIME_COMMAND} bcftools annotate --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns CHROM,POS,REF,ALT,INFO/IS_TR --output-type z ${SAMPLE_ID}_kanpig_tmp2.vcf.gz > ${SAMPLE_ID}_kanpig_tmp3.vcf.gz
    tabix -f ${SAMPLE_ID}_kanpig_tmp3.vcf.gz
    rm -f ${SAMPLE_ID}_annotations.tsv.gz
    
    # Annotating non-TRs
    ${TIME_COMMAND} bcftools view --regions-file ${NOT_TANDEM_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_kanpig_tmp3.vcf.gz > ${SAMPLE_ID}_kanpig_not_tr.vcf.gz
    ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t0\n' ${SAMPLE_ID}_kanpig_not_tr.vcf.gz | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
    tabix -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
    rm -f ${SAMPLE_ID}_kanpig_not_tr.vcf.gz
    ${TIME_COMMAND} bcftools annotate --annotations ${SAMPLE_ID}_annotations.tsv.gz --columns CHROM,POS,REF,ALT,INFO/IS_TR --output-type z ${SAMPLE_ID}_kanpig_tmp3.vcf.gz > ${SAMPLE_ID}_kanpig_tmp4.vcf.gz
    tabix -f ${SAMPLE_ID}_kanpig_tmp4.vcf.gz
    rm -f ${SAMPLE_ID}_annotations.tsv.gz
    
    mv ${SAMPLE_ID}_kanpig_tmp4.vcf.gz ${SAMPLE_ID}_kanpig_annotated.vcf.gz
    mv ${SAMPLE_ID}_kanpig_tmp4.vcf.gz.tbi ${SAMPLE_ID}_kanpig_annotated.vcf.gz.tbi
    rm -rf ./${SAMPLE_ID}_truvari/ ${SAMPLE_ID}_kanpig_tmp*
}




# --------------------------- Main program -----------------------------
set -euxo pipefail

TIME_COMMAND=" "
INFINITY="1000000000"

Download 
CanonizeDipcallVcf ${SAMPLE_ID} ${DIPCALL_VCF_GZ} ${DIPCALL_VCF_GZ}.tbi ${DIPCALL_BED}
CanonizeKanpigVcf ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz ${SAMPLE_ID}_kanpig.vcf.gz.tbi ${DIPCALL_BED}
Annotate ${SAMPLE_ID} ${TANDEM_BED} ${NOT_TANDEM_BED}
