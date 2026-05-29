#!/bin/bash
#
REMOTE_INDIR_BEFORE="gs://fc-7f861a33-ddb4-4b2f-8d10-5679c9df6108/v3/15x/workpackage_1"
REMOTE_INDIR_AFTER="gs://fc-secure-95bbd6eb-6d63-49aa-a980-47f3c1342b1e/scratch/cunial_intersample_vcf/v3/ultralong_annotate_production/second_attempt"
REMOTE_INDIR_TRUTH="gs://fc-7f861a33-ddb4-4b2f-8d10-5679c9df6108/v3/training_resource_for_ultralong_svimasm/individual_vcfs"
SAMPLE_IDS_FILE=""

set -euxo pipefail

INFINITY="1000000000"

rm -f output.csv
while read SAMPLE_ID; do
    # Localizing truth
    gcloud storage cp ${REMOTE_INDIR_TRUTH}/${SAMPLE_ID}_canonized.vcf.gz ./${SAMPLE_ID}_truth.vcf.gz
    gcloud storage cp ${REMOTE_INDIR_TRUTH}/${SAMPLE_ID}_canonized.vcf.gz.tbi ./${SAMPLE_ID}_truth.vcf.gz.tbi

    # Localizing the "before" VCF. Working only on SVTYPEs that are not INS, for 
    # simplicity.
    gcloud storage cp ${REMOTE_INDIR_BEFORE}/${SAMPLE_ID}_ultralong.bcf ./${SAMPLE_ID}_before.bcf
    gcloud storage cp ${REMOTE_INDIR_BEFORE}/${SAMPLE_ID}_ultralong.bcf.csi ./${SAMPLE_ID}_before.bcf.csi
    bcftools filter --include 'SVTYPE=="DEL" || SVTYPE=="DUP" || SVTYPE=="INV"' --output-type z ./${SAMPLE_ID}_before.bcf --output ${SAMPLE_ID}_before.vcf.gz
    rm -f ${SAMPLE_ID}_before.bcf* ; bcftools index -f -t ${SAMPLE_ID}_before.vcf.gz

    # Localizing the "after" VCF. Working only on SVTYPEs that are neither INS
    # nor INSDUP, for simplicity.
    gcloud storage cp ${REMOTE_INDIR_AFTER}/${SAMPLE_ID}_del.vcf.'gz*' .
    gcloud storage cp ${REMOTE_INDIR_AFTER}/${SAMPLE_ID}_dup.vcf.'gz*' .
    gcloud storage cp ${REMOTE_INDIR_AFTER}/${SAMPLE_ID}_inv.vcf.'gz*' .
    bcftools concat --allow-overlaps --remove-duplicates --output-type z ./${SAMPLE_ID}_del.vcf.gz ./${SAMPLE_ID}_dup.vcf.gz ./${SAMPLE_ID}_inv.vcf.gz --output ${SAMPLE_ID}_after.vcf.gz
    rm -f ${SAMPLE_ID}_del.vcf.gz* ${SAMPLE_ID}_dup.vcf.gz* ${SAMPLE_ID}_inv.vcf.gz* ; bcftools index -f -t ${SAMPLE_ID}_after.vcf.gz

    # Computing GT concordance
    truvari bench -b ${SAMPLE_ID}_truth.vcf.gz -c ${SAMPLE_ID}_before.vcf.gz --sizemin 10000 --sizemax ${INFINITY} --sizefilt 1 --refdist 500 --pctseq 0 --pctsize 0.9 --pctovl 0 --pick single -o ./${SAMPLE_ID}_truvari/
    CONCORDANCE_BEFORE=$(grep 'gt_concordance' ./${SAMPLE_ID}_truvari/summary.json | cut -w -f 3)
    rm -rf ./${SAMPLE_ID}_truvari/
    truvari bench -b ${SAMPLE_ID}_truth.vcf.gz -c ${SAMPLE_ID}_after.vcf.gz --sizemin 10000 --sizemax ${INFINITY} --sizefilt 1 --refdist 500 --pctseq 0 --pctsize 0.9 --pctovl 0 --pick single -o ./${SAMPLE_ID}_truvari/
    CONCORDANCE_AFTER=$(grep 'gt_concordance' ./${SAMPLE_ID}_truvari/summary.json | cut -w -f 3)
    rm -rf ./${SAMPLE_ID}_truvari/

    # Next iteration
    echo -e "${SAMPLE_ID},${CONCORDANCE_BEFORE},${CONCORDANCE_AFTER}" >> output.csv
done < ${SAMPLE_IDS_FILE}
