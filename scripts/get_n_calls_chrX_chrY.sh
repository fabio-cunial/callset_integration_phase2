#!/bin/bash
#
REMOTE_DIR="gs://..."
set -x

rm -f n_chrx_f.txt
while read SAMPLE_ID; do
    gcloud storage cp ${REMOTE_DIR}/${SAMPLE_ID}_kanpig.vcf.'gz*' .
    bcftools view --no-header ${SAMPLE_ID}_kanpig.vcf.gz chrX | wc -l >> n_chrx_f.txt
done < female_samples.txt

rm -f n_chrx_m.txt
rm -f n_chry_m.txt
while read SAMPLE_ID; do
    gcloud storage cp ${REMOTE_DIR}/${SAMPLE_ID}_kanpig.vcf.'gz*' .
    bcftools view --no-header ${SAMPLE_ID}_kanpig.vcf.gz chrX | wc -l >> n_chrx_m.txt
    bcftools view --no-header ${SAMPLE_ID}_kanpig.vcf.gz chrY | wc -l >> n_chry_m.txt
done < male_samples.txt
