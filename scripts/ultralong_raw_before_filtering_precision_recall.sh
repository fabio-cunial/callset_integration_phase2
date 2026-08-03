#!/bin/bash
#
COVERAGE="30x"
BUCKET_URI=" "
ROOT_DIR_QUERY="${BUCKET_URI}/v3/${COVERAGE}/workpackage_1"
ROOT_DIR_TRUTH="${BUCKET_URI}/v3/training_resource_for_ultralong_svimasm/ultralong_truth"

set -euxo pipefail


while read SAMPLE_ID; do
    gcloud storage cp ${ROOT_DIR_QUERY}/${SAMPLE_ID}_ultralong.bcf ./query.bcf
    gcloud storage cp ${ROOT_DIR_TRUTH}/${SAMPLE_ID}_svimasm.vcf.gz ./truth.vcf.gz
    
    bcftools index query.bcf
    bcftools index --nrecords query.bcf
    
    bcftools filter --include 'SVTYPE="DEL"' --output-type z query.bcf --output query_del.vcf.gz
    tabix -f query_del.vcf.gz
    bcftools filter --include 'SVTYPE="DEL" && ABS(SVLEN)>10000' --output-type z truth.vcf.gz --output truth_del.vcf.gz
    tabix -f truth_del.vcf.gz
    
    rm -rf ./truvari/
    truvari bench -c query_del.vcf.gz -b truth_del.vcf.gz --pctseq 0 --sizemin 1 --sizefilt 1 --sizemax 2000000 -o ./truvari/
    
    bcftools filter --include 'SVTYPE="INS"' --output-type z query.bcf --output query_ins.vcf.gz
    tabix -f query_ins.vcf.gz
    bcftools filter --include 'SVTYPE="INS" && ABS(SVLEN)>10000' --output-type z truth.vcf.gz --output truth_ins.vcf.gz
    tabix -f truth_ins.vcf.gz
    
    rm -rf ./truvari/
    truvari bench -c query_ins.vcf.gz -b truth_ins.vcf.gz --pctseq 0 --sizemin 1 --sizefilt 1 --sizemax 2000000 -o ./truvari/
done < samples.txt
