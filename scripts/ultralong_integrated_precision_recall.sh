#!/bin/bash
#
COVERAGE="15x"
COHORT_BCF="ultralong_lenient.bcf"
BUCKET_URI=" "
ROOT_DIR_TRUTH="${BUCKET_URI}/v3/training_resource_for_ultralong_svimasm/ultralong_truth"
N_THREADS="2"

set -euxo pipefail


while read SAMPLE_ID; do
    # Preparing the query
    bcftools view --threads ${N_THREADS} --samples ${SAMPLE_ID} --output-type u ${COHORT_BCF} | bcftools +fill-tags --threads ${N_THREADS} --output-type u - -- --tags AC | bcftools view --threads ${N_THREADS} --min-ac 1 --output-type z --write-index --output query.vcf.gz
    bcftools index --threads ${N_THREADS} -t query.vcf.gz
    bcftools index --nrecords query.vcf.gz

    # Preparing the truth
    gcloud storage cp ${ROOT_DIR_TRUTH}/${SAMPLE_ID}_svimasm.vcf.gz ./truth.vcf.gz
    
    bcftools filter --threads ${N_THREADS} --include 'SVTYPE="DEL"' --output-type z query.vcf.gz --output query_del.vcf.gz
    tabix -f query_del.vcf.gz
    bcftools filter --threads ${N_THREADS} --include 'SVTYPE="DEL" && ABS(SVLEN)>10000' --output-type z truth.vcf.gz --output truth_del.vcf.gz
    tabix -f truth_del.vcf.gz
    
    rm -rf ./truvari/
    truvari bench -c query_del.vcf.gz -b truth_del.vcf.gz --pctseq 0 --sizemin 1 --sizefilt 1 --sizemax 2000000 -o ./truvari/
    
    bcftools filter --threads ${N_THREADS} --include 'SVTYPE="INS"' --output-type z query.vcf.gz --output query_ins.vcf.gz
    tabix -f query_ins.vcf.gz
    bcftools filter --threads ${N_THREADS} --include 'SVTYPE="INS" && ABS(SVLEN)>10000' --output-type z truth.vcf.gz --output truth_ins.vcf.gz
    tabix -f truth_ins.vcf.gz
    
    rm -rf ./truvari/
    truvari bench -c query_ins.vcf.gz -b truth_ins.vcf.gz --pctseq 0 --sizemin 1 --sizefilt 1 --sizemax 2000000 -o ./truvari/
done < samples_${COVERAGE}.txt
