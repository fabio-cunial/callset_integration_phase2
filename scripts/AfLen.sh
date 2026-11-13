#!/bin/bash
#
# Draws (allele frequency, SV length) heatmaps for INS and DEL inside and
# outside tandem repeats.
#
REMOTE_BUCKET="..."
REMOTE_TR_BED="${REMOTE_BUCKET}/scratch/cunial_intersample_vcf/v2/plot_hwe/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed"
MIN_OVERLAP_FRACTION="0.9"
N_THREADS="2"
X_BINS="100"
Y_BINS="50"
LOG_X="1"

set -euxo pipefail

gcloud storage cp ${REMOTE_TR_BED} ./tr.bed
N_SAMPLES=$(bcftools query --list-samples v3_main.bcf | wc -l)
bcftools +fill-tags --threads ${N_THREADS} --output-type u v3_main.bcf -- -t AC | bcftools view --threads ${N_THREADS} --drop-genotypes --output-type z --output annotated.vcf.gz
bedtools coverage -a annotated.vcf.gz -b tr.bed > coverage.tsv
(bcftools view --header-only annotated.vcf.gz; awk -F'\t' -v OFS='\t' -v min_overlap_fraction=${MIN_OVERLAP_FRACTION} '/^#/ { print; next } $NF >= min_overlap_fraction { NF -= 4; print }' coverage.tsv) | bcftools query --format '%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/AC\n' > inside_tr.tsv
(bcftools view --header-only annotated.vcf.gz; awk -F'\t' -v OFS='\t' -v min_overlap_fraction=${MIN_OVERLAP_FRACTION} '/^#/ { print; next } $NF < min_overlap_fraction { NF -= 4; print }' coverage.tsv)  | bcftools query --format '%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/AC\n' > outside_tr.tsv
python AfLen.py inside_tr.tsv outside_tr.tsv aflen.png ${N_SAMPLES} ${X_BINS} ${Y_BINS} ${LOG_X}
