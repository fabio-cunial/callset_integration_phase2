#!/bin/bash
#
INPUT_VCF_GZ="truvari_collapsed.bcf"

REPEATMASKER_TSV_GZ="repeatmasker.gz"
CENTROMERES_BED="centromeres.bed"
GAPS_TSV="gaps.tsv"
SEGDUPS_BED="GRCh38_segdups.bed"
CONFIDENT_BED="GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed"

REFERENCE_FAI="GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"

MERES_SLACK_BP="2000"
SEGDUPS_SLACK_BP="2000"
REPEAT_SLACK_BP="200"

set -euxo pipefail

# Building the allowed regions BED
zgrep -E 'LINE|SINE|LTR|Retroposon|Simple_repeat|Satellite|DNA' ${REPEATMASKER_TSV_GZ} | cut -f 6,7,8 | grep -vE '_alt|_fix' | bedtools slop -i - -g ${REFERENCE_FAI} -b ${REPEAT_SLACK_BP} | bedtools sort -i - -g ${REFERENCE_FAI} | bedtools merge -i - | bedtools complement -i - -g ${REFERENCE_FAI} | bedtools sort -i - -g ${REFERENCE_FAI} > good_repeatmasker_sorted.bed
java GetHardfilterBed ${REFERENCE_FAI} ${CENTROMERES_BED} ${GAPS_TSV} ${MERES_SLACK_BP} good_meres1.bed 0.01 good_meres2.bed
bedtools sort -i good_meres1.bed -g ${REFERENCE_FAI} > good_meres1_sorted.bed
bedtools slop -i ${SEGDUPS_BED} -g ${REFERENCE_FAI} -b ${SEGDUPS_SLACK_BP} | bedtools sort -i - -g ${REFERENCE_FAI} | bedtools merge -i - | bedtools complement -i - -g ${REFERENCE_FAI} | bedtools sort -i - -g ${REFERENCE_FAI} > good_segdups_sorted.bed
bedtools intersect -a good_repeatmasker_sorted.bed -b good_meres1_sorted.bed -b good_segdups_sorted.bed -b ${CONFIDENT_BED} > good.bed
rm -f good_repeatmasker_sorted.bed good_meres1.bed good_meres2.bed good_meres1_sorted.bed good_segdups_sorted.bed

# Filtering BNDs
bcftools view --drop-genotypes --output-type v ${INPUT_VCF_GZ} > cleaned.vcf
java BndCanonize cleaned.vcf | bcftools sort --output-type z > left.vcf.gz
tabix -f left.vcf.gz
rm -f cleaned.vcf
java BndSymmetrize left.vcf.gz | bcftools sort --output-type z > right.vcf.gz
tabix -f right.vcf.gz
bcftools view --no-header --regions-file good.bed --regions-overlap pos --output-type v left.vcf.gz | sort -t $'\t' -k 3,3 > left_good.tsv
bcftools view --no-header --regions-file good.bed --regions-overlap pos --output-type v right.vcf.gz | sort -t $'\t' -k 3,3 > right_good.tsv
bcftools view --header-only left.vcf.gz > result.vcf
join -t $'\t' -1 3 -2 3 left_good.tsv right_good.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$2,$3,$1,$4,$5,$6,$7,$8); }' | sort -t $'\t' -k 1,1 -k2,2n >> result.vcf
rm -f left_good.tsv right_good.tsv left.vcf.gz right.vcf.gz
