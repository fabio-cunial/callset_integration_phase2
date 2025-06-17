#!/bin/bash
#
# The AGP file for GRCh38 was downloaded from:
# ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.agp.gz
# Remark: all coordinates in AGP are one-based and inclusive, see e.g.
# https://www.ncbi.nlm.nih.gov/genbank/genome_agp_specification/
#
AGP_FILE="/Users/fcunial/Downloads/hg38.agp"
REF_FAI="/Users/fcunial/Downloads/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"

set -euxo pipefail

rm -f gaps.bed intervals.txt chunk_*

# Extracting all gaps of unknown length or of length >=500 in the canonical
# chromosomes.
awk '{ if (length($1)<=5 && (($5=="N" && $6>=500) || $5=="U")) { printf("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2-1,$3-1+1,$4,$5,$6,$7,$8,$9); } }' ${AGP_FILE} | sort -V > gaps.bed
bedtools complement -L -i gaps.bed -g ${REF_FAI} | awk '{ printf("%s:%s-%s\n",$1,$2,$3); }' > intervals.txt
split -l 1 -d -a 3 intervals.txt chunk_
