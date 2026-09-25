#!/bin/bash
#
# Script for the jupyter notebook that plots ROC curves for the main VCF.
# Subsets to a confident BED and measures precision and recall using dipcall 
# truth.
#
N_THREADS=$1  # Should be equal to the number of physical cores
SAMPLE_IDS_TXT="sample_ids.txt"
BGZIP_COMMAND="/home/jupyter/HPRC_Y2/edit/htslib-1.19.1/bgzip"
TABIX_COMMAND="/home/jupyter/HPRC_Y2/edit/htslib-1.19.1/tabix"
BCFTOOLS_COMMAND="/home/jupyter/HPRC_Y2/edit/bcftools-1.19/bcftools"
TRUVARI_COMMAND="truvari"
INFINITY="1000000000"

set -euo pipefail


# Remark: dipcall VCFs are assumed to contain only INS and DEL, to have already 
# been subset to dipcall's BED, and to only contain calls in [20..10k] bp.
#
# Remark: query VCFs are assumed to only contain calls in [20..10k] bp, and to
# have had their DUPs been converted to INS. They may contain other types, e.g.
# INV.
#
cat << 'END' > bench.sh
#!/bin/bash
SAMPLE_ID=$1

set -euo pipefail

# Subsetting the truth to:
# - standard chromosomes.
${BCFTOOLS_COMMAND} index -f -t truth-vcfs/${SAMPLE_ID}_canonized.vcf.gz
${BCFTOOLS_COMMAND} view --regions $(echo chr{1..22} | tr ' ' ','),chrX,chrY --output-type z truth-vcfs/${SAMPLE_ID}_canonized.vcf.gz --output truth-vcfs/${SAMPLE_ID}_truth.vcf.gz
${BCFTOOLS_COMMAND} index -f -t truth-vcfs/${SAMPLE_ID}_truth.vcf.gz

# Subsetting the query VCFs to:
# - INS/DEL only;
# - the sample's dipcall BED;
# - standard chromosomes.
${BCFTOOLS_COMMAND} view --include 'SVTYPE="INS" || SVTYPE="DEL"' --regions-file dipcall-beds/${SAMPLE_ID}.bed --regions-overlap pos --output-type u merged-vcfs/${SAMPLE_ID}_sv.vcf.gz --output-type z --output merged-vcfs/${SAMPLE_ID}_sv_dipcallbed.vcf.gz
${BCFTOOLS_COMMAND} index -f -t merged-vcfs/${SAMPLE_ID}_sv_dipcallbed.vcf.gz
${BCFTOOLS_COMMAND} view --regions $(echo chr{1..22} | tr ' ' ','),chrX,chrY --output-type z merged-vcfs/${SAMPLE_ID}_sv_dipcallbed.vcf.gz --output merged-vcfs/${SAMPLE_ID}_sv_dipcallbed_chr.vcf.gz
${BCFTOOLS_COMMAND} index -f -t merged-vcfs/${SAMPLE_ID}_sv_dipcallbed_chr.vcf.gz

${BCFTOOLS_COMMAND} view --include 'SVTYPE="INS" || SVTYPE="DEL"' --regions-file dipcall-beds/${SAMPLE_ID}.bed --regions-overlap pos --output-type z scored-vcfs/${SAMPLE_ID}_scored.bcf --output scored-vcfs/${SAMPLE_ID}_scored_dipcallbed.vcf.gz
${BCFTOOLS_COMMAND} index -f -t scored-vcfs/${SAMPLE_ID}_scored_dipcallbed.vcf.gz
${BCFTOOLS_COMMAND} view --regions $(echo chr{1..22} | tr ' ' ','),chrX,chrY --output-type z scored-vcfs/${SAMPLE_ID}_scored_dipcallbed.vcf.gz --output scored-vcfs/${SAMPLE_ID}_scored_dipcallbed_chr.vcf.gz
${BCFTOOLS_COMMAND} index -f -t scored-vcfs/${SAMPLE_ID}_scored_dipcallbed_chr.vcf.gz

# Benchmarking: before kanpig.
rm -rf ./merged-vcfs/${SAMPLE_ID}_truvari/
${TRUVARI_COMMAND} bench -b truth-vcfs/${SAMPLE_ID}_truth.vcf.gz -c merged-vcfs/${SAMPLE_ID}_sv_dipcallbed_chr.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick single -o ./merged-vcfs/${SAMPLE_ID}_truvari/
mv ./merged-vcfs/${SAMPLE_ID}_truvari/summary.json ./merged-vcfs/${SAMPLE_ID}_truvari_all.json
rm -rf ./merged-vcfs/${SAMPLE_ID}_truvari/
${TRUVARI_COMMAND} bench -b truth-vcfs/${SAMPLE_ID}_truth.vcf.gz -c merged-vcfs/${SAMPLE_ID}_sv_dipcallbed_chr.vcf.gz --sizemin 50 --sizemax ${INFINITY} --sizefilt 50 --pctsize 0.9 --pctseq 0.9 --pick single -o ./merged-vcfs/${SAMPLE_ID}_truvari/
mv ./merged-vcfs/${SAMPLE_ID}_truvari/summary.json ./merged-vcfs/${SAMPLE_ID}_truvari_50.json
rm -rf ./merged-vcfs/${SAMPLE_ID}_truvari/
${TRUVARI_COMMAND} bench -b truth-vcfs/${SAMPLE_ID}_truth.vcf.gz -c merged-vcfs/${SAMPLE_ID}_sv_dipcallbed_chr.vcf.gz --sizemin 1 --sizemax 49 --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick single -o ./merged-vcfs/${SAMPLE_ID}_truvari/
mv ./merged-vcfs/${SAMPLE_ID}_truvari/summary.json ./merged-vcfs/${SAMPLE_ID}_truvari_20.json
rm -rf ./merged-vcfs/${SAMPLE_ID}_truvari/

# Benchmarking: after kanpig.
rm -rf ./scored-vcfs/${SAMPLE_ID}_truvari/
${TRUVARI_COMMAND} bench -b truth-vcfs/${SAMPLE_ID}_truth.vcf.gz -c scored-vcfs/${SAMPLE_ID}_scored_dipcallbed_chr.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick single -o ./scored-vcfs/${SAMPLE_ID}_truvari/
mv ./scored-vcfs/${SAMPLE_ID}_truvari/summary.json ./scored-vcfs/${SAMPLE_ID}_truvari_all.json
rm -rf ./scored-vcfs/${SAMPLE_ID}_truvari/
${TRUVARI_COMMAND} bench -b truth-vcfs/${SAMPLE_ID}_truth.vcf.gz -c scored-vcfs/${SAMPLE_ID}_scored_dipcallbed_chr.vcf.gz --sizemin 50 --sizemax ${INFINITY} --sizefilt 50 --pctsize 0.9 --pctseq 0.9 --pick single -o ./scored-vcfs/${SAMPLE_ID}_truvari/
mv ./scored-vcfs/${SAMPLE_ID}_truvari/summary.json ./scored-vcfs/${SAMPLE_ID}_truvari_50.json
rm -rf ./scored-vcfs/${SAMPLE_ID}_truvari/
${TRUVARI_COMMAND} bench -b truth-vcfs/${SAMPLE_ID}_truth.vcf.gz -c scored-vcfs/${SAMPLE_ID}_scored_dipcallbed_chr.vcf.gz --sizemin 1 --sizemax 49 --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick single -o ./scored-vcfs/${SAMPLE_ID}_truvari/
mv ./scored-vcfs/${SAMPLE_ID}_truvari/summary.json ./scored-vcfs/${SAMPLE_ID}_truvari_20.json
rm -rf ./scored-vcfs/${SAMPLE_ID}_truvari/
END
chmod +x bench.sh




# -------------------------------- Main program --------------------------------

export N_THREADS BGZIP_COMMAND TABIX_COMMAND BCFTOOLS_COMMAND TRUVARI_COMMAND SAMPLE_IDS_TXT INFINITY
xargs --arg-file=${SAMPLE_IDS_TXT} --max-lines=1 --max-procs=${N_THREADS} ./bench.sh

# Outputting table
BEFORE_ALL_P=$( cat merged-vcfs/*_truvari_all.json | grep precision | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')
BEFORE_ALL_R=$( cat merged-vcfs/*_truvari_all.json | grep recall | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')
BEFORE_ALL_C=$( cat merged-vcfs/*_truvari_all.json | grep concordance | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')

AFTER_ALL_P=$( cat scored-vcfs/*_truvari_all.json | grep precision | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')
AFTER_ALL_R=$( cat scored-vcfs/*_truvari_all.json | grep recall | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')
AFTER_ALL_C=$( cat scored-vcfs/*_truvari_all.json | grep concordance | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')

BEFORE_50_P=$( cat merged-vcfs/*_truvari_50.json | grep precision | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')
BEFORE_50_R=$( cat merged-vcfs/*_truvari_50.json | grep recall | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')
BEFORE_50_C=$( cat merged-vcfs/*_truvari_50.json | grep concordance | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')

AFTER_50_P=$( cat scored-vcfs/*_truvari_50.json | grep precision | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')
AFTER_50_R=$( cat scored-vcfs/*_truvari_50.json | grep recall | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')
AFTER_50_C=$( cat scored-vcfs/*_truvari_50.json | grep concordance | awk '{print $2}' | tr ',' ' ' | awk '{sum+=$1} END {print sum/NR}')

echo -e "Before kanpig, all:\tprecision=${BEFORE_ALL_P}\trecall=${BEFORE_ALL_R}\tgtc=${BEFORE_ALL_C}"
echo -e "After kanpig, all:\tprecision=${AFTER_ALL_P}\trecall=${AFTER_ALL_R}\tgtc=${AFTER_ALL_C}"
echo -e "Before kanpig, >=50:\tprecision=${BEFORE_50_P}\trecall=${BEFORE_50_R}\tgtc=${BEFORE_50_C}"
echo -e "After kanpig, >=50:\tprecision=${AFTER_50_P}\trecall=${AFTER_50_R}\tgtc=${AFTER_50_C}"
