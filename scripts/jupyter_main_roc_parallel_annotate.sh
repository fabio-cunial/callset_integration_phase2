#!/bin/bash
#
# Script for the jupyter notebook that plots ROC curves for the main VCF.
# Subsets to a confident BED and then annotates TRs and TPs using dipcall truth.
#
N_THREADS=$1  # Should be equal to the number of physical cores
SAMPLE_IDS_TXT="sample_ids.txt"
BGZIP_COMMAND="/home/jupyter/HPRC_Y2/edit/htslib-1.19.1/bgzip"
TABIX_COMMAND="/home/jupyter/HPRC_Y2/edit/htslib-1.19.1/tabix"
BCFTOOLS_COMMAND="/home/jupyter/HPRC_Y2/edit/bcftools-1.19/bcftools"
TRUVARI_COMMAND="truvari"
TR_BED="human_GRCh38_no_alt_analysis_set.trf.bed"
TR_BED_GZ="${TR_BED}.gz"
CONFIDENT_BED="GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed"
INFINITY="1000000000"

set -euo pipefail


cat << 'END' > annotate.sh
#!/bin/bash
#
SAMPLE_ID=$1

set -euo pipefail

# Subsetting to the sample's dipcall BED
${BCFTOOLS_COMMAND} view --regions-file dipcall-beds/${SAMPLE_ID}.bed --output-type z --regions-overlap pos scored-vcfs/${SAMPLE_ID}_scored.bcf --output scored-vcfs/${SAMPLE_ID}_scored_dipcallbed.vcf.gz

# Annotating query records with confident region status
${BCFTOOLS_COMMAND} annotate --no-version -a ${CONFIDENT_BED} -h confident_header.txt \
    -m +CONFIDENT_ONE -c CHROM,FROM,TO --min-overlap :0.9 scored-vcfs/${SAMPLE_ID}_scored_dipcallbed.vcf.gz | \
    grep -v '##INFO=<ID=CONFIDENT_ONE' | \
    sed -E 's/CONFIDENT_ONE/IN_CONFIDENT_BED=1/g' | \
${BCFTOOLS_COMMAND} annotate --no-version -a ${CONFIDENT_BED} \
    -m -CONFIDENT_ZERO -c CHROM,FROM,TO --min-overlap :0.9 | \
    grep -v '##INFO=<ID=CONFIDENT_ZERO' | \
    sed -E 's/CONFIDENT_ZERO/IN_CONFIDENT_BED=0/g' | \
${BCFTOOLS_COMMAND} view -Oz -o scored-vcfs/${SAMPLE_ID}_confident.vcf.gz
${BCFTOOLS_COMMAND} index -f -t scored-vcfs/${SAMPLE_ID}_confident.vcf.gz

# Annotating query records with TR status
${BCFTOOLS_COMMAND} annotate --no-version -a ${TR_BED_GZ} -h tr_header.txt \
    -m +TR_ONE -c CHROM,FROM,TO --min-overlap :0.9 scored-vcfs/${SAMPLE_ID}_confident.vcf.gz | \
    grep -v '##INFO=<ID=TR_ONE' | \
    sed -E 's/TR_ONE/TR=1/g' | \
${BCFTOOLS_COMMAND} annotate --no-version -a ${TR_BED_GZ} \
    -m -TR_ZERO -c CHROM,FROM,TO --min-overlap :0.9 | \
    grep -v '##INFO=<ID=TR_ZERO' | \
    sed -E 's/TR_ZERO/TR=0/g' | \
${BCFTOOLS_COMMAND} view -Oz -o scored-vcfs/${SAMPLE_ID}_confident_tr.vcf.gz
${BCFTOOLS_COMMAND} index -f -t scored-vcfs/${SAMPLE_ID}_confident_tr.vcf.gz
rm -f scored-vcfs/${SAMPLE_ID}_confident.vcf.gz*

# Annotating query records with truth status, using `truvari bench`.
${BCFTOOLS_COMMAND} index -f -t truth-vcfs/${SAMPLE_ID}_canonized.vcf.gz
rm -rf ./scored-vcfs/${SAMPLE_ID}_truvari/
${TRUVARI_COMMAND} bench -b truth-vcfs/${SAMPLE_ID}_canonized.vcf.gz -c scored-vcfs/${SAMPLE_ID}_confident_tr.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick single -o ./scored-vcfs/${SAMPLE_ID}_truvari/
mv ./scored-vcfs/${SAMPLE_ID}_truvari/summary.json ./scored-vcfs/${SAMPLE_ID}_truvari.json
echo '##INFO=<ID=MATCHES_DIPCALL,Number=1,Type=Integer,Description="Matches dipcall">' > scored-vcfs/${SAMPLE_ID}_header.txt
# TPs
${BCFTOOLS_COMMAND} query --format '%CHROM\t%POS\t%REF\t%ALT\t1\n' ./scored-vcfs/${SAMPLE_ID}_truvari/tp-comp.vcf.gz | ${BGZIP_COMMAND} -c > scored-vcfs/${SAMPLE_ID}_tps.tsv.gz
${TABIX_COMMAND} -f -s1 -b2 -e2 scored-vcfs/${SAMPLE_ID}_tps.tsv.gz
${BCFTOOLS_COMMAND} annotate --annotations scored-vcfs/${SAMPLE_ID}_tps.tsv.gz --header-lines scored-vcfs/${SAMPLE_ID}_header.txt --columns CHROM,POS,REF,ALT,MATCHES_DIPCALL --output-type z scored-vcfs/${SAMPLE_ID}_confident_tr.vcf.gz --output scored-vcfs/${SAMPLE_ID}_confident_tr_tp.vcf.gz
${BCFTOOLS_COMMAND} index -f -t scored-vcfs/${SAMPLE_ID}_confident_tr_tp.vcf.gz
# FPs
${BCFTOOLS_COMMAND} query --format '%CHROM\t%POS\t%REF\t%ALT\t0\n' ./scored-vcfs/${SAMPLE_ID}_truvari/fp.vcf.gz | ${BGZIP_COMMAND} -c > scored-vcfs/${SAMPLE_ID}_fps.tsv.gz
${TABIX_COMMAND} -f -s1 -b2 -e2 scored-vcfs/${SAMPLE_ID}_fps.tsv.gz
${BCFTOOLS_COMMAND} annotate --annotations scored-vcfs/${SAMPLE_ID}_fps.tsv.gz --columns CHROM,POS,REF,ALT,MATCHES_DIPCALL --output-type z scored-vcfs/${SAMPLE_ID}_confident_tr_tp.vcf.gz --output scored-vcfs/${SAMPLE_ID}_confident_tr_tp_fp.vcf.gz
${BCFTOOLS_COMMAND} index -f -t scored-vcfs/${SAMPLE_ID}_confident_tr_tp_fp.vcf.gz
rm -rf scored-vcfs/${SAMPLE_ID}_tps.tsv.gz* scored-vcfs/${SAMPLE_ID}_fps.tsv.gz* scored-vcfs/${SAMPLE_ID}_confident_tr.vcf.gz* scored-vcfs/${SAMPLE_ID}_confident_tr_tp.vcf.gz* scored-vcfs/${SAMPLE_ID}_header.txt scored-vcfs/${SAMPLE_ID}_truvari/
END
chmod +x annotate.sh




# -------------------------------- Main program --------------------------------

export N_THREADS BGZIP_COMMAND TABIX_COMMAND BCFTOOLS_COMMAND TRUVARI_COMMAND TR_BED_GZ CONFIDENT_BED SAMPLE_IDS_TXT INFINITY
echo '##INFO=<ID=TR,Number=1,Type=Integer,Description=">=90% of variant interval overlaps with tandem repeat region">' > tr_header.txt
echo '##INFO=<ID=IN_CONFIDENT_BED,Number=1,Type=Integer,Description=">=90% of variant interval overlaps with confident region">' > confident_header.txt
cut -f 1-3 ${TR_BED} | ${BGZIP_COMMAND} -c > ${TR_BED_GZ}
${TABIX_COMMAND} -f -0 -s1 -b2 -e3 ${TR_BED_GZ}
xargs --arg-file=${SAMPLE_IDS_TXT} --max-lines=1 --max-procs=${N_THREADS} ./annotate.sh
rm -f tr_header.txt ${TR_BED_GZ}*
