#!/bin/bash
#
SAMPLE_ID=""
INPUT_VCF_GZ=""
INPUT_BAM=""
INPUT_FAI=""

N_THREADS="4"
CLASSPATH="."
N_BINS="10"
BREAKPOINT_WINDOW_BP="500"

TIME_COMMAND="/usr/bin/time --verbose"


# --------------------------- Steps of the pipeline ----------------------------

# Annotates an input VCF with coverage measures extracted from a BAM.
#
function AnnotateCoverageBins() {
    SAMPLE_ID=$1
    INPUT_VCF=$2
    INPUT_BAM=$3
    N_BINS=$4
    BREAKPOINT_WINDOW_BP=$5

    ${TIME_COMMAND} java -cp ${CLASSPATH} UltralongIntervalGetBins ${INPUT_VCF} ${INPUT_FAI} ${N_BINS} ${BREAKPOINT_WINDOW_BP} > ${SAMPLE_ID}_bins.bed
    ${TIME_COMMAND} samtools bedcov ${SAMPLE_ID}_bins.bed ${INPUT_BAM} > ${SAMPLE_ID}_counts.bed
    rm -f ${SAMPLE_ID}_bins.bed
    ${TIME_COMMAND} java -cp ${CLASSPATH} UltralongIntervalAnnotateCoverage ${SAMPLE_ID}_counts.bed $(( ${N_BINS} + 4 )) | sort -k 1,1 > ${SAMPLE_ID}_tags.tsv
    rm -f ${SAMPLE_ID}_counts.bed
    ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\n' ${INPUT_VCF} | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
    ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_tags.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
    rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_tags.tsv
    tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
    echo '##INFO=<ID=BIN_BEFORE_COVERAGE,Number=1,Type=Float,Description="Coverage of the bin before the call">' > ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=BIN_LEFT_COVERAGE,Number=1,Type=Float,Description="Coverage of the left-breakpoint bin">' >> ${SAMPLE_ID}_header.txt
    COLUMNS='CHROM,POS,~ID,INFO/BIN_BEFORE_COVERAGE,INFO/BIN_LEFT_COVERAGE'
    for i in $( seq 1 ${N_BINS} ); do
        echo '##INFO=<ID=BIN_'${i}'_COVERAGE,Number=1,Type=Float,Description="Coverage of the i-th bin">' >> ${SAMPLE_ID}_header.txt
        COLUMNS=${COLUMNS}',INFO/BIN_'${i}'_COVERAGE'
    done
    echo '##INFO=<ID=BIN_RIGHT_COVERAGE,Number=1,Type=Float,Description="Coverage of the right-breakpoint bin">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=BIN_AFTER_COVERAGE,Number=1,Type=Float,Description="Coverage of the bin after the call">' >> ${SAMPLE_ID}_header.txt
    COLUMNS=${COLUMNS}',INFO/BIN_RIGHT_COVERAGE,INFO/BIN_AFTER_COVERAGE'
    ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
    rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
}


# Annotates an input VCF with clipped alignment measures extracted from a BAM.
#
function AnnotateClippedAlignments() {
    SAMPLE_ID=$1
    INPUT_VCF=$2
    INPUT_BAM=$3
    BREAKPOINT_WINDOW_BP=$4
    
    ${TIME_COMMAND} java -cp ${CLASSPATH} UltralongIntervalGetBins ${INPUT_VCF} ${INPUT_FAI} 0 ${BREAKPOINT_WINDOW_BP} > ${SAMPLE_ID}_bins.bed
    while read CHROM START END ID; do
        ${TIME_COMMAND} samtools view --threads ${N_THREADS} ${INPUT_BAM} ${CHROM}:${START}-${END} > ${ID}.sam
        java -cp ${CLASSPATH} UltralongIntervalGetClips ${ID}.sam ${START} ${END} ${ID}
        rm -f ${ID}.sam
        sort ${ID}_left.clips > ${ID}_left_sorted.clips
        sort ${ID}_right.clips > ${ID}_right_sorted.clips
        rm -f ${ID}_left.clips ${ID}_right.clips
    done < ${SAMPLE_ID}_bins.bed
    rm -f ${SAMPLE_ID}_variants.txt
    ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} | sort | uniq > ${SAMPLE_ID}_variantID_sorted.txt
    rm -f ${SAMPLE_ID}_counts.tsv
    while read ID; do
        LL=$(wc -l < ${ID}_left_left_sorted.clips)
        LR=$(wc -l < ${ID}_left_right_sorted.clips)
        RL=$(wc -l < ${ID}_right_left_sorted.clips)
        RR=$(wc -l < ${ID}_right_right_sorted.clips)
        LL_RL=$(comm -1 -2 ${ID}_left_left_sorted.clips ${ID}_right_left_sorted.clips | wc -l)
        LL_RR=$(comm -1 -2 ${ID}_left_left_sorted.clips ${ID}_right_right_sorted.clips | wc -l)
        LR_RL=$(comm -1 -2 ${ID}_left_right_sorted.clips ${ID}_right_left_sorted.clips | wc -l)
        LR_RR=$(comm -1 -2 ${ID}_left_right_sorted.clips ${ID}_right_right_sorted.clips | wc -l)
        echo "${ID}\t${LL}\t${LR}\t${RL}\t${RR}\t${LL_RL}\t${LL_RR}\t${LR_RL}\t${LR_RR}" >> ${SAMPLE_ID}_counts.tsv
        rm -f ${ID}_*.clips
    done < ${SAMPLE_ID}_variantID_sorted.txt
    rm -f ${SAMPLE_ID}_variantID_sorted.txt
    ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-3 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
    ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
    rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv
    tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
    echo '##INFO=<ID=LL,Number=1,Type=Integer,Description="Left window: number of left-clipped alignemnts.">' > ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR,Number=1,Type=Integer,Description="Left window: number of right-clipped alignemnts.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=RL,Number=1,Type=Integer,Description="Right window: number of left-clipped alignemnts.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=RR,Number=1,Type=Integer,Description="Right window: number of right-clipped alignemnts.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LL_RL,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LL_RR,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR_RL,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR_RR,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    COLUMNS='CHROM,POS,~ID,INFO/LL,INFO/LR,INFO/RL,INFO/RR,INFO/LL_RL,INFO/LL_RR,INFO/LR_RL,INFO/LR_RR'
    ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
    rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
}




# ------------------------------ Main program ----------------------------------

cp ${INPUT_VCF_GZ} in.vcf.gz
cp ${INPUT_VCF_GZ}.tbi in.vcf.gz.tbi

# Annotating DELs
bcftools filter --threads ${N_THREADS} --include 'SVTYPE="DEL"' --output-type v in.vcf.gz --output out.vcf
rm -f in.vcf.gz* ; mv out.vcf in.vcf

AnnotateCoverageBins ${SAMPLE_ID} in.vcf ${INPUT_BAM} ${N_BINS} ${BREAKPOINT_WINDOW_BP}
rm -f in.vcf ; mv ${SAMPLE_ID}_annotated.vcf in.vcf
AnnotateClippedAlignments ${SAMPLE_ID} in.vcf ${INPUT_BAM} ${BREAKPOINT_WINDOW_BP}
rm -f in.vcf ; mv ${SAMPLE_ID}_annotated.vcf in.vcf

bcftools view --threads ${N_THREADS} --output-type z in.vcf --output ${SAMPLE_ID}_del_annotated.vcf.gz
bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_del_annotated.vcf.gz
