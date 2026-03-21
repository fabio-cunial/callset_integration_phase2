#!/bin/bash
#

N_THREADS=$3
INPUT_FAI=$4






# Annotates an input VCF with coverage statistics from a BAM.
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
    ${TIME_COMMAND} java -cp ${CLASSPATH} UltralongIntervalAnnotateCoverage ${SAMPLE_ID}_counts.bed ${N_BINS} | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
    tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
    echo '##INFO=<ID=BIN_MINUS_1_COVERAGE,Number=1,Type=Float,Description="Coverage of the i-th bin">' > ${SAMPLE_ID}_header.txt
    COLUMNS='CHROM,POS,~ID,INFO/BIN_MINUS_1_COVERAGE'
    for i in {1..${N_BINS}}; do
        echo '##INFO=<ID=BIN_'${i}'_COVERAGE,Number=1,Type=Float,Description="Coverage of the i-th bin">' >> ${SAMPLE_ID}_header.txt
        COLUMNS=${COLUMNS}',INFO/BIN_'${i}'_COVERAGE'
    done
    ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
    rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
}


function AnnotateClippedAlignments() {
    
    
    ${TIME_COMMAND} java -cp ${CLASSPATH} UltralongIntervalGetBins ${INPUT_VCF} ${INPUT_FAI} 0 ${BREAKPOINT_WINDOW_BP} > ${SAMPLE_ID}_bins.bed
    while read CHROM START END ID; do
        ${TIME_COMMAND} samtools view --threads ${N_THREADS} --uncompressed ${INPUT_BAM} ${CHROM}:${START}-${END} > ${ID}.sam
        java -cp ${CLASSPATH} UltralongIntervalGetClips ${ID}.sam ${START} ${END} ${ID}
        rm -f ${ID}.sam
        sort ${ID}_left.clips > ${ID}_left_sorted.clips
        sort ${ID}_right.clips > ${ID}_right_sorted.clips
        rm -f ${ID}_left.clips ${ID}_right.clips
    done < ${SAMPLE_ID}_bins.bed
    rm -f ${SAMPLE_ID}_variants.txt
    while read CHROM START END ID; do
        if [[ ${ID%_left} != ${ID} ]]; then
            echo ${ID%_left} >> ${SAMPLE_ID}_variants.txt
        elif [[ ${ID%_right} != ${ID} ]]; then
            echo ${ID%_right} >> ${SAMPLE_ID}_variants.txt
        fi
    done < ${SAMPLE_ID}_bins.bed
    sort ${SAMPLE_ID}_variants.txt | uniq > ${SAMPLE_ID}_variants_sorted.txt
    rm -f ${SAMPLE_ID}_variants.txt
    rm -f annotations.tsv
    while read ID; do
        LL=$(wc -l < ${ID}_left_left_sorted.clips)
        LR=$(wc -l < ${ID}_left_right_sorted.clips)
        RL=$(wc -l < ${ID}_right_left_sorted.clips)
        RR=$(wc -l < ${ID}_right_right_sorted.clips)
        LL_RL=$(comm -1 -2 ${ID}_left_left_sorted.clips ${ID}_right_left_sorted.clips | wc -l)
        LL_RR=$(comm -1 -2 ${ID}_left_left_sorted.clips ${ID}_right_right_sorted.clips | wc -l)
        LR_RL=$(comm -1 -2 ${ID}_left_right_sorted.clips ${ID}_right_left_sorted.clips | wc -l)
        LR_RR=$(comm -1 -2 ${ID}_left_right_sorted.clips ${ID}_right_right_sorted.clips | wc -l)
        echo "${ID}\t${LL}\t${LR}\t${RL}\t${RR}\t${LL_RL}\t${LL_RR}\t${LR_RL}\t${LR_RR}" >> counts.tsv
        rm -f ${ID}_*.clips
    done < ${SAMPLE_ID}_variants_sorted.txt
    ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-3 | sort -k 3,3 > chrom_pos_id.tsv
    join -t $'\t' -1 3 -2 1 chrom_pos_id.tsv counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11); }' > annotations.tsv
    rm -f chrom_pos_id.tsv counts.tsv
    
    
    
    
    
    
}






bcftools filter --threads ${N_THREADS} --include 'SVTYPE=DEL' --output-type v in.vcf.gz --output out.vcf
mv out.vcf in.vcf
