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
ADJACENCY_SLACK_BP="200"

TIME_COMMAND="/usr/bin/time --verbose"


# --------------------------- Steps of the pipeline ----------------------------

# Given an input VCF containing only interval records, the procedure annotates
# each record with coverage measures extracted from a BAM.
#
# Remark: `samtools bedcov` skips reads with any of the following flags set:
# UNMAP, SECONDARY, QCFAIL, DUP
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
    ${TIME_COMMAND} java -cp ${CLASSPATH} UltralongIntervalCreateBedcovAnnotations ${SAMPLE_ID}_counts.bed $(( ${N_BINS} + 4 )) | sort -k 1,1 > ${SAMPLE_ID}_tags.tsv
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


    cat << 'END' > annotate_mapq_secondary.sh
#!/bin/bash
set -euxo pipefail

INPUT_BAM=$1
CHROM=$2
START=$3
END=$4
ID=$5

samtools view ${INPUT_BAM} ${CHROM}:${START}-${END} | awk '{ sum+=$5; count++ } END { print (count>0?sum/count:0) }' > ${ID}_mapq.txt
samtools view --count ${INPUT_BAM} --require-flags 256 ${CHROM}:${START}-${END} > ${ID}_secondary.txt
END
    chmod +x annotate_mapq_secondary.sh


# Given an input VCF containing only interval records, the procedure annotates 
# each record with the average MAPQ and the number of secondary alignments
# (=repeat-induced multi-mappings) over its left and right breakpoint, extracted
# from a BAM.
#
# Remark: we do not collect the number of supplementary alignments, since more
# specific counts are already captured by `AnnotateClippedAlignments()`.
#
function AnnotateMapqSecondary() {
    SAMPLE_ID=$1
    INPUT_VCF=$2
    INPUT_BAM=$3
    BREAKPOINT_WINDOW_BP=$4
    
    ${TIME_COMMAND} java -cp ${CLASSPATH} UltralongIntervalGetBins ${INPUT_VCF} ${INPUT_FAI} 0 ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
    ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_mapq_secondary.sh ${INPUT_BAM}
    rm -f ${SAMPLE_ID}_bins.wsv
    ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} | sort | uniq > ${SAMPLE_ID}_variantID_sorted.txt
    rm -f ${SAMPLE_ID}_counts.tsv
    while read ID; do
        MAPQ_LEFT=$(cat ${ID}_left_mapq.txt)
        MAPQ_RIGHT=$(cat ${ID}_right_mapq.txt)
        SECONDARY_LEFT=$(cat ${ID}_left_secondary.txt)
        SECONDARY_RIGHT=$(cat ${ID}_right_secondary.txt)
        echo -e "${ID}\t${MAPQ_LEFT}\t${MAPQ_RIGHT}\t${SECONDARY_LEFT}\t${SECONDARY_RIGHT}" >> ${SAMPLE_ID}_counts.tsv
    done < ${SAMPLE_ID}_variantID_sorted.txt
    rm -f ${SAMPLE_ID}_variantID_sorted.txt ${ID}_*_mapq.txt ${ID}_*_secondary.txt
    ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-3 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
    ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%f\t%f\t%d\t%d\n",$2,$3,$1,$4,$5,$6,$7); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
    rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv
    tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
    echo '##INFO=<ID=BIN_LEFT_MAPQ,Number=1,Type=Float,Description="Left breakpoint window: avg MAPQ.">' > ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=BIN_RIGHT_MAPQ,Number=1,Type=Float,Description="Right breakpoint window: avg MAPQ.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=BIN_LEFT_SECONDARY,Number=1,Type=Integer,Description="Left breakpoint window: number of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=BIN_RIGHT_SECONDARY,Number=1,Type=Integer,Description="Right breakpoint window: number of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
    COLUMNS='CHROM,POS,~ID,INFO/BIN_LEFT_MAPQ,INFO/BIN_RIGHT_MAPQ,INFO/BIN_LEFT_SECONDARY,INFO/BIN_RIGHT_SECONDARY'
    ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
    rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
}


    cat << 'END' > annotate_clipped_alignments_1.sh
#!/bin/bash
set -euxo pipefail

INPUT_BAM=$1
CLASSPATH=$2
CHROM=$3
START=$4
END=$5
ID=$6

samtools view ${INPUT_BAM} ${CHROM}:${START}-${END} > ${ID}.sam
java -cp ${CLASSPATH} UltralongIntervalGetClips ${ID}.sam ${START} ${END} ${ID}
rm -f ${ID}.sam
sort -k 1,1 ${ID}_leftmaximal.txt > ${ID}_leftmaximal_sorted.txt
sort -k 1,1 ${ID}_rightmaximal.txt > ${ID}_rightmaximal_sorted.txt
rm -f ${ID}_leftmaximal.txt ${ID}_rightmaximal.txt
END
    chmod +x annotate_clipped_alignments_1.sh
    
    
    cat << 'END' > annotate_clipped_alignments_2.sh
#!/bin/bash
set -euxo pipefail

CLASSPATH=$1
ADJACENCY_SLACK_BP=$2
ID=$3

LL=$(wc -l < ${ID}_left_leftmaximal_sorted.txt)
LR=$(wc -l < ${ID}_left_rightmaximal_sorted.txt)
RL=$(wc -l < ${ID}_right_leftmaximal_sorted.txt)
RR=$(wc -l < ${ID}_right_rightmaximal_sorted.txt)
LL_RL=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_leftmaximal_sorted.txt ${LL} 1 ${ID}_right_leftmaximal_sorted.txt ${RL} 1 ${ADJACENCY_SLACK_BP} | tr ',' '\t')
LL_RR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_leftmaximal_sorted.txt ${LL} 1 ${ID}_right_rightmaximal_sorted.txt ${RR} 0 ${ADJACENCY_SLACK_BP} | tr ',' '\t')
LR_RL=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_rightmaximal_sorted.txt ${LR} 0 ${ID}_right_leftmaximal_sorted.txt ${RL} 1 ${ADJACENCY_SLACK_BP} | tr ',' '\t')
LR_RR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_rightmaximal_sorted.txt ${LR} 0 ${ID}_right_rightmaximal_sorted.txt ${RR} 0 ${ADJACENCY_SLACK_BP} | tr ',' '\t')
echo -e "${ID}\t${LL}\t${LR}\t${RL}\t${RR}\t${LL_RL}\t${LL_RR}\t${LR_RL}\t${LR_RR}" >> ${ID}_counts.txt
rm -f ${ID}_*maximal_sorted.txt
END
    chmod +x annotate_clipped_alignments_2.sh


# Given an input VCF containing only interval calls, the procedure annotates it
# with clipped alignment measures extracted from a BAM.
#
function AnnotateClippedAlignments() {
    SAMPLE_ID=$1
    INPUT_VCF=$2
    INPUT_BAM=$3
    BREAKPOINT_WINDOW_BP=$4
    ADJACENCY_SLACK_BP=$5
    
    ${TIME_COMMAND} java -cp ${CLASSPATH} UltralongIntervalGetBins ${INPUT_VCF} ${INPUT_FAI} 0 ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
    ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_1.sh ${INPUT_BAM} ${CLASSPATH}
    rm -f ${SAMPLE_ID}_bins.wsv
    ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} > ${SAMPLE_ID}_variantID.txt
    rm -f *_counts.tsv
    ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_variantID.txt --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_2.sh ${CLASSPATH} ${ADJACENCY_SLACK_BP}
    rm -f ${SAMPLE_ID}_variantID.txt
    cat *_counts.txt | sort -k 1,1 > ${SAMPLE_ID}_counts.tsv
    rm -f *_counts.txt
    ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-3 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
    ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
    rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv
    tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
    echo '##INFO=<ID=LL,Number=1,Type=Integer,Description="Left window: number of left-clipped alignments.">' > ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR,Number=1,Type=Integer,Description="Left window: number of right-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=RL,Number=1,Type=Integer,Description="Right window: number of left-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=RR,Number=1,Type=Integer,Description="Right window: number of right-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LL_RL_1,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LL_RL_2,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LL_RL_3,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LL_RL_4,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LL_RR_1,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LL_RR_2,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LL_RR_3,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LL_RR_4,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR_RL_1,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR_RL_2,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR_RL_3,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR_RL_4,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left window and a left-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR_RR_1,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR_RR_2,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR_RR_3,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    echo '##INFO=<ID=LR_RR_4,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment in the left window and a right-clipped alignment in the right window.">' >> ${SAMPLE_ID}_header.txt
    COLUMNS='CHROM,POS,~ID,INFO/LL,INFO/LR,INFO/RL,INFO/RR,INFO/LL_RL_1,INFO/LL_RL_2,INFO/LL_RL_3,INFO/LL_RL_4,INFO/LL_RR_1,INFO/LL_RR_2,INFO/LL_RR_3,INFO/LL_RR_4,INFO/LR_RL_1,INFO/LR_RL_2,INFO/LR_RL_3,INFO/LR_RL_4,INFO/LR_RR_1,INFO/LR_RR_2,INFO/LR_RR_3,INFO/LR_RR_4'
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
AnnotateMapqSecondary ${SAMPLE_ID} in.vcf ${INPUT_BAM} ${BREAKPOINT_WINDOW_BP}
rm -f in.vcf ; mv ${SAMPLE_ID}_annotated.vcf in.vcf
AnnotateClippedAlignments ${SAMPLE_ID} in.vcf ${INPUT_BAM} ${BREAKPOINT_WINDOW_BP} ${ADJACENCY_SLACK_BP}
rm -f in.vcf ; mv ${SAMPLE_ID}_annotated.vcf in.vcf

bcftools view --threads ${N_THREADS} --output-type z in.vcf --output ${SAMPLE_ID}_del_annotated.vcf.gz
bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_del_annotated.vcf.gz
