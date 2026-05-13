version 1.0


# 
#
workflow SV_Integration_UltralongAnnotate {
    input {
        File chunk_csv
        String remote_outdir
        
        File reference_fa
        File reference_fai
        
        Int ins2dup_bin_length = 100
        Float ins2dup_bin_coverage_ratio = 1.5
        
        Int custom_n_coverage_bins = 10
        Int custom_breakpoint_window_bp = 500
        Int custom_min_clip_length = 200
        Int custom_adjacency_slack_bp = 300

        File feature_extraction_py

        File tr_bed
        File segdup_bed
        Float repeat_overlap_fraction = 0.8
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
        Int preemptible_number = 5
    }
    parameter_meta {
        chunk_csv: "Format: ID,bai,bam,csi,bcf"
    }
    
    call Impl {
        input:
            chunk_csv = chunk_csv,
            remote_outdir = remote_outdir,
            
            reference_fa = reference_fa,
            reference_fai = reference_fai,

            ins2dup_bin_length = ins2dup_bin_length,
            ins2dup_bin_coverage_ratio = ins2dup_bin_coverage_ratio,
            
            custom_n_coverage_bins = custom_n_coverage_bins,
            custom_breakpoint_window_bp = custom_breakpoint_window_bp,
            custom_min_clip_length = custom_min_clip_length,
            custom_adjacency_slack_bp = custom_adjacency_slack_bp,
    
            feature_extraction_py = feature_extraction_py,

            tr_bed = tr_bed,
            segdup_bed = segdup_bed,
            repeat_overlap_fraction = repeat_overlap_fraction,

            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    
    output {
    }
}


# Performance on a 4-core, 8GB VM.
#
# TOOL                                                CPU     RAM     TIME
# BAM download                                                          5m
#
# samtools bedcov                                     80%    900M      15m
# annotate_mapq_secondary.sh                         400%     15M      10s
# bcftools annotate                                  100%     15M      50s
# java UltralongIntervalGetBins                      200%     50M       1s
# java UltralongIntervalCreateBedcovAnnotations      200%     50M       1s
# annotate_clipped_alignments_1.sh                   300%    500M       2m
# annotate_clipped_alignments_2.sh                   300%     50M       2m
# cutefc (1 thread)                                   30%    1.5G      50m
# cutefc (2 threads)                                  50%    1.5G      25m
# samtools view (for extracted BAM, 2 threads)       100%     20M      10m
# cutefc (2 threads, extracted BAM)                  200%    900M       2m               
# feature_extraction.py                              100%     12G       7m
#
# UltralongInsGetIntervals                           200%     50M       1m
# xargs interval_2_breakpoints.sh                    200%    1.5G       1m
# UltralongInsExtractDups                            200%     60M       1s
# truvari collapse                                   100%    100M       1m
#
task Impl {
    input {
        File chunk_csv
        String remote_outdir
        
        File reference_fa
        File reference_fai

        Int ins2dup_bin_length
        Float ins2dup_bin_coverage_ratio
        
        Int custom_n_coverage_bins
        Int custom_breakpoint_window_bp
        Int custom_min_clip_length
        Int custom_adjacency_slack_bp

        File feature_extraction_py

        File tr_bed
        File segdup_bed
        Float repeat_overlap_fraction
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 50
        Int preemptible_number
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        # @param 
        # $2 A row of `chunk_csv`.
        #
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            local ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 2)
            local ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 3)
            local ULTRALONG_CSI=$(echo ${LINE} | cut -d , -f 4)
            local ULTRALONG_BCF=$(echo ${LINE} | cut -d , -f 5)
            
            date 1>&2
            gcloud storage cp ${ALIGNED_BAM} ./${SAMPLE_ID}.bam
            date 1>&2
            gcloud storage cp ${ALIGNED_BAI} ./${SAMPLE_ID}.bam.bai
            gcloud storage cp ${ULTRALONG_BCF} ./${SAMPLE_ID}.bcf
            gcloud storage cp ${ULTRALONG_CSI} ./${SAMPLE_ID}.bcf.csi
            
            # Converting to .vcf.gz for downstream tools
            bcftools view --threads ${N_THREADS} --output-type z ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}.vcf.gz
            rm -f ${SAMPLE_ID}.bcf*
        }
        
        
        # Deletes all files related to the given sample
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}*.bam* ${SAMPLE_ID}*.bcf* ${SAMPLE_ID}*.vcf* ${SAMPLE_ID}*.tsv* ${SAMPLE_ID}*.csv*
        }
        
        
        # Ensures that the VCF is correctly formatted.
        #
        function CanonizeVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            local DEFAULT_QUAL="60"   # Arbitrary
            
            # 1. Removing SVLEN from symbolic ALTs and fixing END.
            # This is necessary for `truvari bench` to work on symbolic records:
            # https://github.com/ACEnglish/truvari/issues/290
            java -cp ~{docker_dir} FixUltralongRecords ${INPUT_VCF_GZ} ~{reference_fai} > ${SAMPLE_ID}_out.vcf
            rm -f ${INPUT_VCF_GZ}* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 2. Cleaning REF, ALT, QUAL, FILTER.
            # - REF and ALT must be uppercase for XGBoost scoring downstream to
            #   work.
            # - We force every record to PASS, to rule out any filter-dependent
            #   effect in downstream tools.
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanRefAltQual ${SAMPLE_ID}_in.vcf ${DEFAULT_QUAL} > ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 3. Making sure IDs will be distinct at the inter-sample level 
            # (they are already distinct at the intra-sample level, thanks to
            # the steps upstream).
            ${TIME_COMMAND} bcftools annotate --set-id ${SAMPLE_ID}'_%ID' --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            bcftools view --threads ${N_THREADS} --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_canonized.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_canonized.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf
        }
        
        
        
        
        # ------------------------- Custom annotations -------------------------

cat << 'END' > samtools_bedcov_thread.sh
#!/bin/bash
set -euxo pipefail

INPUT_BAM=$1
INPUT_BED=$2

samtools bedcov ${INPUT_BED} ${INPUT_BAM} > ${INPUT_BED}_counts.bed
END
        chmod +x samtools_bedcov_thread.sh


        # Given an input VCF containing only interval records, the procedure
        # annotates each record with coverage measures extracted from a BAM.
        #
        # Remark: `samtools bedcov` skips reads with any of the following flags
        # set: UNMAP, SECONDARY, QCFAIL, DUP
        #
        function AnnotateCoverageBins_Interval() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local N_BINS=$4
            local BREAKPOINT_WINDOW_BP=$5

            # Running `samtools bedcov` in parallel, since it can be very slow
            # for large intervals.
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalGetBins ${INPUT_VCF} ~{reference_fai} ${N_BINS} ${BREAKPOINT_WINDOW_BP} > ${SAMPLE_ID}_bins.bed
            split -d -a 5 -l 1 ${SAMPLE_ID}_bins.bed ${SAMPLE_ID}_chunk_
            ls ${SAMPLE_ID}_chunk_* | sort -V > ${SAMPLE_ID}_list.txt
            rm -f ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_list.txt --max-lines=1 --max-procs=${N_THREADS} ./samtools_bedcov_thread.sh ${INPUT_BAM}
            rm -f ${SAMPLE_ID}_counts.bed
            while read -u 4 CHUNK; do
                cat ${CHUNK}_counts.bed >> ${SAMPLE_ID}_counts.bed
            done 4< ${SAMPLE_ID}_list.txt
            rm -f ${SAMPLE_ID}_chunk_* ${SAMPLE_ID}_list.txt 
            
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalCreateBedcovAnnotations ${SAMPLE_ID}_counts.bed $(( ${N_BINS} + 4 )) | sort -k 1,1 > ${SAMPLE_ID}_tags.tsv
            rm -f ${SAMPLE_ID}_counts.bed
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\n' ${INPUT_VCF} | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_tags.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_tags.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_BEFORE_COVERAGE,Number=1,Type=Float,Description="Coverage of the bin before the call">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_LEFT_COVERAGE,Number=1,Type=Float,Description="Coverage of the left-breakpoint bin">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/BIN_BEFORE_COVERAGE,INFO/BIN_LEFT_COVERAGE'
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
        
        
        # Given an input VCF containing only point records, the procedure
        # annotates each record with its the coverage of a window around POS 
        # extracted from a BAM.
        #
        # Remark: `samtools bedcov` skips reads with any of the following flags
        # set: UNMAP, SECONDARY, QCFAIL, DUP
        #
        function AnnotateCoverageBins_Point() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4

            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongPointGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} > ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} samtools bedcov ${SAMPLE_ID}_bins.bed ${INPUT_BAM} > ${SAMPLE_ID}_counts.bed
            rm -f ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongPointCreateBedcovAnnotations ${SAMPLE_ID}_counts.bed ${BREAKPOINT_WINDOW_BP} | sort -k 1,1 > ${SAMPLE_ID}_tags.tsv
            rm -f ${SAMPLE_ID}_counts.bed
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\n' ${INPUT_VCF} | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_tags.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%d\n",$2,$3,$1,$4); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_tags.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_POS,Number=1,Type=Float,Description="Coverage of the bin around the POS breakpoint">' > ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/BIN_POS'
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

samtools view --no-header ${INPUT_BAM} ${CHROM}:${START}-${END} | awk '{ sum+=$5; count++ } END { print (count>0?sum/count:0) }' > ${ID}_mapq.txt
samtools view --count ${INPUT_BAM} --require-flags 256 ${CHROM}:${START}-${END} > ${ID}_secondary.txt
END
        chmod +x annotate_mapq_secondary.sh


        # Given an input VCF containing only interval records, the procedure
        # annotates each record with the average MAPQ and the number of
        # secondary alignments (=repeat-induced multi-mappings) over its left
        # and right breakpoint, extracted from a BAM.
        #
        # Remark: we do not collect the number of supplementary alignments,
        # since more specific counts are already captured by
        # `AnnotateClippedAlignments()`.
        #
        function AnnotateMapqSecondary_Interval() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
    
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalGetBins ${INPUT_VCF} ~{reference_fai} 0 ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_mapq_secondary.sh ${INPUT_BAM}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} | sort | uniq > ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_counts.tsv
            while read -u 5 ID; do
                local MAPQ_LEFT=$(cat ${ID}_left_mapq.txt)
                local MAPQ_RIGHT=$(cat ${ID}_right_mapq.txt)
                local SECONDARY_LEFT=$(cat ${ID}_left_secondary.txt)
                local SECONDARY_RIGHT=$(cat ${ID}_right_secondary.txt)
                echo -e "${ID}\t${MAPQ_LEFT}\t${MAPQ_RIGHT}\t${SECONDARY_LEFT}\t${SECONDARY_RIGHT}" >> ${SAMPLE_ID}_counts.tsv
                rm -f ${ID}_*_mapq.txt ${ID}_*_secondary.txt
            done 5< ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_variantID_sorted.txt
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-3 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%f\t%f\t%d\t%d\n",$2,$3,$1,$4,$5,$6,$7); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_LEFT_MAPQ,Number=1,Type=Float,Description="Left breakpoint window: avg MAPQ.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_RIGHT_MAPQ,Number=1,Type=Float,Description="Right breakpoint window: avg MAPQ.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_LEFT_SECONDARY,Number=1,Type=Integer,Description="Left breakpoint window: number of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_RIGHT_SECONDARY,Number=1,Type=Integer,Description="Right breakpoint window: number of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/BIN_LEFT_MAPQ,INFO/BIN_RIGHT_MAPQ,INFO/BIN_LEFT_SECONDARY,INFO/BIN_RIGHT_SECONDARY'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }
        
        
        # Given an input VCF containing only point records, the procedure
        # annotates each record with the average MAPQ and the number of
        # secondary alignments (=repeat-induced multi-mappings) over its POS
        # breakpoint, extracted from a BAM.
        #
        # Remark: we do not collect the number of supplementary alignments,
        # since more specific counts are already captured by
        # `AnnotateClippedAlignments()`.
        #
        function AnnotateMapqSecondary_Point() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
    
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongPointGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_mapq_secondary.sh ${INPUT_BAM}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} | sort | uniq > ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_counts.tsv
            while read -u 6 ID; do
                local MAPQ_POINT=$(cat ${ID}_point_mapq.txt)
                local SECONDARY_POINT=$(cat ${ID}_point_secondary.txt)
                echo -e "${ID}\t${MAPQ_POINT}\t${SECONDARY_POINT}" >> ${SAMPLE_ID}_counts.tsv
                rm -f ${ID}_*_mapq.txt ${ID}_*_secondary.txt
            done 6< ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_variantID_sorted.txt
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-3 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%f\t%d\n",$2,$3,$1,$4,$5); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_POINT_MAPQ,Number=1,Type=Float,Description="Breakpoint window: avg MAPQ.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_SECONDARY,Number=1,Type=Integer,Description="Breakpoint window: number of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/BIN_POINT_MAPQ,INFO/BIN_POINT_SECONDARY'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }

        cat << 'END' > annotate_clipped_alignments_1.sh
#!/bin/bash
set -euxo pipefail

INPUT_BAM=$1
CLASSPATH=$2
MIN_CLIP_LENGTH=$3
CHROM=$4
START=$5
END=$6
ID=$7

samtools view --no-header ${INPUT_BAM} ${CHROM}:${START}-${END} > ${ID}.sam
java -cp ${CLASSPATH} UltralongIntervalGetClips ${ID}.sam ${START} ${END} ${MIN_CLIP_LENGTH} ${ID}
rm -f ${ID}.sam
sort -k 1,1 ${ID}_leftmaximal.txt > ${ID}_leftmaximal_sorted.txt
sort -k 1,1 ${ID}_rightmaximal.txt > ${ID}_rightmaximal_sorted.txt
rm -f ${ID}_leftmaximal.txt ${ID}_rightmaximal.txt
END
        chmod +x annotate_clipped_alignments_1.sh
    
    
        cat << 'END' > annotate_clipped_alignments_2_interval.sh
#!/bin/bash
set -euxo pipefail

CLASSPATH=$1
ADJACENCY_SLACK_BP=$2
ID=$3

LL=$(wc -l < ${ID}_left_leftmaximal_sorted.txt)
LR=$(wc -l < ${ID}_left_rightmaximal_sorted.txt)
RL=$(wc -l < ${ID}_right_leftmaximal_sorted.txt)
RR=$(wc -l < ${ID}_right_rightmaximal_sorted.txt)
LL_RL=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_leftmaximal_sorted.txt ${LL} 1 ${ID}_right_leftmaximal_sorted.txt ${RL} 1 ${ADJACENCY_SLACK_BP} 0 | tr ',' '\t')
LL_RR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_leftmaximal_sorted.txt ${LL} 1 ${ID}_right_rightmaximal_sorted.txt ${RR} 0 ${ADJACENCY_SLACK_BP} 0 | tr ',' '\t')
LR_RL=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_rightmaximal_sorted.txt ${LR} 0 ${ID}_right_leftmaximal_sorted.txt ${RL} 1 ${ADJACENCY_SLACK_BP} 0 | tr ',' '\t')
LR_RR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_left_rightmaximal_sorted.txt ${LR} 0 ${ID}_right_rightmaximal_sorted.txt ${RR} 0 ${ADJACENCY_SLACK_BP} 0 | tr ',' '\t')
echo -e "${ID}\t${LL}\t${LR}\t${RL}\t${RR}\t${LL_RL}\t${LL_RR}\t${LR_RL}\t${LR_RR}" >> ${ID}_counts.txt
rm -f ${ID}_*maximal_sorted.txt
END
        chmod +x annotate_clipped_alignments_2_interval.sh

        
        cat << 'END' > annotate_clipped_alignments_2_point.sh
#!/bin/bash
set -euxo pipefail

CLASSPATH=$1
ADJACENCY_SLACK_BP=$2
ID=$3

PL=$(wc -l < ${ID}_point_leftmaximal_sorted.txt)
PR=$(wc -l < ${ID}_point_rightmaximal_sorted.txt)
PL_PL=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_point_leftmaximal_sorted.txt ${PL} 1 ${ID}_point_leftmaximal_sorted.txt ${PL} 1 ${ADJACENCY_SLACK_BP} 1 | tr ',' '\t')
PL_PR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_point_leftmaximal_sorted.txt ${PL} 1 ${ID}_point_rightmaximal_sorted.txt ${PR} 0 ${ADJACENCY_SLACK_BP} 1 | tr ',' '\t')
PR_PR=$(java -cp ${CLASSPATH} UltralongIntervalIntersectClips ${ID}_point_rightmaximal_sorted.txt ${PR} 0 ${ID}_point_rightmaximal_sorted.txt ${PR} 0 ${ADJACENCY_SLACK_BP} 1 | tr ',' '\t')
echo -e "${ID}\t${PL}\t${PR}\t${PL_PL}\t${PL_PR}\t${PR_PR}" >> ${ID}_counts.txt
rm -f ${ID}_*maximal_sorted.txt
END
        chmod +x annotate_clipped_alignments_2_point.sh


        # Given an input VCF containing only interval calls, the procedure
        # annotates it with clipped alignment measures extracted from a BAM.
        #
        function AnnotateClippedAlignments_Interval() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
            local ADJACENCY_SLACK_BP=$5
            local MIN_CLIP_LENGTH=$6

            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalGetBins ${INPUT_VCF} ~{reference_fai} 0 ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_1.sh ${INPUT_BAM} ~{docker_dir} ${MIN_CLIP_LENGTH}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} > ${SAMPLE_ID}_variantID.txt
            rm -f *_counts.tsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_variantID.txt --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_2_interval.sh ~{docker_dir} ${ADJACENCY_SLACK_BP}
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
            local COLUMNS='CHROM,POS,~ID,INFO/LL,INFO/LR,INFO/RL,INFO/RR,INFO/LL_RL_1,INFO/LL_RL_2,INFO/LL_RL_3,INFO/LL_RL_4,INFO/LL_RR_1,INFO/LL_RR_2,INFO/LL_RR_3,INFO/LL_RR_4,INFO/LR_RL_1,INFO/LR_RL_2,INFO/LR_RL_3,INFO/LR_RL_4,INFO/LR_RR_1,INFO/LR_RR_2,INFO/LR_RR_3,INFO/LR_RR_4'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }


        # Given an input VCF containing only point records, the procedure
        # annotates it with clipped alignment measures extracted from a BAM.
        #
        function AnnotateClippedAlignments_Point() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
            local ADJACENCY_SLACK_BP=$5
            local MIN_CLIP_LENGTH=$6
    
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongPointGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_1.sh ${INPUT_BAM} ~{docker_dir} ${MIN_CLIP_LENGTH}
            rm -f ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} bcftools query --format '%ID\n' ${INPUT_VCF} > ${SAMPLE_ID}_variantID.txt
            rm -f *_counts.tsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_variantID.txt --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_2_point.sh ~{docker_dir} ${ADJACENCY_SLACK_BP}
            rm -f ${SAMPLE_ID}_variantID.txt
            cat *_counts.txt | sort -k 1,1 > ${SAMPLE_ID}_counts.tsv
            rm -f *_counts.txt
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-3 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=PL,Number=1,Type=Integer,Description="Number of left-clipped alignments.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PR,Number=1,Type=Integer,Description="Number of right-clipped alignments.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PL_1,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment and a left-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PL_2,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment and a left-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PL_3,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment and a left-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PL_4,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment and a left-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PR_1,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PR_2,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PR_3,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PL_PR_4,Number=1,Type=Integer,Description="Number of reads with a left-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PR_PR_1,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PR_PR_2,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PR_PR_3,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=PR_PR_4,Number=1,Type=Integer,Description="Number of reads with a right-clipped alignment and a right-clipped alignment in the same window.">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/PL,INFO/PR,INFO/PL_PL_1,INFO/PL_PL_2,INFO/PL_PL_3,INFO/PL_PL_4,INFO/PL_PR_1,INFO/PL_PR_2,INFO/PL_PR_3,INFO/PL_PR_4,INFO/PR_PR_1,INFO/PR_PR_2,INFO/PR_PR_3,INFO/PR_PR_4'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }


        function AnnotateCustom_NotIns() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            
            mv ${INPUT_VCF} ${SAMPLE_ID}_in.vcf

            AnnotateCoverageBins_Interval ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_n_coverage_bins} ~{custom_breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            AnnotateMapqSecondary_Interval ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            AnnotateClippedAlignments_Interval ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp} ~{custom_adjacency_slack_bp} ~{custom_min_clip_length}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf

            mv ${SAMPLE_ID}_in.vcf ${INPUT_VCF}
        }

        
        function AnnotateCustom_Ins() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            
            mv ${INPUT_VCF} ${SAMPLE_ID}_in.vcf
            
            AnnotateCoverageBins_Point ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            AnnotateMapqSecondary_Point ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            AnnotateClippedAlignments_Point ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{custom_breakpoint_window_bp} ~{custom_adjacency_slack_bp} ~{custom_min_clip_length}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf

            mv ${SAMPLE_ID}_in.vcf ${INPUT_VCF}
        }
        
        
        
        
        # ------------------- Annotations from Kalra et al. --------------------
        
        # Runs a slightly modified version of the code from:
        #
        # Kalra, Paulin, Sedlazeck. "A systematic assessment of machine
        # learning for structural variant filtering." bioRxiv (2026): 2026-01.
        #
        # The script is run in parallel over VCF chunks, since it can be very
        # slow.
        #
        function FeatureExtraction() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local ALIGNMENTS_BAM=$3
            
            bcftools view --header-only ${INPUT_VCF} > ${SAMPLE_ID}_header.txt
            local N_RECORDS=$(bcftools view --no-header ${INPUT_VCF} | wc -l)
            local N_RECORDS_PER_THREAD=$(( (${N_RECORDS} + ${N_THREADS} - 1) / ${N_THREADS} ))
            bcftools view --no-header ${INPUT_VCF} | split -l ${N_RECORDS_PER_THREAD} - ${SAMPLE_ID}_chunk_
            for FILE in ${SAMPLE_ID}_chunk_*; do
                cat ${SAMPLE_ID}_header.txt ${FILE} > ${FILE}.vcf
                rm -f ${FILE}
                ${TIME_COMMAND} python ~{feature_extraction_py} ${FILE}.vcf ${ALIGNMENTS_BAM} ~{reference_fa} ${FILE}_features.csv 1>&2 &
            done
            wait
            cat ${SAMPLE_ID}_chunk_*_features.csv > ${SAMPLE_ID}_features.csv
            rm -f ${SAMPLE_ID}_chunk_*_features.csv ${SAMPLE_ID}_header.txt
            head -n 10 ${SAMPLE_ID}_features.csv 1>&2 || echo "0"
            tail -n +2 ${SAMPLE_ID}_features.csv | tr ',' '\t' | cut -f 1,2,6,8-19 | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=FEX_DEPTH_RATIO,Number=1,Type=Float,Description="depth_ratio from feature_extraction">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_DEPTH_MAD,Number=1,Type=Float,Description="depth_mad from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_AB,Number=1,Type=Float,Description="ab from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_CN_SLOP,Number=1,Type=Float,Description="cn_slop from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_MQ_DROP,Number=1,Type=Float,Description="mq_drop from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_CLIP_FRAC,Number=1,Type=Float,Description="clip_frac from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_SPLIT_READS,Number=1,Type=Integer,Description="split_reads from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_READ_LEN_MED,Number=1,Type=Float,Description="read_len_med from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_STRAND_BIAS,Number=1,Type=Float,Description="strand_bias from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_GC_FRAC,Number=1,Type=Float,Description="gc_frac from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_HOMOPOLYMER_MAX,Number=1,Type=Integer,Description="homopolymer_max from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=FEX_LCR_MASK,Number=1,Type=Integer,Description="lcr_mask from feature_extraction">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/FEX_DEPTH_RATIO,INFO/FEX_DEPTH_MAD,INFO/FEX_AB,INFO/FEX_CN_SLOP,INFO/FEX_MQ_DROP,INFO/FEX_CLIP_FRAC,INFO/FEX_SPLIT_READS,INFO/FEX_READ_LEN_MED,INFO/FEX_STRAND_BIAS,INFO/FEX_GC_FRAC,INFO/FEX_HOMOPOLYMER_MAX,INFO/FEX_LCR_MASK'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }
        
        
        
        
        # -------------------- Annotations from genotypers ---------------------
        
        # Uses all available cores
        #
        function Cutefc() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local ALIGNMENTS_BAM=$3
            local EXTRACT_BAM=$4

            local SLACK_BP=1000  # Arbitrary

            if [ ${EXTRACT_BAM} -eq 1 ]; then
                java -cp ~{docker_dir} UltralongGetIntervals ${INPUT_VCF} ~{reference_fai} ${SLACK_BP} > ${SAMPLE_ID}_regions.bed
                ${TIME_COMMAND} samtools view --threads ${N_THREADS} --bam --regions-file ${SAMPLE_ID}_regions.bed ${ALIGNMENTS_BAM} --output ${SAMPLE_ID}_extracted.bam
                ${TIME_COMMAND} samtools index --threads ${N_THREADS} ${SAMPLE_ID}_extracted.bam
                ALIGNMENTS_BAM=${SAMPLE_ID}_extracted.bam
            fi
            
            mkdir ./cutefc_dir/
            ${TIME_COMMAND} cuteFC --threads ${N_THREADS} --genotype --max_size -1 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -Ivcf ${INPUT_VCF} ${ALIGNMENTS_BAM} ~{reference_fa} ${SAMPLE_ID}_cutefc.vcf ./cutefc_dir
            rm -rf ./cutefc_dir
            bcftools query --format '%CHROM\t%POS\t%ID\t[%GT]\t[%GQ]\t[%DR]\t[%DV]\t[%PL]\t%INFO/CIPOS\t%INFO/CILEN\t%INFO/RE\t%INFO/STRAND\n' ${SAMPLE_ID}_cutefc.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                GT_COUNT=-1; \
                if ($4=="0/0" || $4=="0|0" || $4=="./."  || $4==".|." || $4=="./0" || $4==".|0" || $4=="0/." || $4=="0|." || $4=="0" || $4==".") GT_COUNT=0; \
                else if ($4=="0/1" || $4=="0|1" || $4=="1/0" || $4=="1|0" || $4=="./1" || $4==".|1" || $4=="1/." || $4=="1|." || $4=="1") GT_COUNT=1; \
                else if ($4=="1/1" || $4=="1|1") GT_COUNT=2; \
                \
                if ($5==".") GQ=-1; \
                else GQ=$5; \
                \
                if ($6==".") DR=-1; \
                else DR=$6; \
                \
                if ($7==".") DV=-1; \
                else DV=$7; \
                \
                PL_1=-1; PL_2=-1; PL_3=-1; \
                p=0; \
                for (i=1; i<=length($8); i++) { \
                    if (substr($8,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    PL_1=substr($8,1,p-1); \
                    q=0; \
                    for (i=p+1; i<=length($8); i++) { \
                        if (substr($8,i,1)==",") { q=i; break; } \
                    } \
                    if (q>0) { \
                        PL_2=substr($8,p+1,q-1-p); \
                        PL_3=substr($8,q+1); \
                    } \
                    else { PL_2=substr($8,p+1); } \
                } \
                else { PL_1=$8; }
                \
                CIPOS_1=-1; CIPOS_2=-1; \
                p=0; \
                for (i=1; i<=length($9); i++) { \
                    if (substr($9,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    CIPOS_1=substr($9,1,p-1); \
                    CIPOS_2=substr($9,p+1); \
                } \
                else { CIPOS_1=$9; } \
                \
                CILEN_1=-1; CILEN_2=-1; \
                p=0; \
                for (i=1; i<=length($10); i++) { \
                    if (substr($10,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    CILEN_1=substr($10,1,p-1); \
                    CILEN_2=substr($10,p+1); \
                } \
                else { CILEN_1=$10; } \
                \
                if ($11==".") RE=-1; \
                else RE=$11; \
                \
                STRAND=-1; \
                if ($12=="--") STRAND=0; \
                else if ($12=="-+") STRAND=1; \
                else if ($12=="+-") STRAND=2; \
                else if ($12=="++") STRAND=3; \
                \
                printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n",$1,$2,$3,GT_COUNT,GQ,DR,DV,PL_1,PL_2,PL_3,CIPOS_1,CIPOS_2,CILEN_1,CILEN_2,RE,STRAND); \
            }' | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_cutefc.vcf
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=CUTEFC_GT_COUNT,Number=1,Type=Integer,Description="Cutefc GT converted to an integer in {0,1,2}.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_GQ,Number=1,Type=Integer,Description="Genotype quality according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_DR,Number=1,Type=Integer,Description="High-quality reference reads according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_DV,Number=1,Type=Integer,Description="High-quality variant reads according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_PL_1,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_PL_2,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_PL_3,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CIPOS_1,Number=1,Type=Integer,Description="Confidence interval around POS for imprecise variants according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CIPOS_2,Number=1,Type=Integer,Description="Confidence interval around POS for imprecise variants according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CILEN_1,Number=1,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_CILEN_2,Number=1,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_RE,Number=1,Type=Integer,Description="Number of read support this record according to cutefc">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=CUTEFC_STRAND,Number=1,Type=Integer,Description="Cutefc strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--) converted to an integer in {0,1,2,3}.">' >> ${SAMPLE_ID}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/CUTEFC_GT_COUNT,INFO/CUTEFC_GQ,INFO/CUTEFC_DR,INFO/CUTEFC_DV,INFO/CUTEFC_PL_1,INFO/CUTEFC_PL_2,INFO/CUTEFC_PL_3,INFO/CUTEFC_CIPOS_1,INFO/CUTEFC_CIPOS_2,INFO/CUTEFC_CILEN_1,INFO/CUTEFC_CILEN_2,INFO/CUTEFC_RE,INFO/CUTEFC_STRAND'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt

            if [ ${EXTRACT_BAM} -eq 1 ]; then
                rm -f ${SAMPLE_ID}_extracted.bam*
            fi
        }
        
        


        # ------------------------ Repeat annotations --------------------------

        function VcfToBed_Start() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
      
            bcftools query --format '%CHROM\t%POS\t%INFO/SVLEN\t%ID\n' ${INPUT_VCF} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%d\t%s\n",$1,$2,$2,$4); }' | sort -k1,1 -k2,2n > ${SAMPLE_ID}_start.bed
        }


        function VcfToBed_StartEndInterval() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
      
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%INFO/SVLEN\t%ID\n' ${INPUT_VCF} > ${SAMPLE_ID}_matrix.tsv
            ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%d\t%s\n",$1,$2,$2,$4); }' ${SAMPLE_ID}_matrix.tsv | sort -k1,1 -k2,2n > ${SAMPLE_ID}_start.bed
            ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%d\t%s\n",$1, $2 + $3 - 1, $2 + $3 - 1, $4); }' ${SAMPLE_ID}_matrix.tsv | sort -k1,1 -k2,2n > ${SAMPLE_ID}_end.bed
            ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%d\t%s\n",$1,$2, $2 + $3 - 1, $4); }' ${SAMPLE_ID}_matrix.tsv | sort -k1,1 -k2,2n > ${SAMPLE_ID}_interval.bed
            rm -f ${SAMPLE_ID}_matrix.tsv
        }


        function AnnotateTrack_Point() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local POINT_BED=$3
            local POINT_ID=$4
            local TRACK_BED=$5
            local TRACK_ID=$6

            ${TIME_COMMAND} bedtools intersect -wa -u -a ${POINT_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t1\n",$1,$2,$4); }' > ${SAMPLE_ID}_${POINT_ID}_track.tsv
            ${TIME_COMMAND} bedtools intersect -wa -v -a ${POINT_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t0\n",$1,$2,$4); }' >> ${SAMPLE_ID}_${POINT_ID}_track.tsv
            sort -k 1,1 -k 2,2n ${SAMPLE_ID}_${POINT_ID}_track.tsv | bgzip > ${SAMPLE_ID}_${POINT_ID}_track.tsv.gz
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_${POINT_ID}_track.tsv.gz
            echo '##INFO=<ID='${POINT_ID}'_'${TRACK_ID}',Number=1,Type=Integer,Description="'${POINT_ID}' breakpoint is contained in a '${TRACK_ID}'">' > ${SAMPLE_ID}_header.txt
            COLUMNS='CHROM,POS,~ID,INFO/'${POINT_ID}'_'${TRACK_ID}
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_${POINT_ID}_track.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_${POINT_ID}_track* ${SAMPLE_ID}_header.txt
        }


        function AnnotateTrack_Interval() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INTERVAL_BED=$3
            local TRACK_BED=$4
            local TRACK_ID=$5
            local OVERLAP_FRACTION=$6

            ${TIME_COMMAND} bedtools intersect -wa -u -f ${OVERLAP_FRACTION} -a ${INTERVAL_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t1\n",$1,$2,$4); }' > ${SAMPLE_ID}_interval_track.tsv
            ${TIME_COMMAND} bedtools intersect -wa -v -f ${OVERLAP_FRACTION} -a ${INTERVAL_BED} -b ${TRACK_BED} | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t0\n",$1,$2,$4); }' >> ${SAMPLE_ID}_interval_track.tsv
            sort -k 1,1 -k 2,2n ${SAMPLE_ID}_interval_track.tsv | bgzip > ${SAMPLE_ID}_interval_track.tsv.gz
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_interval_track.tsv.gz
            echo '##INFO=<ID=INTERVAL_'${TRACK_ID}',Number=1,Type=Integer,Description="Interval overlaps the '${TRACK_ID}' track by at least '${OVERLAP_FRACTION}'">' > ${SAMPLE_ID}_header.txt
            COLUMNS='CHROM,POS,~ID,INFO/INTERVAL_'${TRACK_ID}
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_interval_track.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_interval_track* ${SAMPLE_ID}_header.txt
        }


        
        
        # ---------------- Handling alternative representations ----------------
        # Remarks: 
        #
        # 1. Some INS records may correspond to a DUP. Such records could be 
        # located at the beginning/end of their DUP, or in the middle of the 
        # DUP. INS records located at DUP breakpoints would see peculiar BAM 
        # features at their POS, while those located in the middle of their DUP 
        # would likely see no significant BAM features.
        #
        # 2. Assume that we keep such records in the INS VCF, without any 
        # special treatment.
        # 2.1 Assume that we then use a stringent match against both the INS and
        # the DUP truth VCF (with --dup-to-ins) to mark the true INS records. 
        # This would mark as true only INS records at the start (and perhaps at 
        # the end) of their DUP. The model would then learn to mark as true the 
        # INS records that have the BAM features of the true INS records, as 
        # well as those that have the BAM features of DUP breakpoints. This is 
        # fine, but INS records located in the middle of their DUP would not 
        # have such features and would likely be marked as false at test time:
        # this may result in losing some DUPs completely.
        # 2.2 Assume instead that we use a lenient match proportional to SVLEN
        # to mark the true INS records. This would mark as true also INS in the
        # middle of their DUP. This is fine, but since such INS have no BAM 
        # features, the model's precision may decrease.
        #
        # 3. We try to rewrite INS records as DUP using a simple heuristic that 
        # detects a longest interval with increased BAM coverage and position
        # and length compatible with the INS.
        # 3.1 Note that an INS could encode a complex DUP that corresponds to a
        # permutation of the source interval, and coverage inside the source 
        # interval may be uneven depending on how many times each source 
        # fragment occurs in the INS sequence. Our method puts in the DUP class
        # both simple and complex DUPs.
        # 3.2 We could have used `truvari anno remap`, rather than depth: this 
        # might have given more accurate breakpoints, and it might have allowed 
        # us to separate simple DUPs from complex DUPs. However, we want to 
        # perform the same classification at testing time, and we observed that
        # mapping some INS can take a large amount of RAM and thus make the task
        # expensive.
        #
        # 4. None of these problems affect the main (non-ultralong) VCF, since
        # features come from Kanpig which simply takes paths in the variation
        # graph, regardless of how they are represented in the VCF, and since 
        # the records that are marked as true by `truvari bench` are likely to 
        # capture the Kanpig features of true variants.
        #
        # 5. Alternative representations should have been handled in
        # `SV_Integration_Workpackage1.wdl` before truvari collapse. We do it
        # here just to reuse the intra-sample VCFs we already have from that
        # workflow.


        cat << 'END' > interval_2_breakpoints.sh
#!/bin/bash
set -euxo pipefail

CLASSPATH=$1
INPUT_BAM=$2
BIN_LENGTH=$3
BIN_COVERAGE_RATIO=$4
REGION=$5
ID=$6

samtools depth -aa -r ${REGION} ${INPUT_BAM} -o ${ID}_depth.tsv
java -cp ${CLASSPATH} UltralongDepthGetBreakpoints ${ID}_depth.tsv $(wc -l < ${ID}_depth.tsv) ${BIN_LENGTH} ${BIN_COVERAGE_RATIO} > ${ID}_breakpoints.tsv
rm -f ${ID}_depth.tsv
END
        chmod +x interval_2_breakpoints.sh


        function Ins2Dup() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local INPUT_BAM=$3
            local BIN_LENGTH=$4
            local BIN_COVERAGE_RATIO=$5

            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "SVTYPE=\"INS\"" --output-type v ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_ins.vcf
            local N_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_ins.vcf | wc -l)
            if [ ${N_INS} -eq 0 ]; then
                ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "SVTYPE!=\"INS\"" --output-type v ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_not_ins.vcf
                return
            fi
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "SVTYPE!=\"INS\"" --output-type z ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_not_ins.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_not_ins.vcf.gz
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongInsGetIntervals ${SAMPLE_ID}_ins.vcf ~{reference_fai} > ${SAMPLE_ID}_ins_intervals.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_ins_intervals.wsv --max-lines=1 --max-procs=${N_THREADS} ./interval_2_breakpoints.sh ~{docker_dir} ${INPUT_BAM} ${BIN_LENGTH} ${BIN_COVERAGE_RATIO}
            rm -f ${SAMPLE_ID}_ins_intervals.wsv
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongInsExtractDups ${SAMPLE_ID}_ins.vcf . ${SAMPLE_ID}_ins_ins.vcf ${SAMPLE_ID}_ins_dup.vcf
            rm -f *_breakpoints.tsv
            local N_INS_DUP=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_ins_dup.vcf | wc -l)
            if [ ${N_INS_DUP} -eq 0 ]; then
                rm -f ${SAMPLE_ID}_ins_dup.vcf ${SAMPLE_ID}_ins_ins.vcf
                gunzip ${SAMPLE_ID}_not_ins.vcf.gz
                return
            fi
            rm -f ${SAMPLE_ID}_ins.vcf ; mv ${SAMPLE_ID}_ins_ins.vcf ${SAMPLE_ID}_ins.vcf
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_ins_dup.vcf --output ${SAMPLE_ID}_ins_dup.vcf.gz
            rm -f ${SAMPLE_ID}_ins_dup.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins_dup.vcf.gz
            # Merging `ins_dup` and `not_ins`.
            # Remark: truvari collapse is run with the same parameters as in
            # `SV_Integration_Workpackage1.wdl`. 
            # Remark: truvari needs `bcftools merge` and it does not work with 
            # `bcftools concat`.
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type z ${SAMPLE_ID}_not_ins.vcf.gz ${SAMPLE_ID}_ins_dup.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_not_ins.vcf.gz* ${SAMPLE_ID}_ins_dup.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_not_ins.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_not_ins.vcf.gz
            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_not_ins.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0 --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_not_ins.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_not_ins.vcf
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        samtools --version 1>&2
        bcftools --version 1>&2
        truvari --help 1>&2
        cuteFC --version 1>&2
        truvari --help 1>&2
        df -h 1>&2
        
        while read -u 3 LINE; do
            # Skipping the sample if it has already been processed
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi

            # 1. Canonizing the VCF and handling alternative representations
            LocalizeSample ${SAMPLE_ID} ${LINE}
            df -h 1>&2
            CanonizeVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz
            Ins2Dup ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}.bam ~{ins2dup_bin_length} ~{ins2dup_bin_coverage_ratio}
            rm -f ${SAMPLE_ID}_canonized.vcf.gz

            # # 2. Adding custom annotations
            # N_NOT_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_not_ins.vcf | wc -l)
            # if [ ${N_NOT_INS} -gt 0 ]; then
            #     AnnotateCustom_NotIns ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf
            # fi
            # N_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_ins.vcf | wc -l)
            # if [ ${N_INS} -gt 0 ]; then
            #     AnnotateCustom_Ins ${SAMPLE_ID} ${SAMPLE_ID}_ins.vcf
            # fi

            # # 3. Adding repeat annotations
            # # 3.1 Not INS
            # if [ ${N_NOT_INS} -gt 0 ]; then
            #     VcfToBed_StartEndInterval ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf
            #     AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_start.bed "START" ~{tr_bed} "TR"
            #     rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
            #     AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_end.bed "END" ~{tr_bed} "TR"
            #     rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
            #     AnnotateTrack_Interval ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_interval.bed ~{tr_bed} "TR" ~{repeat_overlap_fraction}
            #     rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
            #     AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_start.bed "START" ~{segdup_bed} "SD"
            #     rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
            #     AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_end.bed "END" ~{segdup_bed} "SD"
            #     rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
            #     AnnotateTrack_Interval ${SAMPLE_ID} ${SAMPLE_ID}_not_ins.vcf ${SAMPLE_ID}_interval.bed ~{segdup_bed} "SD" ~{repeat_overlap_fraction}
            #     rm -f ${SAMPLE_ID}_not_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_not_ins.vcf
            # fi
            # # 3.2 INS
            # if [ ${N_INS} -gt 0 ]; then
            #     VcfToBed_Start ${SAMPLE_ID} ${SAMPLE_ID}_ins.vcf
            #     AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_ins.vcf ${SAMPLE_ID}_start.bed "START" ~{tr_bed} "TR"
            #     rm -f ${SAMPLE_ID}_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_ins.vcf
            #     AnnotateTrack_Point ${SAMPLE_ID} ${SAMPLE_ID}_ins.vcf ${SAMPLE_ID}_start.bed "START" ~{segdup_bed} "SD"
            #     rm -f ${SAMPLE_ID}_ins.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_ins.vcf
            #     rm -f ${SAMPLE_ID}_start.bed
            # fi

            # 4. Merging INS and non-INS VCFs
            ${TIME_COMMAND} bgzip --compress-level 1 ${SAMPLE_ID}_not_ins.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_not_ins.vcf.gz
            ${TIME_COMMAND} bgzip --compress-level 1 ${SAMPLE_ID}_ins.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins.vcf.gz
            ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type v ${SAMPLE_ID}_not_ins.vcf.gz ${SAMPLE_ID}_ins.vcf.gz --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_not_ins.vcf* ${SAMPLE_ID}_ins.vcf* ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf

            # # 5. Adding annotations from Kalra et al.
            # FeatureExtraction ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam
            # rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf

            # 6. Adding annotations from genotypers
            # Cutefc ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam 1
            # rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf






            cp ${SAMPLE_ID}_in.vcf original.vcf

            Cutefc ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam 0
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            bcftools query --format '%CHROM\t%POS\t%ID\t%INFO/CUTEFC_GT_COUNT\t%INFO/CUTEFC_GQ\t%INFO/CUTEFC_DR\t%INFO/CUTEFC_DV\t%INFO/CUTEFC_PL_1\t%INFO/CUTEFC_PL_2\t%INFO/CUTEFC_PL_3\t%INFO/CUTEFC_CIPOS_1\t%INFO/CUTEFC_CIPOS_2\t%INFO/CUTEFC_CILEN_1\t%INFO/CUTEFC_CILEN_2\t%INFO/CUTEFC_RE\t%INFO/CUTEFC_STRAND\n' ${SAMPLE_ID}_in.vcf > original.tsv
            cat original.tsv 1>&2

            Cutefc ${SAMPLE_ID} original.vcf ${SAMPLE_ID}.bam 1
            mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            bcftools query --format '%CHROM\t%POS\t%ID\t%INFO/CUTEFC_GT_COUNT\t%INFO/CUTEFC_GQ\t%INFO/CUTEFC_DR\t%INFO/CUTEFC_DV\t%INFO/CUTEFC_PL_1\t%INFO/CUTEFC_PL_2\t%INFO/CUTEFC_PL_3\t%INFO/CUTEFC_CIPOS_1\t%INFO/CUTEFC_CIPOS_2\t%INFO/CUTEFC_CILEN_1\t%INFO/CUTEFC_CILEN_2\t%INFO/CUTEFC_RE\t%INFO/CUTEFC_STRAND\n' ${SAMPLE_ID}_in.vcf > extracted.tsv
            cat extracted.tsv 1>&2

            diff --brief original.tsv extracted.tsv


            



            
            # Splitting by type and uploading
            bcftools filter --include "SVTYPE=\"DEL\"" --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_del.vcf.gz &
            bcftools filter --include "SVTYPE=\"DUP\"" --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_dup.vcf.gz &
            bcftools filter --include "SVTYPE=\"INV\"" --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_inv.vcf.gz &
            bcftools filter --include "SVTYPE=\"INS\"" --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_ins.vcf.gz &
            wait
            bcftools index -f -t ${SAMPLE_ID}_del.vcf.gz &
            bcftools index -f -t ${SAMPLE_ID}_dup.vcf.gz &
            bcftools index -f -t ${SAMPLE_ID}_inv.vcf.gz &
            bcftools index -f -t ${SAMPLE_ID}_ins.vcf.gz &
            wait
            gcloud storage mv ${SAMPLE_ID}_'del.vcf.gz*' ${SAMPLE_ID}_'inv.vcf.gz*' ${SAMPLE_ID}_'dup.vcf.gz*' ${SAMPLE_ID}_'ins.vcf.gz*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
            ls -laht 1>&2
        done 3< ~{chunk_csv}
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
