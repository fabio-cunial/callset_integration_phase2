version 1.0


# 
#
workflow UltralongAnnotate {
    input {
        File chunk_csv
        String remote_outdir
        
        File reference_fa
        File reference_fai
        
        File feature_extraction_py
        
        Int n_coverage_bins = 10
        Int breakpoint_window_bp = 500
        Int min_clip_length = 200
        Int adjacency_slack_bp = 300
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
        Int preemptible_number = 10
    }
    parameter_meta {
        chunk_csv: "Format: ID,bai,bam,csi,bcf"
        ultralong_training_resource_vcf_gz: "This should be a training resource specifically built for ultralong records. The standard training resource excludes such records."
    }
    
    call Impl {
        input:
            chunk_csv = chunk_csv,
            remote_outdir = remote_outdir,
            
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            
            feature_extraction_py = feature_extraction_py,
            
            n_coverage_bins = n_coverage_bins,
            breakpoint_window_bp = breakpoint_window_bp,
            min_clip_length = min_clip_length,
            adjacency_slack_bp = adjacency_slack_bp,
    
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    
    output {
    }
}


# Performance on a 4-core, 8GB VM:
#
# TOOL                                                CPU     RAM     TIME
# BAM download                                                          5m
#
# samtools bedcov                                     70%    500M       3m
# annotate_mapq_secondary.sh                         400%     15M      10s
# bcftools annotate                                  100%     15M      50s
# truvari bench                                      100%    200M      10s
# java UltralongIntervalGetBins                      200%     50M       1s
# java UltralongIntervalCreateBedcovAnnotations      200%     50M       1s
# annotate_clipped_alignments_1.sh                   400%    200M       1m
# annotate_clipped_alignments_2.sh                   400%     50M       1m
#
# sniffles 2.7.3                                      30%    500M      40m
# cutefc                                              30%    1.5G      50m
# lrcaller left                                       30%    200M      50m
# lrcaller right                                     300%    200M       4s
# feature_extraction.py                              100%    3.5G       2m
#
task Impl {
    input {
        File chunk_csv
        String remote_outdir
        
        File reference_fa
        File reference_fai
        
        File feature_extraction_py
        
        Int n_coverage_bins
        Int breakpoint_window_bp
        Int min_clip_length
        Int adjacency_slack_bp
        
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
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        # @param 
        # $2 A row of `chunk_csv`.
        #
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 2)
            ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 3)
            ULTRALONG_CSI=$(echo ${LINE} | cut -d , -f 4)
            ULTRALONG_BCF=$(echo ${LINE} | cut -d , -f 5)
            
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
            
            # Checking the integrity of the BAM
            ${TIME_COMMAND} samtools quickcheck -v ${SAMPLE_ID}.bam || { echo "ERROR: BAM corrupted"; exit 1; }
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
            
            DEFAULT_QUAL="60"   # Arbitrary
            
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
            
            # 3. Making sure IDs are distinct at the inter-sample level (they
            # are already distinct at the intra-sample level, thanks to the 
            # steps upstream).
            ${TIME_COMMAND} bcftools annotate --set-id ${SAMPLE_ID}'_%ID' --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            bcftools view --threads ${N_THREADS} --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_canonized.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_canonized.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf
        }
        
        
        
        
        # ------------------------- Custom annotations -------------------------
        
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

            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalGetBins ${INPUT_VCF} ~{reference_fai} ${N_BINS} ${BREAKPOINT_WINDOW_BP} > ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} samtools bedcov ${SAMPLE_ID}_bins.bed ${INPUT_BAM} > ${SAMPLE_ID}_counts.bed
            rm -f ${SAMPLE_ID}_bins.bed
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalCreateBedcovAnnotations ${SAMPLE_ID}_counts.bed $(( ${N_BINS} + 4 )) | sort -k 1,1 > ${SAMPLE_ID}_tags.tsv
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
            COLUMNS='CHROM,POS,~ID,INFO/BIN_POS'
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
            while read -u 4 ID; do
                MAPQ_LEFT=$(cat ${ID}_left_mapq.txt)
                MAPQ_RIGHT=$(cat ${ID}_right_mapq.txt)
                SECONDARY_LEFT=$(cat ${ID}_left_secondary.txt)
                SECONDARY_RIGHT=$(cat ${ID}_right_secondary.txt)
                echo -e "${ID}\t${MAPQ_LEFT}\t${MAPQ_RIGHT}\t${SECONDARY_LEFT}\t${SECONDARY_RIGHT}" >> ${SAMPLE_ID}_counts.tsv
                rm -f ${ID}_*_mapq.txt ${ID}_*_secondary.txt
            done 4< ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_variantID_sorted.txt
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
            while read -u 4 ID; do
                MAPQ_POINT=$(cat ${ID}_point_mapq.txt)
                SECONDARY_POINT=$(cat ${ID}_point_secondary.txt)
                echo -e "${ID}\t${MAPQ_POINT}\t${SECONDARY_POINT}" >> ${SAMPLE_ID}_counts.tsv
                rm -f ${ID}_*_mapq.txt ${ID}_*_secondary.txt
            done 4< ${SAMPLE_ID}_variantID_sorted.txt
            rm -f ${SAMPLE_ID}_variantID_sorted.txt
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF} | cut -f 1-3 | sort -k 3,3 > ${SAMPLE_ID}_chrom_pos_id.tsv
            ${TIME_COMMAND} join -t $'\t' -1 3 -2 1 ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t%f\t%d\n",$2,$3,$1,$4,$5); }' | sort -k 1,1 -k 2,2n | bgzip > ${SAMPLE_ID}_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_chrom_pos_id.tsv ${SAMPLE_ID}_counts.tsv
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=BIN_POINT_MAPQ,Number=1,Type=Float,Description="Breakpoint window: avg MAPQ.">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=BIN_POINT_SECONDARY,Number=1,Type=Integer,Description="Breakpoint window: number of secondary alignments.">' >> ${SAMPLE_ID}_header.txt
            COLUMNS='CHROM,POS,~ID,INFO/BIN_POINT_MAPQ,INFO/BIN_POINT_SECONDARY'
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

            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalGetBins ${INPUT_VCF} ~{reference_fai} 0 ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_1.sh ${INPUT_BAM} ~{docker_dir} ~{min_clip_length}
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
            COLUMNS='CHROM,POS,~ID,INFO/LL,INFO/LR,INFO/RL,INFO/RR,INFO/LL_RL_1,INFO/LL_RL_2,INFO/LL_RL_3,INFO/LL_RL_4,INFO/LL_RR_1,INFO/LL_RR_2,INFO/LL_RR_3,INFO/LL_RR_4,INFO/LR_RL_1,INFO/LR_RL_2,INFO/LR_RL_3,INFO/LR_RL_4,INFO/LR_RR_1,INFO/LR_RR_2,INFO/LR_RR_3,INFO/LR_RR_4'
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
    
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongPointGetBins ${INPUT_VCF} ~{reference_fai} ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_bins.wsv --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_1.sh ${INPUT_BAM} ~{docker_dir} ~{min_clip_length}
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
            COLUMNS='CHROM,POS,~ID,INFO/PL,INFO/PR,INFO/PL_PL_1,INFO/PL_PL_2,INFO/PL_PL_3,INFO/PL_PL_4,INFO/PL_PR_1,INFO/PL_PR_2,INFO/PL_PR_3,INFO/PL_PR_4,INFO/PR_PR_1,INFO/PR_PR_2,INFO/PR_PR_3,INFO/PR_PR_4'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }
        
        
        
        
        # ------------------- Annotations from Kalra et al. --------------------
        
        # Runs a slightly modified version of the code from:
        #
        # Kalra, Paulin, Sedlazeck. "A systematic assessment of machine
        # learning for structural variant filtering." bioRxiv (2026): 2026-01.
        #
        function FeatureExtraction() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local ALIGNMENTS_BAM=$3
            
            ${TIME_COMMAND} python ~{feature_extraction_py} ${INPUT_VCF} ${ALIGNMENTS_BAM} ~{reference_fa} 1>&2
            head -n 10 features.csv 1>&2 || echo "0"
            tail -n +2 features.csv | tr ',' '\t' | cut -f 1,2,6,8-19 | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
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
            COLUMNS='CHROM,POS,~ID,INFO/FEX_DEPTH_RATIO,INFO/FEX_DEPTH_MAD,INFO/FEX_AB,INFO/FEX_CN_SLOP,INFO/FEX_MQ_DROP,INFO/FEX_CLIP_FRAC,INFO/FEX_SPLIT_READS,INFO/FEX_READ_LEN_MED,INFO/FEX_STRAND_BIAS,INFO/FEX_GC_FRAC,INFO/FEX_HOMOPOLYMER_MAX,INFO/FEX_LCR_MASK'
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_annotated.vcf
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }
        
        
        
        
        # -------------------- Annotations from genotypers ---------------------
        
        # Single-core
        #
        function Sniffles() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local ALIGNMENTS_BAM=$3
            
            ${TIME_COMMAND} sniffles --threads 1 --input ${ALIGNMENTS_BAM} --reference ~{reference_fa} --genotype-vcf ${INPUT_VCF} --vcf ${SAMPLE_ID}_sniffles.vcf 1>&2
            bcftools query --format '%CHROM\t%POS\t%ID\t[%GT]\t[%GQ]\t[%DR]\t[%DV]\n' ${SAMPLE_ID}_sniffles.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
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
                printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\n",$1,$2,$3,GT_COUNT,GQ,DR,DV); \
            }' | bgzip -c > ${SAMPLE_ID}_annotations_sniffles.tsv.gz
            rm -f ${SAMPLE_ID}_sniffles.vcf
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations_sniffles.tsv.gz
            echo '##INFO=<ID=SNIFFLES_GT_COUNT,Number=1,Type=Integer,Description="Sniffles GT converted to an integer in {0,1,2}.">' > ${SAMPLE_ID}_header_sniffles.txt
            echo '##INFO=<ID=SNIFFLES_GQ,Number=1,Type=Integer,Description="Genotype quality according to sniffles">' >> ${SAMPLE_ID}_header_sniffles.txt
            echo '##INFO=<ID=SNIFFLES_DR,Number=1,Type=Integer,Description="Number of reference reads according to sniffles">' >> ${SAMPLE_ID}_header_sniffles.txt
            echo '##INFO=<ID=SNIFFLES_DV,Number=1,Type=Integer,Description="Number of variant reads according to sniffles">' >> ${SAMPLE_ID}_header_sniffles.txt
            echo 'CHROM,POS,~ID,INFO/SNIFFLES_GT_COUNT,INFO/SNIFFLES_GQ,INFO/SNIFFLES_DR,INFO/SNIFFLES_DV' > ${SAMPLE_ID}_columns_sniffles.txt
        }
        
        
        # Single-core
        #
        function Cutefc() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local ALIGNMENTS_BAM=$3
                
            mkdir ./cutefc_dir/
            ${TIME_COMMAND} cuteFC --threads 1 --genotype --max_size -1 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -Ivcf ${INPUT_VCF} ${ALIGNMENTS_BAM} ~{reference_fa} ${SAMPLE_ID}_cutefc.vcf ./cutefc_dir
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
            }' | bgzip -c > ${SAMPLE_ID}_annotations_cutefc.tsv.gz
            rm -f ${SAMPLE_ID}_cutefc.vcf
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations_cutefc.tsv.gz
            echo '##INFO=<ID=CUTEFC_GT_COUNT,Number=1,Type=Integer,Description="Cutefc GT converted to an integer in {0,1,2}.">' > ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_GQ,Number=1,Type=Integer,Description="Genotype quality according to cutefc">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_DR,Number=1,Type=Integer,Description="High-quality reference reads according to cutefc">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_DV,Number=1,Type=Integer,Description="High-quality variant reads according to cutefc">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_PL_1,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_PL_2,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_PL_3,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_CIPOS_1,Number=1,Type=Integer,Description="Confidence interval around POS for imprecise variants according to cutefc">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_CIPOS_2,Number=1,Type=Integer,Description="Confidence interval around POS for imprecise variants according to cutefc">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_CILEN_1,Number=1,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends according to cutefc">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_CILEN_2,Number=1,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends according to cutefc">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_RE,Number=1,Type=Integer,Description="Number of read support this record according to cutefc">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo '##INFO=<ID=CUTEFC_STRAND,Number=1,Type=Integer,Description="Cutefc strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--) converted to an integer in {0,1,2,3}.">' >> ${SAMPLE_ID}_header_cutefc.txt
            echo 'CHROM,POS,~ID,INFO/CUTEFC_GT_COUNT,INFO/CUTEFC_GQ,INFO/CUTEFC_DR,INFO/CUTEFC_DV,INFO/CUTEFC_PL_1,INFO/CUTEFC_PL_2,INFO/CUTEFC_PL_3,INFO/CUTEFC_CIPOS_1,INFO/CUTEFC_CIPOS_2,INFO/CUTEFC_CILEN_1,INFO/CUTEFC_CILEN_2,INFO/CUTEFC_RE,INFO/CUTEFC_STRAND' > ${SAMPLE_ID}_columns_cutefc.txt
        }
        
        
        # This is a separate script just to discard the STDERR generated by
        # lrcaller (which is too large and full of warnings about CIGAR parsing)
        # while keeping the STDERR of TIME_COMMAND.
        #
        cat << 'END' > lrcaller.sh
#!/bin/bash
set -euxo pipefail

SAMPLE_ID=$1
INPUT_VCF_GZ=$2
ALIGNMENTS_BAM=$3
BREAKPOINT=$4
REFERENCE_FA=$5

if [ ${BREAKPOINT} -eq 0 ]; then
    BREAKPOINT_FLAG=" "
else
    BREAKPOINT_FLAG="--right_breakpoint"
fi
lrcaller --number_of_threads 1 --dyn-w-size ${BREAKPOINT_FLAG} --fa ${REFERENCE_FA} ${ALIGNMENTS_BAM} ${INPUT_VCF_GZ} ${SAMPLE_ID}_lrcaller_${BREAKPOINT}_out.vcf 2> /dev/null
END
        chmod +x lrcaller.sh
        
        
        # Single-core.
        #
        # Remark: LRcaller can crash on some INS.
        #
        function Lrcaller() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local ALIGNMENTS_BAM=$3
            local BREAKPOINT=$4
            
            bcftools view --threads ${N_THREADS} --drop-genotypes --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_lrcaller_${BREAKPOINT}.vcf
            ${TIME_COMMAND} ./lrcaller.sh ${SAMPLE_ID} ${SAMPLE_ID}_lrcaller_${BREAKPOINT}.vcf ${ALIGNMENTS_BAM} ${BREAKPOINT} ~{reference_fa}
            rm -f ${SAMPLE_ID}_lrcaller_${BREAKPOINT}.vcf ; mv ${SAMPLE_ID}_lrcaller_${BREAKPOINT}_out.vcf ${SAMPLE_ID}_lrcaller_${BREAKPOINT}.vcf
            grep '^[^#]' ${SAMPLE_ID}_lrcaller_${BREAKPOINT}.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                printf("%s",$1); \
                for (i=2; i<=3; i++) printf("\t%s",$i); \
                for (i=10; i<=14; i++) { \
                    gsub(/[:,]/,"\t",$i); \
                    printf("\t%s",$i); \
                } \
                printf("\n"); \
            }' | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                GT_COUNT=-1; \
                if ($4=="0/0" || $4=="0|0" || $4=="./."  || $4==".|." || $4=="./0" || $4==".|0" || $4=="0/." || $4=="0|." || $4=="0" || $4==".") GT_COUNT=0; \
                else if ($4=="0/1" || $4=="0|1" || $4=="1/0" || $4=="1|0" || $4=="./1" || $4==".|1" || $4=="1/." || $4=="1|." || $4=="1") GT_COUNT=1; \
                else if ($4=="1/1" || $4=="1|1") GT_COUNT=2; \
                $4=GT_COUNT; \
                \
                GT_COUNT=-1; \
                if ($14=="0/0" || $14=="0|0" || $14=="./."  || $14==".|." || $14=="./0" || $14==".|0" || $14=="0/." || $14=="0|." || $14=="0" || $14==".") GT_COUNT=0; \
                else if ($14=="0/1" || $14=="0|1" || $14=="1/0" || $14=="1|0" || $14=="./1" || $14==".|1" || $14=="1/." || $14=="1|." || $14=="1") GT_COUNT=1; \
                else if ($14=="1/1" || $14=="1|1") GT_COUNT=2; \
                $14=GT_COUNT; \
                \
                GT_COUNT=-1; \
                if ($24=="0/0" || $24=="0|0" || $24=="./."  || $24==".|." || $24=="./0" || $24==".|0" || $24=="0/." || $24=="0|." || $24=="0" || $24==".") GT_COUNT=0; \
                else if ($24=="0/1" || $24=="0|1" || $24=="1/0" || $24=="1|0" || $24=="./1" || $24==".|1" || $24=="1/." || $24=="1|." || $24=="1") GT_COUNT=1; \
                else if ($24=="1/1" || $24=="1|1") GT_COUNT=2; \
                $24=GT_COUNT; \
                \
                GT_COUNT=-1; \
                if ($34=="0/0" || $34=="0|0" || $34=="./."  || $34==".|." || $34=="./0" || $34==".|0" || $34=="0/." || $34=="0|." || $34=="0" || $34==".") GT_COUNT=0; \
                else if ($34=="0/1" || $34=="0|1" || $34=="1/0" || $34=="1|0" || $34=="./1" || $34==".|1" || $34=="1/." || $34=="1|." || $34=="1") GT_COUNT=1; \
                else if ($34=="1/1" || $34=="1|1") GT_COUNT=2; \
                $34=GT_COUNT; \
                \
                GT_COUNT=-1; \
                if ($44=="0/0" || $44=="0|0" || $44=="./."  || $44==".|." || $44=="./0" || $44==".|0" || $44=="0/." || $44=="0|." || $44=="0" || $44==".") GT_COUNT=0; \
                else if ($44=="0/1" || $44=="0|1" || $44=="1/0" || $44=="1|0" || $44=="./1" || $44==".|1" || $44=="1/." || $44=="1|." || $44=="1") GT_COUNT=1; \
                else if ($44=="1/1" || $44=="1|1") GT_COUNT=2; \
                $44=GT_COUNT; \
                \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' | bgzip -c > ${SAMPLE_ID}_annotations_lrcaller_${BREAKPOINT}.tsv.gz
            rm -f ${SAMPLE_ID}_lrcaller_${BREAKPOINT}.vcf
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations_lrcaller_${BREAKPOINT}.tsv.gz
            if [ ${BREAKPOINT} -eq 0 ]; then
                SUFFIX="left"
            else
                SUFFIX="right"
            fi
            
            echo '##INFO=<ID=LRCALLER_GTCOUNT1_'${SUFFIX}',Number=1,Type=Integer,Description="Genotype count">' > ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_GTCOUNT2_'${SUFFIX}',Number=1,Type=Integer,Description="Genotype count">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_GTCOUNT3_'${SUFFIX}',Number=1,Type=Integer,Description="Genotype count">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_GTCOUNT4_'${SUFFIX}',Number=1,Type=Integer,Description="Genotype count">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_GTCOUNT5_'${SUFFIX}',Number=1,Type=Integer,Description="Genotype count">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_AD11_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_AD12_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_AD13_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_AD21_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_AD22_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_AD23_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_AD31_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_AD32_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_AD33_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_AD41_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_AD42_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_AD43_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_AD51_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_AD52_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_AD53_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_VA11_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_VA12_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_VA13_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_VA21_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_VA22_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_VA23_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_VA31_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_VA32_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_VA33_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_VA41_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_VA42_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_VA43_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_VA51_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_VA52_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_VA53_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_PL11_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_PL12_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_PL13_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_PL21_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_PL22_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_PL23_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_PL31_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_PL32_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_PL33_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_PL41_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_PL42_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_PL43_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            echo '##INFO=<ID=LRCALLER_PL51_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_PL52_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            echo '##INFO=<ID=LRCALLER_PL53_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ${SAMPLE_ID}_header_lrcaller_${BREAKPOINT}.txt
            
            if [ ${BREAKPOINT} -eq 0 ]; then
                echo 'CHROM,POS,~ID,INFO/LRCALLER_GTCOUNT1_'${SUFFIX}',INFO/LRCALLER_AD11_'${SUFFIX}',INFO/LRCALLER_AD12_'${SUFFIX}',INFO/LRCALLER_AD13_'${SUFFIX}',INFO/LRCALLER_VA11_'${SUFFIX}',INFO/LRCALLER_VA12_'${SUFFIX}',INFO/LRCALLER_VA13_'${SUFFIX}',INFO/LRCALLER_PL11_'${SUFFIX}',INFO/LRCALLER_PL12_'${SUFFIX}',INFO/LRCALLER_PL13_'${SUFFIX}',INFO/LRCALLER_GTCOUNT2_'${SUFFIX}',INFO/LRCALLER_AD21_'${SUFFIX}',INFO/LRCALLER_AD22_'${SUFFIX}',INFO/LRCALLER_AD23_'${SUFFIX}',INFO/LRCALLER_VA21_'${SUFFIX}',INFO/LRCALLER_VA22_'${SUFFIX}',INFO/LRCALLER_VA23_'${SUFFIX}',INFO/LRCALLER_PL21_'${SUFFIX}',INFO/LRCALLER_PL22_'${SUFFIX}',INFO/LRCALLER_PL23_'${SUFFIX}',INFO/LRCALLER_GTCOUNT3_'${SUFFIX}',INFO/LRCALLER_AD31_'${SUFFIX}',INFO/LRCALLER_AD32_'${SUFFIX}',INFO/LRCALLER_AD33_'${SUFFIX}',INFO/LRCALLER_VA31_'${SUFFIX}',INFO/LRCALLER_VA32_'${SUFFIX}',INFO/LRCALLER_VA33_'${SUFFIX}',INFO/LRCALLER_PL31_'${SUFFIX}',INFO/LRCALLER_PL32_'${SUFFIX}',INFO/LRCALLER_PL33_'${SUFFIX}',INFO/LRCALLER_GTCOUNT4_'${SUFFIX}',INFO/LRCALLER_AD41_'${SUFFIX}',INFO/LRCALLER_AD42_'${SUFFIX}',INFO/LRCALLER_AD43_'${SUFFIX}',INFO/LRCALLER_VA41_'${SUFFIX}',INFO/LRCALLER_VA42_'${SUFFIX}',INFO/LRCALLER_VA43_'${SUFFIX}',INFO/LRCALLER_PL41_'${SUFFIX}',INFO/LRCALLER_PL42_'${SUFFIX}',INFO/LRCALLER_PL43_'${SUFFIX}',INFO/LRCALLER_GTCOUNT5_'${SUFFIX}',INFO/LRCALLER_AD51_'${SUFFIX}',INFO/LRCALLER_AD52_'${SUFFIX}',INFO/LRCALLER_AD53_'${SUFFIX}',INFO/LRCALLER_VA51_'${SUFFIX}',INFO/LRCALLER_VA52_'${SUFFIX}',INFO/LRCALLER_VA53_'${SUFFIX}',INFO/LRCALLER_PL51_'${SUFFIX}',INFO/LRCALLER_PL52_'${SUFFIX}',INFO/LRCALLER_PL53_'${SUFFIX} > ${SAMPLE_ID}_columns_lrcaller_0.txt
            else
                echo 'CHROM,POS,~ID,INFO/LRCALLER_GTCOUNT1_'${SUFFIX}',INFO/LRCALLER_AD11_'${SUFFIX}',INFO/LRCALLER_AD12_'${SUFFIX}',INFO/LRCALLER_AD13_'${SUFFIX}',INFO/LRCALLER_VA11_'${SUFFIX}',INFO/LRCALLER_VA12_'${SUFFIX}',INFO/LRCALLER_VA13_'${SUFFIX}',INFO/LRCALLER_PL11_'${SUFFIX}',INFO/LRCALLER_PL12_'${SUFFIX}',INFO/LRCALLER_PL13_'${SUFFIX}',INFO/LRCALLER_GTCOUNT2_'${SUFFIX}',INFO/LRCALLER_AD21_'${SUFFIX}',INFO/LRCALLER_AD22_'${SUFFIX}',INFO/LRCALLER_AD23_'${SUFFIX}',INFO/LRCALLER_VA21_'${SUFFIX}',INFO/LRCALLER_VA22_'${SUFFIX}',INFO/LRCALLER_VA23_'${SUFFIX}',INFO/LRCALLER_PL21_'${SUFFIX}',INFO/LRCALLER_PL22_'${SUFFIX}',INFO/LRCALLER_PL23_'${SUFFIX}',INFO/LRCALLER_GTCOUNT3_'${SUFFIX}',INFO/LRCALLER_AD31_'${SUFFIX}',INFO/LRCALLER_AD32_'${SUFFIX}',INFO/LRCALLER_AD33_'${SUFFIX}',INFO/LRCALLER_VA31_'${SUFFIX}',INFO/LRCALLER_VA32_'${SUFFIX}',INFO/LRCALLER_VA33_'${SUFFIX}',INFO/LRCALLER_PL31_'${SUFFIX}',INFO/LRCALLER_PL32_'${SUFFIX}',INFO/LRCALLER_PL33_'${SUFFIX}',INFO/LRCALLER_GTCOUNT4_'${SUFFIX}',INFO/LRCALLER_AD41_'${SUFFIX}',INFO/LRCALLER_AD42_'${SUFFIX}',INFO/LRCALLER_AD43_'${SUFFIX}',INFO/LRCALLER_VA41_'${SUFFIX}',INFO/LRCALLER_VA42_'${SUFFIX}',INFO/LRCALLER_VA43_'${SUFFIX}',INFO/LRCALLER_PL41_'${SUFFIX}',INFO/LRCALLER_PL42_'${SUFFIX}',INFO/LRCALLER_PL43_'${SUFFIX}',INFO/LRCALLER_GTCOUNT5_'${SUFFIX}',INFO/LRCALLER_AD51_'${SUFFIX}',INFO/LRCALLER_AD52_'${SUFFIX}',INFO/LRCALLER_AD53_'${SUFFIX}',INFO/LRCALLER_VA51_'${SUFFIX}',INFO/LRCALLER_VA52_'${SUFFIX}',INFO/LRCALLER_VA53_'${SUFFIX}',INFO/LRCALLER_PL51_'${SUFFIX}',INFO/LRCALLER_PL52_'${SUFFIX}',INFO/LRCALLER_PL53_'${SUFFIX} > ${SAMPLE_ID}_columns_lrcaller_1.txt
            fi
        }
        
        
        # Sequentially adds the annotations from every genotyper to a VCF.
        #
        function TransferGenotypersAnnotations() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local LRCALLER_USE_RIGHT=$3
            
            COLUMNS=$(cat ${SAMPLE_ID}_columns_sniffles.txt | head -n 1)
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations_sniffles.tsv.gz --header-lines ${SAMPLE_ID}_header_sniffles.txt --columns ${COLUMNS} --output-type v ${INPUT_VCF} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            COLUMNS=$(cat ${SAMPLE_ID}_columns_cutefc.txt | head -n 1)
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations_cutefc.tsv.gz --header-lines ${SAMPLE_ID}_header_cutefc.txt --columns ${COLUMNS} --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            COLUMNS=$(cat ${SAMPLE_ID}_columns_lrcaller_0.txt | head -n 1)
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations_lrcaller_0.tsv.gz --header-lines ${SAMPLE_ID}_header_lrcaller_0.txt --columns ${COLUMNS} --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            if [ ${LRCALLER_USE_RIGHT} -eq 1 ]; then
                COLUMNS=$(cat ${SAMPLE_ID}_columns_lrcaller_1.txt | head -n 1)
                ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_annotations_lrcaller_1.tsv.gz --header-lines ${SAMPLE_ID}_header_lrcaller_1.txt --columns ${COLUMNS} --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            fi

            mv ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}_annotated.vcf
            
            rm -f ${SAMPLE_ID}_annotations_sniffles.tsv.gz* ${SAMPLE_ID}_header_sniffles.txt ${SAMPLE_ID}_columns_sniffles.txt
            rm -f ${SAMPLE_ID}_annotations_cutefc.tsv.gz* ${SAMPLE_ID}_header_cutefc.txt ${SAMPLE_ID}_columns_cutefc.txt
            rm -f ${SAMPLE_ID}_annotations_lrcaller_*.tsv.gz* ${SAMPLE_ID}_header_lrcaller_*.txt ${SAMPLE_ID}_columns_lrcaller_*.txt
        }
        
        
        
        
        # -------------------- Pipelines for each SV type ----------------------
        
        function AnnotateDels() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'SVTYPE="DEL"' --output-type v ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_out.vcf
            mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 1. Adding custom features (fast).
            AnnotateCoverageBins_Interval ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{n_coverage_bins} ~{breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            AnnotateMapqSecondary_Interval ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            AnnotateClippedAlignments_Interval ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{breakpoint_window_bp} ~{adjacency_slack_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            
            # 2. Adding features from Kalra et al. (fast).
            FeatureExtraction ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            
            # 3. Adding features from genotypers (slow).
            Sniffles ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam &
            Cutefc ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam &
            Lrcaller ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam 0 &
            Lrcaller ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam 1 &
            wait
            TransferGenotypersAnnotations ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf 1
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf            
            
            # Uploading
            bcftools view --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_del.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_del.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf
            gcloud storage mv ${SAMPLE_ID}_del.vcf.'gz*' ~{remote_outdir}/
        }
        
        
        function AnnotateIns() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'SVTYPE="INS"' --output-type v ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_out.vcf
            mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 1. Adding custom features (fast).
            AnnotateCoverageBins_Point ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            AnnotateMapqSecondary_Point ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{breakpoint_window_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            AnnotateClippedAlignments_Point ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam ~{breakpoint_window_bp} ~{adjacency_slack_bp}
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            
            # 2. Adding features from Kalra et al. (fast).
            FeatureExtraction ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            
            # 3. Adding features from genotypers (slow).
            Sniffles ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam &
            Cutefc ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam &
            Lrcaller ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf ${SAMPLE_ID}.bam 0 &
            wait
            TransferGenotypersAnnotations ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf 0
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_annotated.vcf ${SAMPLE_ID}_in.vcf
            
            # Uploading
            bcftools view --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_ins.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf
            gcloud storage mv ${SAMPLE_ID}_ins.vcf.'gz*' ~{remote_outdir}/
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        samtools --version 1>&2
        bcftools --version 1>&2
        truvari --help 1>&2
        sniffles --version 1>&2
        cuteFC --version 1>&2
        lrcaller --version 1>&2
        df -h 1>&2
        
        # Processing samples
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            # Skipping the sample if it has already been processed
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            LocalizeSample ${SAMPLE_ID} ${LINE}
            df -h 1>&2
            CanonizeVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz
            #AnnotateDels ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz
            AnnotateIns ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz
            
            
# -----------> Each genotyper should be called just once on all records, since its annotations do not depend on SVTYPE. So we should first annotate each SVTYPE with our custom code in isolation per type, then concatenate the annotated VCFs, and then run the genotypers on all calls.
            
            
            # Other SV types are omitted for now
            
            # Next iteration
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
