version 1.0


# 
#
workflow UltralongAnnotate {
    input {
        File chunk_csv
        String remote_outdir
        
        File reference_fai
        File ultralong_training_resource_bed
        File ultralong_training_resource_vcf_gz
        File ultralong_training_resource_tbi
        
        Int n_coverage_bins = 10
        Int breakpoint_window_bp = 500
        Int min_clip_length = 200
        Int adjacency_slack_bp = 300
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong"
        Int preemptible_number = 4
    }
    parameter_meta {
        chunk_csv: "Format: ID,bai,bam,csi,bcf"
        ultralong_training_resource_vcf_gz: "This should be a training resource specifically built for ultralong records. The standard training resource excludes such records."
    }
    
    call Impl {
        input:
            chunk_csv = chunk_csv,
            remote_outdir = remote_outdir,
            
            reference_fai = reference_fai,
            ultralong_training_resource_bed = ultralong_training_resource_bed,
            ultralong_training_resource_vcf_gz = ultralong_training_resource_vcf_gz,
            ultralong_training_resource_tbi = ultralong_training_resource_tbi,
    
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


# TOOL                                                CPU     RAM     TIME
# samtools bedcov           
# samtools view
# bcftools annotate
# truvari bench    
# java UltralongIntervalGetBins
# java UltralongIntervalCreateBedcovAnnotations
# java UltralongIntervalGetClips
# java UltralongIntervalIntersectClips
#
task Impl {
    input {
        File chunk_csv
        String remote_outdir
        
        File reference_fai
        File ultralong_training_resource_bed
        File ultralong_training_resource_vcf_gz
        File ultralong_training_resource_tbi
        
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
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        
        
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
            
            # Converting to .vcf.gz for the genotypers
            bcftools view --threads ${N_THREADS} --output-type z ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}.vcf.gz
            rm -f ${SAMPLE_ID}.bcf*
            
            # Checking the integrity of the BAM
            ${TIME_COMMAND} samtools quickcheck -v ${SAMPLE_ID}.bam || { echo "ERROR: BAM corrupted"; exit 1; }
        }
        
        
        # Deletes all and only the files downloaded by `LocalizeSample()`.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}*.bam* ${SAMPLE_ID}*.bcf* ${SAMPLE_ID}*.vcf.gz*
        }
        
        
        # Ensures that the VCF is correctly formatted
        #
        function CanonizeVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            # 1. Cleaning REF, ALT, QUAL, FILTER.
            # - REF and ALT must be uppercase for XGBoost scoring downstream to
            #   work.
            # - We force every record to PASS, to rule out any filter-dependent
            #   effect in downstream tools.
            DEFAULT_QUAL="60"   # Arbitrary
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanRefAltQual ${INPUT_VCF_GZ} ${DEFAULT_QUAL} > ${SAMPLE_ID}_out.vcf
            rm -f ${INPUT_VCF_GZ}* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 2.1 Removing END, since its values may be inconsistent and make
            # GATK crash downstream.
            # 2.2 Making sure IDs are distinct at the inter-sample level (they
            # are already distinct at the intra-sample level thanks to the 
            # steps upstream).
            ${TIME_COMMAND} bcftools annotate --remove INFO/END --set-id ${SAMPLE_ID}'_%ID' --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            bcftools view --threads ${N_THREADS} --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_canonized.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_canonized.vcf.gz
        }
        
        
        # Given an input VCF containing only interval records, the procedure
        # annotates each record with coverage measures extracted from a BAM.
        #
        # Remark: `samtools bedcov` skips reads with any of the following flags
        # set: UNMAP, SECONDARY, QCFAIL, DUP
        #
        function AnnotateCoverageBins() {
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
        function AnnotateMapqSecondary() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local INPUT_BAM=$3
            local BREAKPOINT_WINDOW_BP=$4
    
            ${TIME_COMMAND} java -cp ~{docker_dir} UltralongIntervalGetBins ${INPUT_VCF} ~{reference_fai} 0 ${BREAKPOINT_WINDOW_BP} | tr '\t' ' ' > ${SAMPLE_ID}_bins.wsv
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
                rm -f ${ID}_*_mapq.txt ${ID}_*_secondary.txt
            done < ${SAMPLE_ID}_variantID_sorted.txt
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


        # Given an input VCF containing only interval calls, the procedure
        # annotates it with clipped alignment measures extracted from a BAM.
        #
        function AnnotateClippedAlignments() {
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
            ${TIME_COMMAND} xargs --arg-file=${SAMPLE_ID}_variantID.txt --max-lines=1 --max-procs=${N_THREADS} ./annotate_clipped_alignments_2.sh ~{docker_dir} ${ADJACENCY_SLACK_BP}
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
        
        
        # Removes SVLEN from symbolic ALTs, in order not to interfere with
        # `truvari bench`.
        #
        function ResetAlts() {
            local INPUT_VCF=$1
            local OUTPUT_VCF=$2
            
            date 1>&2
            ( bcftools view --header-only ${INPUT_VCF} ; bcftools view --no-header ${INPUT_VCF} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                if (substr($0,1,1)!="#" && substr($5,1,1)=="<") $5 = substr($5,1,4) ">"; \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' ) > ${OUTPUT_VCF}
            date 1>&2
        }
        
        
        # Given an interval-only input VCF, the procedure compresses it and
        # computes its records with a stringent match to the training resource.
        #
        # Remark: sequence similarity is not used to decide a match.
        #
        function GetTrainingRecords() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local SUFFIX=$3
            
            bcftools view --output-type z ${INPUT_VCF} --output ${SAMPLE_ID}_${SUFFIX}.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${SUFFIX}.vcf.gz
            rm -f ${INPUT_VCF}
            
            ${TIME_COMMAND} truvari bench -b ~{ultralong_training_resource_vcf_gz} -c ${SAMPLE_ID}_${SUFFIX}.vcf.gz --includebed ~{ultralong_training_resource_bed} --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0 --pick single -o ./truvari_${SAMPLE_ID}/
            
            mv truvari_${SAMPLE_ID}/tp-comp.vcf.gz ${SAMPLE_ID}_${SUFFIX}_training.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${SUFFIX}_training.vcf.gz
            rm -rf truvari_${SAMPLE_ID}/
        }
        
        
        # DEL-focused workflow.
        #
        function ProcessDels() {
            local INPUT_VCF_GZ=$1
            local SAMPLE_ID=$2
            
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'SVTYPE="DEL"' --output-type v ${INPUT_VCF_GZ} --output out.vcf
            ResetAlts out.vcf in.vcf
            
            # Annotating
            AnnotateCoverageBins ${SAMPLE_ID} in.vcf ${SAMPLE_ID}.bam ~{n_coverage_bins} ~{breakpoint_window_bp}
            rm -f in.vcf ; mv ${SAMPLE_ID}_annotated.vcf in.vcf
            AnnotateMapqSecondary ${SAMPLE_ID} in.vcf ${SAMPLE_ID}.bam ~{breakpoint_window_bp}
            rm -f in.vcf ; mv ${SAMPLE_ID}_annotated.vcf in.vcf
            AnnotateClippedAlignments ${SAMPLE_ID} in.vcf ${SAMPLE_ID}.bam ~{breakpoint_window_bp} ~{adjacency_slack_bp}
            rm -f in.vcf ; mv ${SAMPLE_ID}_annotated.vcf in.vcf
            
            # Computing training records
            GetTrainingRecords ${SAMPLE_ID} in.vcf del
            
            # Uploading
            gcloud storage cp ${SAMPLE_ID}_del.vcf.'gz*' ${SAMPLE_ID}_del_training.vcf.'gz*' ~{remote_outdir}/
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        samtools --version 1>&2
        bcftools --version 1>&2
        truvari --help 1>&2
        df -h 1>&2
        
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            # Skipping the sample if it has already been processed
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            LocalizeSample ${SAMPLE_ID} ${LINE}
            CanonizeVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz
            ProcessDels ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}
            # Other SV types are omitted for now
            
            # Next iteration
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done < ~{chunk_csv}
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
