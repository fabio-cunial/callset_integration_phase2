version 1.0

#
workflow InvestigateComplexDuplicationsPrime {
    input {
        File chunk_tsv
        String remote_outdir

        File reference_fa
        File reference_fai

        Float truvari_anno_remap_cov_threshold = 0.8

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong_remap:latest"
    }
    parameter_meta {
        chunk_tsv: "Format: `sample_id, ultralong_bcf, hap1_fa, hap2_fa, truth_vcf_gz, truth_tbi`, where the truth VCFs are assumed to contain all SV types (to allow comparing query INS to truth INS and DUP)."
    }

    call Impl {
        input:
            chunk_tsv = chunk_tsv,
            remote_outdir = remote_outdir,

            reference_fa = reference_fa,
            reference_fai = reference_fai,

            truvari_anno_remap_cov_threshold = truvari_anno_remap_cov_threshold,

            docker_image = docker_image
    }
}


# TOOL                                                CPU     RAM     TIME
# 
# 
# 
task Impl {
    input {
        File chunk_tsv
        String remote_outdir

        File reference_fa
        File reference_fai

        Float truvari_anno_remap_cov_threshold

        String docker_image
        Int n_cpu = 3
        Int mem_gb = 25
        Int disk_size_gb = 50
        Int preemptible_number = 0
    }

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        EFFECTIVE_RAM_MB=$(( ~{mem_gb}*1024 - 500 ))
        



        # --------------------- Steps of the pipeline --------------------------

        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            local ULTRALONG_BCF=$(echo ${LINE} | cut -d , -f 2)
            local HAP1_FA=$(echo ${LINE} | cut -d , -f 3)
            local HAP2_FA=$(echo ${LINE} | cut -d , -f 4)
            local TRUTH_VCF_GZ=$(echo ${LINE} | cut -d , -f 5)
            local TRUTH_TBI=$(echo ${LINE} | cut -d , -f 6)

            gcloud storage cp ${ULTRALONG_BCF} ${SAMPLE_ID}.bcf
            gcloud storage cp ${TRUTH_VCF_GZ} ${SAMPLE_ID}_truth.vcf.gz
            gcloud storage cp ${TRUTH_TBI} ${SAMPLE_ID}_truth.vcf.gz.tbi
        }


        function DelocalizeSample() {
            local SAMPLE_ID=$1

            rm -f ${SAMPLE_ID}.bcf
            rm -f ${SAMPLE_ID}_truth.vcf.gz
            rm -f ${SAMPLE_ID}_truth.vcf.gz.tbi
        }




        # -------------------------- Main program ------------------------------
    
        INFINITY="3000000000"
        MIN_SV_LENGTH="5000"  # Arbitrary
        TRUVARI_BENCH_SIMILARITY_THRESHOLD="0.7"  # Arbitrary lower bound

        # Preparing the reference for `truvari anno remap`.
        #
        # Remark: the reference should be minimap-indexed only once and reused
        # by every sample.
        mkdir ./ref_files
        mv ~{reference_fa} ./ref_files/
        mv ~{reference_fai} ./ref_files/
        REF_FA_BASENAME=$(basename ~{reference_fa})

        # Processing each sample in the chunk file
        # Remark: we assume that INSDUPs have been converted back to INS.
        cat ~{chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE || [ -n "${LINE}" ]; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            LocalizeSample ${SAMPLE_ID} ${LINE}
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'INFO/SVTYPE=="INS"' --output-type z ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}_ins.vcf.gz
            rm -f ${SAMPLE_ID}.bcf* ; tabix -f ${SAMPLE_ID}_ins.vcf.gz
            ${TIME_COMMAND} truvari anno remap --threads ${N_THREADS} --aligner minimap2 --min-length 1 --max-length 3000000000 --cov-threshold ~{truvari_anno_remap_cov_threshold} -r ./ref_files/${REF_FA_BASENAME} ${SAMPLE_ID}_ins.vcf.gz -o ${SAMPLE_ID}_ins_annotated.vcf.gz
            rm -f ${SAMPLE_ID}_ins.vcf.gz*
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'INFO/remap_classification=="tandem" || INFO/remap_classification=="tandem_complex" || INFO/remap_classification=="tandem_inverted"' --output-type z ${SAMPLE_ID}_ins_annotated.vcf.gz --output ${SAMPLE_ID}_ins_remapped.vcf.gz
            rm -f ${SAMPLE_ID}_ins_annotated.vcf.gz* ; tabix -f ${SAMPLE_ID}_ins_remapped.vcf.gz
            N_REMAPPED_INS=$( bcftools index --nrecords ${SAMPLE_ID}_ins_remapped.vcf.gz )
            if [ ${N_REMAPPED_INS} -eq 0 ]; then
                rm -f ${SAMPLE_ID}_ins_remapped.vcf.gz*
                DelocalizeSample ${SAMPLE_ID}
                continue
            fi
            rm -rf ./${SAMPLE_ID}_truvari/
            # --chunksize must be >= --refdist
            ${TIME_COMMAND} truvari bench --reference ./ref_files/${REF_FA_BASENAME} --refdist ${INFINITY} --chunksize ${INFINITY} --sizemin ${MIN_SV_LENGTH} --sizefilt ${MIN_SV_LENGTH} --sizemax ${INFINITY} --dup-to-ins --pctovl 0 --pctseq ${TRUVARI_BENCH_SIMILARITY_THRESHOLD} --max-resolve ${INFINITY} --pctsize ${TRUVARI_BENCH_SIMILARITY_THRESHOLD} -b ${SAMPLE_ID}_truth.vcf.gz -c ${SAMPLE_ID}_ins_remapped.vcf.gz -o ./${SAMPLE_ID}_truvari/
            rm -f ${SAMPLE_ID}_ins_remapped.vcf.gz*
            ${TIME_COMMAND} bcftools query --format '%ID\t%INFO/SVLEN\t%INFO/TruScore\t%INFO/PctSeqSimilarity\t%INFO/PctSizeSimilarity\t%INFO/PctRecOverlap\t%INFO/StartDistance\t%INFO/EndDistance\t%INFO/SizeDiff\n' ${SAMPLE_ID}_truvari/tp-comp.vcf.gz > ${SAMPLE_ID}_ins_remapped_tp.tsv
            rm -rf ${SAMPLE_ID}_truvari/
            N_FOUND_INS=$( wc -l < ${SAMPLE_ID}_ins_remapped_tp.tsv )
            PERCENT=$( bc <<< "scale=2; (${N_FOUND_INS}*100)/${N_REMAPPED_INS}" )
            echo "${SAMPLE_ID} has ${N_FOUND_INS} \"true\" remapped sequences (ignoring distance) out of ${N_REMAPPED_INS} (${PERCENT}%)" 1>&2
            gcloud storage mv ${SAMPLE_ID}_ins_remapped_tp.tsv ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
        done 3< chunk.csv
    >>>
    
    output {
    }

    runtime {
        cpu: n_cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " +  disk_size_gb + " HDD"
        preemptible: preemptible_number
        docker: docker_image
    }
}
