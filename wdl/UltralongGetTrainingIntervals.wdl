version 1.0


# Given an existing set of per-sample, filtered and annotated query VCFs, and
# a set of corresponding filtered and canonized dipcall VCFs and BEDs, the
# program selects, for each sample, the query records that match dipcall in the
# confident BED.
#
# Remark: sequence similarity is not used to decide a match.
#
workflow UltralongGetTrainingIntervals {
    input {
        File samples_tsv
        String suffix = "del"
        
        String remote_indir_annotated
        String remote_indir_truth
        String remote_outdir
        
        Int truvari_refdist = 10000000
        Float truvari_pctovl = 0.8
        Float truvari_pctsize = 0.8
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        samples_tsv: "Format: ID, DIPCALL_BED"
        remote_indir_annotated: "Without final slash. Contains per-sample annotated VCFs."
        remote_indir_truth: "Without final slash. Contains per-sample canonized and filtered dipcall VCFs."
        truvari_refdist: "For interval records we set it to a large value to disable the check."
        truvari_pctovl: "For interval records we enable this (it is disabled by default by truvari)."
        truvari_pctsize: "A non-stringent value"
    }
    
    call Impl {
        input:
            samples_tsv = samples_tsv,
            suffix = suffix,
            remote_indir_annotated = remote_indir_annotated,
            remote_indir_truth = remote_indir_truth,
            remote_outdir = remote_outdir,
            truvari_refdist = truvari_refdist,
            truvari_pctsize = truvari_pctsize,
            truvari_pctovl = truvari_pctovl,
            docker_image = docker_image
    }
    
    output {
    }
}


#
task Impl {
    input {
        File samples_tsv
        String suffix
        
        String remote_indir_annotated
        String remote_indir_truth
        String remote_outdir
        
        Int truvari_refdist
        Float truvari_pctsize
        Float truvari_pctovl
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 50
        Int preemptible_number = 0
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
        
        function GetTrainingIntervalsThread() {
            local CHUNK_CSV=$1
            
            while read -u 3 LINE; do
                SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
                DIPCALL_BED=$(echo ${LINE} | cut -d , -f 2)
                
                # Downloading the annotated VCF, skipping the sample if it was
                # not annotated.
                TEST=$( gsutil ls ~{remote_indir_annotated}/${SAMPLE_ID}_~{suffix}.vcf.gz || echo "1" )
                if [ ${TEST} == "1" ]; then
                    continue
                fi
                gcloud storage cp ~{remote_indir_annotated}/${SAMPLE_ID}_~{suffix}.vcf.gz ./${SAMPLE_ID}_query.vcf.gz
                gcloud storage cp ~{remote_indir_annotated}/${SAMPLE_ID}_~{suffix}.vcf.gz.tbi ./${SAMPLE_ID}_query.vcf.gz.tbi
                
                # Downloading and filtering the truth VCF
                gcloud storage cp ~{remote_indir_truth}/${SAMPLE_ID}_canonized.vcf.'gz*' .
                SVTYPE=""
                if [ ~{suffix} = "del" ]; then
                    SVTYPE="DEL"
                fi
                bcftools filter --include 'SVTYPE=="'${SVTYPE}'"' --output-type z ${SAMPLE_ID}_canonized.vcf.gz --output ${SAMPLE_ID}_truth.vcf.gz
                rm -f ${SAMPLE_ID}_canonized.vcf.gz*
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_truth.vcf.gz
                gcloud storage cp ${DIPCALL_BED} ./${SAMPLE_ID}_truth.bed
                
                # Computing matches
                if [ ~{truvari_refdist} -gt 1000 ]; then
                    # To avoid ERROR:root:--chunksize must be >= --refdist
                    CHUNKSIZE_FLAG="--chunksize ~{truvari_refdist}"
                else
                    CHUNKSIZE_FLAG=""
                fi
                ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_truth.vcf.gz -c ${SAMPLE_ID}_query.vcf.gz --includebed ${SAMPLE_ID}_truth.bed --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} ${CHUNKSIZE_FLAG} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
                mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_~{suffix}_training.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_~{suffix}_training.vcf.gz
                rm -rf ${SAMPLE_ID}_query.vcf.gz* ${SAMPLE_ID}_truth.vcf.gz* ${SAMPLE_ID}_truth.bed ${SAMPLE_ID}_truvari/
                
                # Uploading
                gcloud storage mv ${SAMPLE_ID}_~{suffix}_training.vcf.'gz*' ~{remote_outdir}
            done 3< ${CHUNK_CSV}
        }

        

        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        samtools --version 1>&2
        bcftools --version 1>&2
        truvari --help 1>&2
        df -h 1>&2
        
        cat ~{samples_tsv} | tr '\t' ',' > samples.csv
        N_ROWS=$(wc -l < samples.csv)
        if [ ${N_ROWS} -gt ${N_THREADS} ]; then
            N_ROWS_PER_THREAD=$(( ${N_ROWS} / ${N_THREADS} ))
            split -l ${N_ROWS_PER_THREAD} -d -a 4 samples.csv chunk_
        else
            mv samples.csv chunk_0
        fi
        for FILE in $(ls chunk_*); do
            GetTrainingIntervalsThread ${FILE} &
        done
        wait
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
