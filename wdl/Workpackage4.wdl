version 1.0


# Splits every xgboost-scored single-sample VCF into ~100 pieces in order to run
# bcftools merge in parallel.
#
workflow Workpackage4 {
    input {
        File sv_integration_chunk_tsv
        File split_for_bcftools_merge_csv
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 20
    }
    parameter_meta {
        split_for_bcftools_merge_csv: "A partition that covers all chromosomes. Every line is a 0-based, half-open, consecutive chunk of a chromosome. Lines are assumed to be sorted."
    }
    
    call Workpackage4Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}


#
task Workpackage4Impl {
    input {
        File sv_integration_chunk_tsv
        File split_for_bcftools_merge_csv
        String remote_indir
        String remote_outdir
        
        Int n_cpu
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/root"
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local REMOTE_DIR=$2
            
            while : ; do
                TEST=$(gcloud storage cp ${REMOTE_DIR}/${SAMPLE_ID}_scored.vcf.'gz*' . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${REMOTE_DIR}/${SAMPLE_ID}_scored.vcf.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        }
        
        
        # Deletes all files and directories related to the sample
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -rf ./${SAMPLE_ID}_*
        }
        
        
        function Filter() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            i="0"
            while read INTERVAL; do
                echo ${INTERVAL} | tr ',' '\t' > ${SAMPLE_ID}.bed
                ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --regions-file ${SAMPLE_ID}.bed --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > ${SAMPLE_ID}_chunk_${i}.vcf.gz
                tabix -f ${SAMPLE_ID}_chunk_${i}.vcf.gz
                i=$(( ${i} + 1 ))
            done < ~{split_for_bcftools_merge_csv}
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${SAMPLE_ID}_chunk_'*'.vcf.'gz*' ~{remote_outdir}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading chunks for sample ${SAMPLE_ID}. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------

        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            LocalizeSample ${SAMPLE_ID} ~{remote_indir}
            Filter ${SAMPLE_ID} ${SAMPLE_ID}_scored.vcf.gz
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done < chunk.csv
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
