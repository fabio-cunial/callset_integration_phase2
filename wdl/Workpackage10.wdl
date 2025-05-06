version 1.0


# 
#
workflow Workpackage10 {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 50
    }
    parameter_meta {
        sv_integration_chunk_tsv: "A subset of the rows of table `sv_integration_hg38`, without the header."
    }
    
    call Workpackage10Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}


# Performance on 10'070 samples, 15x, GRCh38:
#
# CPU     RAM     TIME
# 
#
task Workpackage10Impl {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir
        
        Int n_cpu
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int compression_level = 1
    
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
            local LINE=$2
            
            while : ; do
                TEST=$(gsutil -m cp ~{remote_indir}/${SAMPLE_ID}_gts.txt . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${SAMPLE_ID}_gts.txt>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        }
        
        
        # Deletes all the files related to a sample
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}_*
        }
        
        
        # Sorts by the first field, which is the row number in the inter-sample
        # VCF.
        #
        function Sort() {
            local SAMPLE_ID=$1
            
            echo "${SAMPLE_ID}" > ${SAMPLE_ID}_sorted.txt
            tail -n +2 ${SAMPLE_ID}_gts.txt > ${SAMPLE_ID}_tmp1.txt
            ${TIME_COMMAND} sort --numeric-sort -t '\t' -k 1,1 ${SAMPLE_ID}_tmp1.txt > ${SAMPLE_ID}_tmp2.txt
            rm -f ${SAMPLE_ID}_tmp1.txt
            cut -f 2 ${SAMPLE_ID}_tmp2.txt >> ${SAMPLE_ID}_sorted.txt
            rm -f ${SAMPLE_ID}_tmp2.txt
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            LocalizeSample ${SAMPLE_ID} ${LINE}
            Sort ${SAMPLE_ID}
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${SAMPLE_ID}_sorted.txt ~{remote_outdir}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading the sorted GT file. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
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
