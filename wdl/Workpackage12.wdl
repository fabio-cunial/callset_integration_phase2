version 1.0


# Transfers some format fields from a chunk of the truvari-collapsed VCF to the
# corresponding chunk of the kanpig-regenotyped VCF.
#
workflow Workpackage12 {
    input {
        String chunk_id
        String remote_indir
        String remote_outdir
    }
    parameter_meta {
        remote_indir: "Containing both input chunks."
    }
    
    call Workpackage11Impl {
        input:
            chunk_id = chunk_id,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


#
task Workpackage11Impl {
    input {
        String chunk_id
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 1
        Int ram_size_gb = 4
        Int disk_size_gb = 50
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        # Localizing chunks
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{remote_indir}/chunk_old_~{chunk_id} . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading chunk. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{remote_indir}/chunk_new_~{chunk_id} . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading chunk. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        # Copying format fields from old to new
        ${TIME_COMMAND} java -cp ~{docker_dir} CopyFormatFast chunk_old_~{chunk_id} chunk_new_~{chunk_id} chunk_out_~{chunk_id} chunk_hwe_~{chunk_id}
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp chunk_out_~{chunk_id} ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading chunk. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp chunk_hwe_~{chunk_id} ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading chunk. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
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
