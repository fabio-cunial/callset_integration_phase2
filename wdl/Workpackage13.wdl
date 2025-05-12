version 1.0


# 
#
workflow Workpackage13 {
    input {
        String remote_indir_header
        String remote_indir_chunks
        String remote_outdir
        
        Int n_cpu = 8
        Int ram_size_gb = 8
        Int disk_size_gb = 3000
    }
    parameter_meta {
    }
    
    call Workpackage13Impl {
        input:
            remote_indir_header = remote_indir_header,
            remote_indir_chunks = remote_indir_chunks,
            remote_outdir = remote_outdir,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}


# Performance on 10'070 samples, 15x, GRCh38, stringent (_S) and lenient (_L):
#
# TOOL                      CPU_S   RAM_S   TIME_S  CPU_L   RAM_L   TIME_L
# cat | bgzip                               1h                      2h     
# tabix                                     1h                      1.5h
#
task Workpackage13Impl {
    input {
        String remote_indir_header
        String remote_indir_chunks
        String remote_outdir
        
        Int n_cpu
        Int ram_size_gb
        Int disk_size_gb
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
        
        
        # Localizing header and chunks
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir_header}/header.txt . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading header. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir_header}/fields_all.txt . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading header. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir_chunks}/'chunk_out_*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading chunks. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        # Concatenating
        CHUNK_FILES=$(ls chunk_out_* | sort --version-sort | tr '\n' ' ')
        date
        cat header.txt fields_all.txt ${CHUNK_FILES} | bgzip -@ ${N_THREADS} --compress-level 2 > final.vcf.gz
        date
        tabix -f final.vcf.gz
        date
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp final.vcf.'gz*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading VCF. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        rm -f final.vcf.gz* header.txt fields_all.txt ${CHUNK_FILES}
        
        # Localizing HWE counts
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir_chunks}/'chunk_hwe_*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading HWE counts. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        CHUNK_FILES=$(ls chunk_hwe_* | sort --version-sort | tr '\n' ' ')
        ${TIME_COMMAND} cat ${CHUNK_FILES} > gt_counts.txt
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp gt_counts.txt ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading HWE counts. Trying again..."
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
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
