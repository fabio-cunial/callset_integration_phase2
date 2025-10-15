version 1.0


#
workflow Bam2Fastq {
    input {
        String sample_id
        String suffix
        
        File aligned_bam
        File aligned_bai
        
        String remote_output_dir
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            sample_id = sample_id,
            suffix = suffix,
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            remote_output_dir = remote_output_dir
    }
    
    output {
    }
}


# Performance with 16 cores and 32GB of RAM on a 66x, 158GB BAM:
#
# COMMAND           CPU     RAM     TIME
# samtools fastq    15%     15M     1h
#
task Impl {
    input {
        String sample_id
        String suffix
        
        File aligned_bam
        File aligned_bai
        
        String remote_output_dir

        Int n_cores = 2
        Int mem_gb = 4
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(aligned_bam,"GB")) )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        # Converting
        ${TIME_COMMAND} samtools fastq -@ ${N_THREADS} -n ~{aligned_bam} > ~{sample_id}_~{suffix}.fastq.gz
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{sample_id}_~{suffix}.fastq.gz ~{remote_output_dir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}
