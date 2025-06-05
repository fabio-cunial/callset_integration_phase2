version 1.0


#
workflow DownloadAssembly {
    input {
        String sample_id
        String hap1_address
        String hap2_address
        String billing_project = "broad-dsp-lrma"
        String remote_output_dir
    }
    parameter_meta {
    }
    
    call DownloadImpl {
        input:
            sample_id = sample_id,
            hap1_address = hap1_address,
            hap2_address = hap2_address,
            billing_project = billing_project,
            remote_output_dir = remote_output_dir
    }
    
    output {
    }
}


task DownloadImpl {
    input {
        String sample_id
        String hap1_address
        String hap2_address
        String remote_output_dir
        String billing_project

        Int n_cores = 4
        Int mem_gb = 4
        Int disk_size_gb = 10
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        FILE_NAME=$(basename ~{hap1_address})
        FILE_NAME="~{sample_id}_hap1_${FILE_NAME}"
        echo "~{remote_output_dir}/${FILE_NAME}" > hap1.txt
        if [[ ~{hap1_address} == gs://* ]]; then
            gsutil ${GSUTIL_UPLOAD_THRESHOLD} -u ~{billing_project} -m cp ~{hap1_address} ~{remote_output_dir}/${FILE_NAME} &
        elif [[ ~{hap1_address} == s3://* ]]; then
            aws s3 --no-sign-request cp ~{hap1_address} - | gsutil cp - ~{remote_output_dir}/${FILE_NAME} &
        else
            curl ~{hap1_address} | gsutil cp - ~{remote_output_dir}/${FILE_NAME} &
        fi
        
        FILE_NAME=$(basename ~{hap2_address})
        FILE_NAME="~{sample_id}_hap2_${FILE_NAME}"
        echo "~{remote_output_dir}/${FILE_NAME}" > hap2.txt
        if [[ ~{hap2_address} == gs://* ]]; then
            gsutil ${GSUTIL_UPLOAD_THRESHOLD} -u ~{billing_project} -m cp ~{hap2_address} ~{remote_output_dir}/${FILE_NAME} &
        elif [[ ~{hap2_address} == s3://* ]]; then
            aws s3 --no-sign-request cp ~{hap2_address} - | gsutil cp - ~{remote_output_dir}/${FILE_NAME} &
        else
            curl ~{hap2_address} | gsutil cp - ~{remote_output_dir}/${FILE_NAME} &
        fi
        
        wait
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
