version 1.0


#
workflow DownloadAssembly {
    input {
        String hap1_address
        String hap2_address
        String billing_project
        String remote_output_dir
    }
    parameter_meta {
    }
    
    call DownloadImpl {
        input:
            hap1_address = hap1_address,
            hap2_address = hap2_address,
            billing_project = billing_project,
            remote_output_dir = remote_output_dir
    }
    
    output {
        String hap1_downloaded_address = DownloadImpl.hap1_downloaded_address
        String hap2_downloaded_address = DownloadImpl.hap2_downloaded_address
    }
}


task DownloadImpl {
    input {
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
        
        
        if [[ ~{hap1_address} == gs://* ]]; then
            gsutil ${GSUTIL_UPLOAD_THRESHOLD} -u ~{billing_project} -m cp ~{hap1_address} ~{remote_output_dir}/ &
        elif [[ ~{hap1_address} == s3://* ]]; then
            aws s3 --no-sign-request cp ~{hap1_address} - | gsutil cp - ~{remote_output_dir}/ &
        else
            curl ~{hap1_address} | gsutil cp - ~{remote_output_dir}/ &
        fi
        if [[ ~{hap2_address} == gs://* ]]; then
            gsutil ${GSUTIL_UPLOAD_THRESHOLD} -u ~{billing_project} -m cp ~{hap2_address} ~{remote_output_dir}/ &
        elif [[ ~{hap2_address} == s3://* ]]; then
            aws s3 --no-sign-request cp ~{hap2_address} - | gsutil cp - ~{remote_output_dir}/ &
        else
            curl ~{hap2_address} | gsutil cp - ~{remote_output_dir}/ &
        fi
        wait
        echo "~{remote_output_dir}/$(basename ~{hap1_address})" > hap1.txt
        echo "~{remote_output_dir}/$(basename ~{hap2_address})" > hap2.txt
    >>>
    
    output {
        String hap1_downloaded_address = read_string(work_dir+"hap1.txt")
        String hap2_downloaded_address = read_string(work_dir+"hap2.txt")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
