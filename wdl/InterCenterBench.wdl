version 1.0


#
workflow InterCenterBench {
    input {
        File center1_vcf_gz
        File center1_tbi
        File center2_vcf_gz
        File center2_tbi
        File center3_vcf_gz
        File center3_tbi
    }
    parameter_meta {
    }

    call Bench as Bench1 {
        input:
            center1_vcf_gz = center1_vcf_gz,
            center1_tbi = center1_tbi,
            center2_vcf_gz = center2_vcf_gz,
            center2_tbi = center2_tbi
    }
    call Bench as Bench2 {
        input:
            center1_vcf_gz = center1_vcf_gz,
            center1_tbi = center1_tbi,
            center2_vcf_gz = center3_vcf_gz,
            center2_tbi = center3_tbi
    }
    call Bench as Bench3 {
        input:
            center1_vcf_gz = center2_vcf_gz,
            center1_tbi = center2_tbi,
            center2_vcf_gz = center3_vcf_gz,
            center2_tbi = center3_tbi
    }
    
    output {
    }
}


#
task Bench {
    input {
        File center1_vcf_gz
        File center1_tbi
        File center2_vcf_gz
        File center2_tbi
        
        Int n_cpu = 4
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
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
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        ${TIME_COMMAND} truvari bench -b ~{center1_vcf_gz} -c ~{center2_vcf_gz} -o ./truvari/
        mv ./truvari/summary.json .
    >>>
    
    output {
        File summary = work_dir + "/summary.json"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}
