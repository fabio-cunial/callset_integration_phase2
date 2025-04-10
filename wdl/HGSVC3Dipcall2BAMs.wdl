version 1.0


# 
#
workflow HGSVC3Dipcall2BAMs {
    input {
        String sample_id
        File tar_gz
    }
    parameter_meta {
    }
    
    call HGSVC3Dipcall2BAMsImpl {
        input:
            sample_id = sample_id,
            tar_gz = tar_gz
    }
    
    output {
        File hap1_bam = HGSVC3Dipcall2BAMsImpl.hap1_bam
        File hap2_bam = HGSVC3Dipcall2BAMsImpl.hap2_bam
    }
}


#
task HGSVC3Dipcall2BAMsImpl {
    input {
        String sample_id
        File tar_gz
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
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        DIR_NAME=$(basename ~{tar_gz})
        DIR_NAME=${DIR_NAME%.tar.gz}
        tar -xzf ~{tar_gz}
        mv ${DIR_NAME}/*.hap1.bam ~{sample_id}_hap1.bam
        mv ${DIR_NAME}/*.hap2.bam ~{sample_id}_hap2.bam
    >>>
    
    output {
        File hap1_bam = work_dir + "/" + sample_id + "_hap1.bam"
        File hap2_bam = work_dir + "/" + sample_id + "_hap2.bam"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: 1
        memory: "8GB"
        disks: "local-disk 500 HDD"
        preemptible: 0
    }
}
