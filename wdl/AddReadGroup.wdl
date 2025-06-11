version 1.0


# Adds an RG to a PacBio BAM.
#
workflow AddReadGroup {
    input {
        String sample_id
        File input_bam
        File input_bai

        Int n_cpus = 8
        Int ram_size_gb = 4
    }
    parameter_meta {
    }
    
    call AddReadGroupImpl {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb
    }
    
    output {
        File output_bam = AddReadGroupImpl.output_bam
        File output_bai = AddReadGroupImpl.output_bai
    }
}


# Performance with 64 cores and 64 GB of RAM.
#
# TASK                      % CPU       RAM     TIME
# samtools addreplacerg     500%        50M     3m
#
task AddReadGroupImpl {
    input {
        String sample_id
        File input_bam
        File input_bai

        Int n_cpus
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_bam, "GB"))
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
        FAKE_RG="@RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:~{sample_id}\tPM:SEQUEL"
        
        ${TIME_COMMAND} samtools addreplacerg -@ ${N_THREADS} -r ${FAKE_RG} -m overwrite_all --no-PG -o ~{sample_id}_rg.bam ~{input_bam}
        ${TIME_COMMAND} samtools index -@ ${N_THREADS} ~{sample_id}_rg.bam
    >>>
    
    output {
        File output_bam = work_dir + "/" + sample_id + "_rg.bam"
        File output_bai = work_dir + "/" + sample_id + "_rg.bam.bai"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
