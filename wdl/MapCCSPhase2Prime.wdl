version 1.0


# Like `MapCCSPhase2.wdl`, but uses minimap2 directly to avoid a pbmm2 crash
# with some FASTQs.
#
workflow MapCCSPhase2Prime {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        Int n_cpus
        Int ram_size_gb
        Int disk_gb
    }
    parameter_meta {
    }
    
    call MapCCSImpl {
        input:
            sample_id = sample_id,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            reads_fastq_gz = reads_fastq_gz,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb,
            disk_gb = disk_gb
    }
    
    output {
        File output_bam = MapCCSImpl.output_bam
        File output_bai = MapCCSImpl.output_bai
    }
}


#
task MapCCSImpl {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        Int n_cpus
        Int ram_size_gb
        Int disk_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = ceil(size(reads_fastq_gz, "GB")) + ceil(size(reference_fa, "GB")) + disk_gb
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
        
        # From `https://github.com/PacificBiosciences/pbmm2#
        # what-are-parameter-sets-and-how-can-i-override-them`
        PBMM2_CCS_PARAMS="-k 19 -w 19 -O 6,26 -E 2,1 -A 1 -B 4 -z 400,50 -r 2000 -g 5000"
        ${TIME_COMMAND} ~{docker_dir}/minimap2/minimap2 ${PBMM2_CCS_PARAMS} -ayYL --MD --eqx --cs -t ${N_THREADS} -K4G ~{reference_fa} ~{reads_fastq_gz} > out.sam
        ${TIME_COMMAND} samtools sort -@4 -m4G --no-PG -o out.bam out.sam
        rm -f out.sam
        ${TIME_COMMAND} samtools calmd -@ ${N_THREADS} --no-PG -b out.bam ~{reference_fa} > ~{sample_id}.bam
        ${TIME_COMMAND} samtools index -@ ${N_THREADS} ~{sample_id}.bam
    >>>
    
    output {
        File output_bam = work_dir + "/" + sample_id + ".bam"
        File output_bai = work_dir + "/" + sample_id + ".bam.bai"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
