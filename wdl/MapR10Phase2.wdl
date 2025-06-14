version 1.0


# Reproduces the AoU Phase 2 production pipeline at:
#
# https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/ONTFlowcellWGSuBAM:sh_ont_fc?tab=info
#
# https://github.com/broadinstitute/long-read-pipelines/blob/sh_ont_fc/wdl/pipelines/ONT/Preprocessing/ONTFlowcellWGSuBAM.wdl
#
# https://github.com/broadinstitute/long-read-pipelines/blob/sh_ont_fc/wdl/tasks/Alignment/AlignONTWGSuBAM.wdl
#
# https://github.com/broadinstitute/long-read-pipelines/blob/sh_ont_fc/wdl/tasks/Alignment/AlignReads.wdl
#
workflow MapR10Phase2 {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        
    }
    parameter_meta {
    }
    
    call RemoveDuplicatedReads {
        input:
            sample_id = sample_id,
            reads_fastq_gz = reads_fastq_gz
    }
    call Minimap2 {
        input:
            sample_id = sample_id,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            reads_fastq_gz = RemoveDuplicatedReads.out_fastq_gz
    }
    call Sam2Bam {
        input:
            sample_id = sample_id,
            input_bam = Minimap2.output_bam,
            input_bai = Minimap2.output_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai
    }
    
    output {
        File output_bam = Sam2Bam.output_bam
        File output_bai = Sam2Bam.output_bai
    }
}


# ONT readsets may contain exact duplicates. 
#
# Remark: this step does not reproduce the AoU pipeline exactly, since we work
# with FASTQs whereas the pipeline worked with BAMs.
#
# Performance with 4 cores and 128GB of RAM.
#
# TASK                      % CPU       RAM     TIME
# seqkit rmdup              350%        1.1G    30m
# 
task RemoveDuplicatedReads {
    input {
        String sample_id
        File reads_fastq_gz
        
        Int n_cpus = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(reads_fastq_gz, "GB"))
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
        
        ${TIME_COMMAND} ~{docker_dir}/seqkit rmdup --by-seq -o out.fastq.gz ~{reads_fastq_gz}
    >>>
    
    output {
        File out_fastq_gz = work_dir + "/out.fastq.gz"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}



# Performance with 32 cores and 64 GB of RAM.
#
# TASK                      % CPU       RAM     TIME
# minimap2                  ???         51G     6.3h
# samtools view             
# minimap2 | samtools view  ???         73G     15h
#
task Minimap2 {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        
        Int n_cpus = 64
        Int ram_size_gb = 128
        Int disk_gb = 1000
        
        String docker = "us.gcr.io/broad-dsp-lrma/lr-minimap2:2.26-gcloud"
    }
    parameter_meta {
    }
    
    Int disk_size_gb = ceil(size(reads_fastq_gz, "GB")) + ceil(size(reference_fa, "GB")) + disk_gb
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        FAKE_RG="@RG\tID:default\tPL:ONT\tDS:READTYPE=UNKNOWN\tPU:default\tSM:~{sample_id}\tPM:ONT"
        
        ls -laht
        df -h
        minimap2 -ayYL --MD --eqx --cs -x map-ont -t ${N_THREADS} ~{reference_fa} -R ${FAKE_RG} ~{reads_fastq_gz} > out.sam
        # -K4G 
        ls -laht
        df -h
        samtools view -@ ${N_THREADS} -b out.sam > out.bam
        samtools index -@ ${N_THREADS} out.bam
        ls -laht
        df -h
    >>>
    
    output {
        File output_bam = work_dir + "/out.bam"
        File output_bai = work_dir + "/out.bam.bai"
    }
    runtime {
        docker: docker
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}



#
task Sam2Bam {
    input {
        String sample_id
        File input_bam
        File input_bai
        File reference_fa
        File reference_fai
        
        Int n_cpus = 64
        Int ram_size_gb = 128
        Int disk_gb = 1000
    }
    parameter_meta {
    }
    
    Int disk_size_gb = ceil(size(input_bam, "GB")) + disk_gb
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ls -laht
        df -h
        mkdir ./tmp/
        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} -T ./tmp/prefix --no-PG -O BAM ~{input_bam} > out.bam
        ls -laht
        df -h
        rm -f ~{input_bam}
        df -h
        ${TIME_COMMAND} samtools calmd -@ ${N_THREADS} --no-PG -b out.bam ~{reference_fa} > ~{sample_id}.bam
        ls -laht
        df -h
        ${TIME_COMMAND} samtools index -@ ${N_THREADS} ~{sample_id}.bam
        ls -laht
        df -h
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
        preemptible: 2
    }
}
