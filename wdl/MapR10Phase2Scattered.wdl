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
workflow MapR10Phase2Scattered {
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
            reads_fastq_gz = reads_fastq_gz
    }
    scatter (i in range(length(RemoveDuplicatedReads.chunks_fastq_gz))) {
        call Minimap2 {
            input:
                sample_id = sample_id,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reads_fastq_gz = RemoveDuplicatedReads.chunks_fastq_gz[i]
        }
        call Sam2Bam {
            input:
                sample_id = sample_id,
                input_bam = Minimap2.output_bam,
                reference_fa = reference_fa,
                reference_fai = reference_fai
        }
    }
    call MergeChunks {
        input:
            sample_id = sample_id,
            input_bam = Sam2Bam.output_bam,
            input_bai = Sam2Bam.output_bai
    }
    
    
    output {
        File output_bam = MergeChunks.output_bam
        File output_bai = MergeChunks.output_bai
    }
}


# ONT readsets may contain exact duplicates. 
#
# Remark: this step does not reproduce the AoU pipeline exactly, since we work
# with FASTQs whereas the pipeline worked with BAMs.
#
# Performance with 8 cores and 8GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# seqkit rmdup              200%        1.1G    1h
# seqkit split              190%        100M    1h
# 
task RemoveDuplicatedReads {
    input {
        File reads_fastq_gz
        Int n_chunks = 8
        
        Int n_cpus = 4
        Int ram_size_gb = 4
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(reads_fastq_gz, "GB"))
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
        rm -f ~{reads_fastq_gz}
        ${TIME_COMMAND} ~{docker_dir}/seqkit split2 --by-part ~{n_chunks} --by-part-prefix chunk --extension .gz --out-dir . out.fastq.gz
    >>>
    
    output {
        Array[File] chunks_fastq_gz = glob(work_dir+"/chunk*")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}



# Performance on one 7.5x chunk with 64 cores and 128 GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# minimap2                  ???         70G     1h
# samtools view             ???         ???     ~1h?
#
task Minimap2 {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        
        Int n_cpus = 64
        Int ram_size_gb = 128
        Int disk_gb = 500
        
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
        minimap2 -t ${N_THREADS} -K4G -x map-ont -ayYL --MD --eqx --cs ~{reference_fa} -R ${FAKE_RG} ~{reads_fastq_gz} > out.sam
        ls -laht
        df -h
        samtools view -@ ${N_THREADS} -b out.sam > out.bam
        ls -laht
        df -h
    >>>
    
    output {
        File output_bam = work_dir + "/out.bam"
    }
    runtime {
        docker: docker
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}


# Performance on one 7.5x chunk with 64 cores and 128 GB of RAM:
#
# TASK                      CPU%        RAM     TIME
# samtools sort             500%        51G     30m
# samtools calmd            600%        500M    15m
# samtools index            1900%       70M     1m
#
task Sam2Bam {
    input {
        String sample_id
        File input_bam
        File reference_fa
        File reference_fai
        
        Int n_cpus = 8
        Int ram_size_gb = 64
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(input_bam, "GB"))
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
        ls -laht
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


# Performance on 30x with 64 cores and 128 GB of RAM:
#
# TASK                      CPU %       RAM     TIME
# samtools merge            500%        60M     1h
# samtools index            100%        150M    1h
#
task MergeChunks {
    input {
        String sample_id
        Array[File] input_bam
        Array[File] input_bai
        
        Int n_cpus = 8
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(input_bam, "GB"))
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
        while read FILE; do
            echo ${FILE} >> list.txt
        done < ~{write_lines(input_bam)}
        ${TIME_COMMAND} samtools merge -@ ${N_THREADS} -b list.txt -p -c --no-PG -o ~{sample_id}_merged.bam
        ${TIME_COMMAND} samtools index -@ ${N_THREADS} ~{sample_id}_merged.bam
        ls -laht
        df -h
    >>>
    
    output {
        File output_bam = work_dir + "/" + sample_id + "_merged.bam"
        File output_bai = work_dir + "/" + sample_id + "_merged.bam.bai"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}
