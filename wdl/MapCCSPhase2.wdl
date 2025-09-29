version 1.0


# Reproduces the AoU Phase 2 production pipeline at:
#
# https://app.terra.bio/#workspaces/allofus-drc-wgs-lr-prod/AoU_DRC_WGS_LongReads_BI_CCS/workflows/allofus-drc-wgs-lr-prod/IngestSMRTcode.GRCh38
# https://github.com/broadinstitute/long-read-pipelines/blob/sh_ingest_singlerg/wdl/tasks/Alignment/AlignHiFiUBAM.wdl
# https://github.com/broadinstitute/long-read-pipelines/blob/sh_ingest_singlerg/wdl/tasks/Utility/PBUtils.wdl
#
workflow MapCCSPhase2 {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        Int n_cpus
        Int ram_size_gb
        Int disk_gb
        String docker = "us.gcr.io/broad-dsp-lrma/lr-smrttools:12.0.0.176214"
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
            disk_gb = disk_gb,
            docker = docker
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
        String docker
    }
    parameter_meta {
    }
    
    Int disk_size_gb = ceil(size(reads_fastq_gz, "GB")) + ceil(size(reference_fa, "GB")) + disk_gb
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Remark: pbmm2 automatically uses all cores.
        pbmm2 align --preset CCS --sort --sample ~{sample_id} ~{reference_fa} ~{reads_fastq_gz} out.bam
        samtools calmd -@ ${N_THREADS} --no-PG -b out.bam ~{reference_fa} > ~{sample_id}.bam
        samtools index -@ ${N_THREADS} ~{sample_id}.bam
    >>>
    
    output {
        File output_bam = work_dir + "/" + sample_id + ".bam"
        File output_bai = work_dir + "/" + sample_id + ".bam.bai"
    }
    runtime {
        docker: docker
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
