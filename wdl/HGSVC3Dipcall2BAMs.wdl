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
        File hap1_bai = HGSVC3Dipcall2BAMsImpl.hap1_bai
        File hap2_bam = HGSVC3Dipcall2BAMsImpl.hap2_bam
        File hap2_bai = HGSVC3Dipcall2BAMsImpl.hap2_bai
    }
}


# Example content of dipcall's `.tar.gz` output:
#
# HG00096.fna.dip.bed
# HG00096.fna.dip.vcf.gz
# HG00096.fna.hap1.bam
# HG00096.fna.hap1.bed
# HG00096.fna.hap1.paf.gz
# HG00096.fna.hap1.paf.gz.log
# HG00096.fna.hap1.sam.gz.log
# HG00096.fna.hap1.var.gz
# HG00096.fna.hap1.var.gz.vst
# HG00096.fna.hap2.bam
# HG00096.fna.hap2.bed
# HG00096.fna.hap2.paf.gz
# HG00096.fna.hap2.paf.gz.log
# HG00096.fna.hap2.sam.gz.log
# HG00096.fna.hap2.var.gz
# HG00096.fna.hap2.var.gz.vst
# HG00096.fna.pair.vcf.gz
# HG00096.fna.tmp.be
#
task HGSVC3Dipcall2BAMsImpl {
    input {
        String sample_id
        File tar_gz
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        DIR_NAME=$(basename ~{tar_gz})
        DIR_NAME=${DIR_NAME%.tar.gz}
        tar -xzf ~{tar_gz}
        mv ${DIR_NAME}/*.hap1.bam ~{sample_id}_hap1.bam
        samtools index --threads ${N_THREADS} --bai ~{sample_id}_hap1.bam
        mv ${DIR_NAME}/*.hap2.bam ~{sample_id}_hap2.bam
        samtools index --threads ${N_THREADS} --bai ~{sample_id}_hap2.bam
    >>>
    
    output {
        File hap1_bam = work_dir + "/" + sample_id + "_hap1.bam"
        File hap1_bai = work_dir + "/" + sample_id + "_hap1.bam.bai"
        File hap2_bam = work_dir + "/" + sample_id + "_hap2.bam"
        File hap2_bai = work_dir + "/" + sample_id + "_hap2.bam.bai"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: 16
        memory: "16GB"
        disks: "local-disk 500 HDD"
        preemptible: 0
    }
}
