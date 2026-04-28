version 1.0


# 
#
workflow SvimAsm {
    input {
        String sample_id
        
        File hap1_bam
        File hap1_bai
        File hap2_bam
        File hap2_bai

        File reference_fa
        File reference_fai
        
        String remote_output_dir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            sample_id = sample_id,
            hap1_bam = hap1_bam,
            hap1_bai = hap1_bai,
            hap2_bam = hap2_bam,
            hap2_bai = hap2_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            remote_output_dir = remote_output_dir,
            docker_image = docker_image
    }
    
    output {
    }
}


#
task Impl {
    input {
        String sample_id
        
        File hap1_bam
        File hap1_bai
        File hap2_bam
        File hap2_bai

        File reference_fa
        File reference_fai
        
        String remote_output_dir
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 32
        Int disk_size_gb = 20
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        
        svim-asm --version 2>&1 || echo "1"

        ${TIME_COMMAND} svim-asm diploid ./svim/ ~{hap1_bam} ~{hap2_bam} ~{reference_fa}
        ${TIME_COMMAND} bcftools view --output-type z ./svim/variants.vcf --output ./~{sample_id}_svim.vcf.gz
        bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svim.vcf.gz

        gcloud storage cp ~{sample_id}_svim.vcf.'gz*' ~{remote_output_dir}/
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
