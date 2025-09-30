version 1.0


# Remark: the workflow normalizes multiallelic sites and then filters SVs.
#
workflow Dipcall2SVs {
    input {
        String sample_id
        File input_vcf_gz
        Int min_sv_length
    }
    parameter_meta {
    }
    
    call Dipcall2SVsImpl {
        input:
            sample_id = sample_id,
            input_vcf_gz = input_vcf_gz,
            min_sv_length = min_sv_length
    }
    
    output {
        File sv_vcf_gz = Dipcall2SVsImpl.sv_vcf_gz
        File sv_tbi = Dipcall2SVsImpl.sv_tbi
    }
}


task Dipcall2SVsImpl {
    input {
        String sample_id
        File input_vcf_gz
        Int min_sv_length
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_vcf_gz,"GB"))
    Int ram_size_gb = 4
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        bcftools norm --multiallelics - --output-type v ~{input_vcf_gz} > input.vcf
        rm -f ~{input_vcf_gz}
        ${TIME_COMMAND} java -cp ~{docker_dir} Dipcall2VCF input.vcf ~{min_sv_length} ~{sample_id}_sv.vcf
        ${TIME_COMMAND} bgzip ~{sample_id}_sv.vcf
        tabix -f ~{sample_id}_sv.vcf.gz
    >>>
    
    output {
        File sv_vcf_gz = sample_id + "_sv.vcf.gz"
        File sv_tbi = sample_id + "_sv.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: 1
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
