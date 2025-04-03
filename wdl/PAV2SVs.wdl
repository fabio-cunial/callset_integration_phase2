version 1.0


#
workflow PAV2SVs {
    input {
        File input_vcf_gz
        Int min_sv_length
        Int compression_level = 1
    }
    parameter_meta {
    }
    
    call PAV2SVsImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            min_sv_length = min_sv_length,
            compression_level = compression_level
    }
    
    output {
        File sv_vcf_gz = PAV2SVsImpl.sv_vcf_gz
        File sv_tbi = PAV2SVsImpl.sv_tbi
    }
}


task PAV2SVsImpl {
    input {
        File input_vcf_gz
        Int min_sv_length
        Int compression_level
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_vcf_gz,"GB"))
    Int ram_size_gb = 4
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
        
        # Removing multiallelic records, if any.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type v ~{input_vcf_gz} > input.vcf
        
        # Keeping only long calls
        ${TIME_COMMAND} java -cp ~{docker_dir} PAV2SVs input.vcf ~{min_sv_length} sv.vcf snp.vcf
        ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} sv.vcf
        tabix -f sv.vcf.gz
    >>>
    
    output {
        File sv_vcf_gz = work_dir + "/sv.vcf.gz"
        File sv_tbi = work_dir + "/sv.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: 4
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
