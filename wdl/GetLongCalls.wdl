version 1.0


# Given a cohort VCF, returns a VCF that contains only long calls.
#
workflow GetLongCalls {
    input {
        File cohort_vcf_gz
        File cohort_tbi
        Int min_sv_length = 0
    }
    parameter_meta {
    }

    call GetLongCallsImpl {
        input:
            cohort_vcf_gz = cohort_vcf_gz,
            cohort_tbi = cohort_tbi,
            min_sv_length = min_sv_length
    }
    
    output {
        File out_vcf_gz = GetLongCallsImpl.out_vcf_gz
        File out_tbi = GetLongCallsImpl.out_tbi
    }
}


#
task GetLongCallsImpl {
    input {
        File cohort_vcf_gz
        File cohort_tbi
        Int min_sv_length
        
        Int n_cpu = 64
        Int ram_size_gb = 64
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 5*ceil(size(cohort_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "SVLEN>=~{min_sv_length} || SVLEN<=-~{min_sv_length}" --output-type z ~{cohort_vcf_gz} > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
        bcftools index --nrecords out.vcf.gz.tbi
    >>>
    
    output {
        File out_vcf_gz = work_dir + "/out.vcf.gz"
        File out_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
