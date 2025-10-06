version 1.0


#
workflow SubsetToAncestry {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File sample_ids
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            vcf_gz = truvari_collapsed_vcf_gz,
            tbi = truvari_collapsed_tbi,
            sample_ids = sample_ids
    }
    
    output {
        File out_vcf_gz = Impl.out_vcf_gz
        File out_tbi = Impl.out_tbi
        File selected_samples = Impl.selected_samples
    }
}


# Remark: the task removes records that do not occur in the selected samples.
#
task Impl {
    input {
        File vcf_gz
        File tbi
        File sample_ids
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        sample_ids: "One sample per line. Not necessarily sorted."
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 3*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        cut -f 1 ~{sample_ids} | sort > desired_samples.txt
        bcftools view --header-only ~{vcf_gz} | tail -n 1 | tr '\t' '\n' | tail -n +10 | sort > present_samples.txt
        comm -1 -2 desired_samples.txt present_samples.txt > selected_samples.txt
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file selected_samples.txt --output-type z ~{vcf_gz} > subset.vcf.gz
        ${TIME_COMMAND} tabix -f subset.vcf.gz
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z subset.vcf.gz > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>

    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
        File selected_samples = "selected_samples.txt"
    }

    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
