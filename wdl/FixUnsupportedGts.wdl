version 1.0


# 
#
workflow FixUnsupportedGts {
    input {
        File ped_tsv
        
        File intersample_vcf_gz
        File intersample_tbi
        
        File script_java
        Int min_ad = 2
        Int ad_index = 7
    }
    parameter_meta {
    }
    
    call SubsetToSamples {
        input:
            cohort_vcf_gz = intersample_vcf_gz,
            cohort_tbi = intersample_tbi,
            ped_tsv = ped_tsv
    }
    call Impl {
        input:
            intersample_vcf_gz = SubsetToSamples.out_vcf_gz,
            intersample_tbi = SubsetToSamples.out_tbi,
            script_java = script_java,
            min_ad = min_ad,
            ad_index = ad_index
    }
    
    output {
        File fixed_vcf_gz = Impl.fixed_vcf_gz
        File fixed_tbi = Impl.fixed_tbi
    }
}



# Restricts the cohort VCF to records that occur in a given set of samples, to
# speed up the following steps.
#
task SubsetToSamples {
    input {
        File cohort_vcf_gz
        File cohort_tbi
        File ped_tsv
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(cohort_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        
        cut -f 2 ~{ped_tsv} > tmp.txt
        cut -f 3 ~{ped_tsv} >> tmp.txt
        cut -f 4 ~{ped_tsv} >> tmp.txt
        sort tmp.txt | uniq > list.txt
        rm -f tmp.txt
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file list.txt --output-type z ~{cohort_vcf_gz} > tmp1.vcf.gz
        ${TIME_COMMAND} tabix -f tmp1.vcf.gz
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1.vcf.gz > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}



#
task Impl {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        File script_java
        
        Int min_ad
        Int ad_index
        
        Int n_cpu = 4
        Int ram_size_gb = 16
        Int disk_size_gb = 500
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
        
        
        if [ ~{defined(script_java)} ]; then
            mv ~{script_java} ./FixUnsupportedGts.java
            javac FixUnsupportedGts.java
        else
            mv ~{docker_dir}/FixUnsupportedGts.java .
            javac FixUnsupportedGts.java
        fi
        ${TIME_COMMAND} ( java -Xmx$(( ~{ram_size_gb} - 2 ))G FixUnsupportedGts ~{intersample_vcf_gz} ~{min_ad} ~{ad_index} | bcftools view --output-type z > fixed.vcf.gz )
        ${TIME_COMMAND} tabix -f fixed.vcf.gz
        ls -laht
    >>>
    
    output {
        File fixed_vcf_gz = "fixed.vcf.gz"
        File fixed_tbi = "fixed.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_squish"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
