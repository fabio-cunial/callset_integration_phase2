version 1.0


#
workflow GetCompositeSvs {
    input {
        File cohort_vcf_gz
        File cohort_tbi
        
        File tandem_bed
        File reference_fai
        
        Int max_distance
        Int min_calls
    }
    parameter_meta {
    }
    
    call ExcludeTRs {
        input:
            cohort_vcf_gz = cohort_vcf_gz,
            cohort_tbi = cohort_tbi,
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    call Impl {
        input:
            cohort_vcf_gz = ExcludeTRs.filtered_vcf,
            cohort_tbi = ExcludeTRs.filtered_tbi,
            max_distance = max_distance,
            min_calls = min_calls
    }
    
    output {
        File out_txt = Impl.out_txt
    }
}


# Performance with 16 cores and 32GB of RAM on the 07 BI v1 VCF:
#
# COMMAND           CPU     RAM     TIME
# bcftools view     
#
#
task ExcludeTRs {
    input {
        File cohort_vcf_gz
        File cohort_tbi
        
        File tandem_bed
        File reference_fai
        
        Int n_cpu = 16
        Int ram_size_gb = 32
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(tandem_bed,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        # Complementing the TR track
        ${TIME_COMMAND} bedtools sort -i ~{tandem_bed} -faidx ~{reference_fai} > sorted.bed
        ${TIME_COMMAND} bedtools complement -i sorted.bed -L -g ~{reference_fai} > complement.bed
        
        # Subsetting the VCF
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --targets complement.bed --targets-overlap pos --output-type z ~{cohort_vcf_gz} > filtered.vcf.gz
        tabix -f filtered.vcf.gz
    >>>
    
    output {
        File filtered_vcf = "filtered.vcf.gz"
        File filtered_tbi = "filtered.vcf.gz.tbi"
        File complement_bed = "complement.bed"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Performance with 16 cores and 32GB of RAM on a 66x, 158GB BAM:
#
# COMMAND           CPU     RAM     TIME
# GetCompositeSvs   
#
task Impl {
    input {
        File cohort_vcf_gz
        File cohort_tbi
        
        Int max_distance
        Int min_calls

        Int n_cores = 2
        Int mem_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*( ceil(size(cohort_vcf_gz,"GB")) )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        EFFECTIVE_RAM_GB=$(( ~{mem_gb} - 2 ))
        
                
        bcftools index --nrecords ~{cohort_tbi}
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G GetCompositeSvs ~{cohort_vcf_gz} ~{max_distance} ~{min_calls} > out.txt
        ls -laht
    >>>
    
    output {
        File out_txt = "out.txt"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 2
    }
}
