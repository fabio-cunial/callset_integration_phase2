version 1.0


#
workflow GetCompositeSvs {
    input {
        File cohort_bcf
        File cohort_csi
        
        File tandem_bed
        File reference_fai
        Int slop_bp
        
        Int max_distance
        Int min_calls
        Int min_sv_length

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
    }
    
    call ExcludeTRs {
        input:
            cohort_bcf = cohort_bcf,
            cohort_csi = cohort_csi,
            tandem_bed = tandem_bed,
            reference_fai = reference_fai,
            slop_bp = slop_bp,
            docker_image = docker_image
    }
    call Impl {
        input:
            cohort_vcf_gz = ExcludeTRs.filtered_vcf_gz,
            cohort_tbi = ExcludeTRs.filtered_tbi,
            max_distance = max_distance,
            min_calls = min_calls,
            min_sv_length = min_sv_length,
            docker_image = docker_image
    }
    
    output {
        File out_bed = Impl.out_bed
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
        File cohort_bcf
        File cohort_csi
        
        File tandem_bed
        File reference_fai
        Int slop_bp
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil( size(cohort_bcf,"GB") )
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        # Complementing the TR track
        ${TIME_COMMAND} bedtools sort -i ~{tandem_bed} -faidx ~{reference_fai} > sorted.bed
        ${TIME_COMMAND} bedtools slop -i sorted.bed -g ~{reference_fai} -b ~{slop_bp} > sorted_slop.bed
        ${TIME_COMMAND} bedtools complement -i sorted_slop.bed -L -g ~{reference_fai} > complement.bed
        
        # Subsetting the VCF
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --targets-file complement.bed --targets-overlap pos --output-type z ~{cohort_bcf} --output filtered.vcf.gz
        tabix -f filtered.vcf.gz
    >>>
    
    output {
        File filtered_vcf_gz = "filtered.vcf.gz"
        File filtered_tbi = "filtered.vcf.gz.tbi"
        File complement_bed = "complement.bed"
    }
    runtime {
        docker: docker_image
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
        Int min_sv_length

        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        
        bcftools index --nrecords ~{cohort_tbi}
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G GetCompositeSvsPrime ~{cohort_vcf_gz} ~{max_distance} ~{min_calls} ~{min_sv_length} out.bed
        ls -laht
    >>>
    
    output {
        File out_bed = "out.bed"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
