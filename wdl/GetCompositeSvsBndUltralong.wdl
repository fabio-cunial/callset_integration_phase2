version 1.0


# Finds clusters of ultralong calls and BND calls that occur in the same sample.
#
workflow GetCompositeSvsBndUltralong {
    input {
        File bnd_bcf
        File ultralong_bcf
        File ultralong_csi
        
        Int max_distance = 10000

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        bnd_bcf: "Must be canonized"
    }
    
    call Concat {
        input:
            bnd_bcf = bnd_bcf,
            ultralong_bcf = ultralong_bcf,
            ultralong_csi = ultralong_csi,
            docker_image = docker_image
    }
    call Impl {
        input:
            cohort_0_vcf_gz = Concat.out_0_vcf_gz,
            cohort_0_tbi = Concat.out_0_tbi,
            cohort_1_vcf_gz = Concat.out_1_vcf_gz,
            cohort_1_tbi = Concat.out_1_tbi,
            max_distance = max_distance,
            docker_image = docker_image
    }
    
    output {
        File out_0_bed = Impl.out_0_bed
        File out_1_bed = Impl.out_1_bed
    }
}


# Creates a symmetrized duplicate of every BND record in the BND VCF and merges
# it with the ultralong VCF.
#
# COMMAND           CPU     RAM     TIME
# bcftools concat   
#
task Concat {
    input {
        File bnd_bcf
        File ultralong_bcf
        File ultralong_csi
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil( size(bnd_bcf,"GB") + size(ultralong_bcf,"GB") )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_MB=$(( ~{ram_size_gb}*1024 - 500 ))

        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --output-type z ~{bnd_bcf} --output bnd.vcf.gz
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M BndSymmetrize bnd.vcf.gz 0 | bcftools sort --output-type b --output bnd_symmetrized_0.bcf
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M BndSymmetrize bnd.vcf.gz 1 | bcftools sort --output-type b --output bnd_symmetrized_1.bcf
        ${TIME_COMMAND} bcftools index -f bnd_symmetrized_0.bcf
        ${TIME_COMMAND} bcftools index -f bnd_symmetrized_1.bcf
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --output-type z bnd_symmetrized_0.bcf ~{ultralong_bcf} --output concat_0.vcf.gz
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --output-type z bnd_symmetrized_1.bcf ~{ultralong_bcf} --output concat_1.vcf.gz
        ${TIME_COMMAND} tabix -f concat_0.vcf.gz
        ${TIME_COMMAND} tabix -f concat_1.vcf.gz
    >>>
    
    output {
        File out_0_vcf_gz = "concat_0.vcf.gz"
        File out_0_tbi = "concat_0.vcf.gz.tbi"
        File out_1_vcf_gz = "concat_1.vcf.gz"
        File out_1_tbi = "concat_1.vcf.gz.tbi"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# COMMAND                       CPU     RAM     TIME
# GetCompositeSvsPrime   
#
task Impl {
    input {
        File cohort_0_vcf_gz
        File cohort_0_tbi
        File cohort_1_vcf_gz
        File cohort_1_tbi
        
        Int max_distance

        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(cohort_0_vcf_gz,"GB")) + ceil(size(cohort_1_vcf_gz,"GB")) )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_MB=$(( ~{ram_size_gb}*1024 - 500 ))
        
        bcftools index --nrecords ~{cohort_0_tbi}
        bcftools index --nrecords ~{cohort_1_tbi}
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M GetCompositeSvsPrime ~{cohort_0_vcf_gz} 0 ~{max_distance} 2 1 out_0.bed
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M GetCompositeSvsPrime ~{cohort_1_vcf_gz} 1 ~{max_distance} 2 1 out_1.bed

        # Sorting by: nCalls > nTypes > maxLength
        ${TIME_COMMAND} sort -k4,4nr -k5,5nr -k11,11nr out_0.bed > out_0_sorted.bed
        mv out_0_sorted.bed out_0.bed
        ${TIME_COMMAND} sort -k4,4nr -k5,5nr -k11,11nr out_1.bed > out_1_sorted.bed
        mv out_1_sorted.bed out_1.bed
        
        ls -laht
    >>>
    
    output {
        File out_0_bed = "out_0.bed"
        File out_1_bed = "out_1.bed"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
