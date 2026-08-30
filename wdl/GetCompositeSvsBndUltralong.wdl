version 1.0


# Finds clusters of ultralong calls and intra-chromosomal BND calls that occur 
# in the same sample, treating BNDs as intervals (inter-chromosomal BNDs are 
# discarded).
#
workflow GetCompositeSvsBndUltralong {
    input {
        File bnd_bcf
        File bnd_csi
        File ultralong_bcf
        File ultralong_csi
        
        Int max_distance

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        bnd_bcf: "Must be canonized"
    }
    
    call Concat {
        input:
            bnd_bcf = bnd_bcf,
            bnd_csi = bnd_csi,
            ultralong_bcf = ultralong_bcf,
            ultralong_csi = ultralong_csi,
            docker_image = docker_image
    }
    call Impl {
        input:
            cohort_vcf_gz = Concat.out_vcf_gz,
            cohort_tbi = Concat.out_tbi,
            max_distance = max_distance,
            docker_image = docker_image
    }
    
    output {
        File out_bed = Impl.out_bed
    }
}


# COMMAND           CPU     RAM     TIME
# bcftools concat   
#
task Concat {
    input {
        File bnd_bcf
        File bnd_csi
        File ultralong_bcf
        File ultralong_csi
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil( size(bnd_bcf,"GB") + size(ultralong_bcf,"GB") )
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --output-type z --output concat.vcf.gz ~{bnd_bcf} ~{ultralong_bcf}
        tabix -f concat.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = "concat.vcf.gz"
        File out_tbi = "concat.vcf.gz.tbi"
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
        File cohort_vcf_gz
        File cohort_tbi
        
        Int max_distance

        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(cohort_vcf_gz,"GB")) )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_MB=$(( ~{ram_size_gb}*1024 - 500 ))
        
        bcftools index --nrecords ~{cohort_tbi}
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M GetCompositeSvsPrime ~{cohort_vcf_gz} ~{max_distance} 2 1 out.bed

        # Sorting by: nCalls > nTypes > maxLength
        ${TIME_COMMAND} sort -k4,4nr -k5,5nr -k11,11nr out.bed > out_sorted.bed
        mv out_sorted.bed out.bed
        
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
