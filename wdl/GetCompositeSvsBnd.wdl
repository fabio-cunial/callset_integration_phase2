version 1.0


# Finds clusters of inter-chromosomal BND breakpoints that occur in the same
# sample.
#
workflow GetCompositeSvsBnd {
    input {
        File bnd_bcf
        Int max_distance = 10000

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        bnd_bcf: "Must be canonized"
    }
    
    call Impl {
        input:
            bnd_bcf = bnd_bcf,
            max_distance = max_distance,
            docker_image = docker_image
    }
    
    output {
        File out_bed = Impl.out_bed
    }
}


# COMMAND                       CPU     RAM     TIME
# BndSymmetrize
# GetCompositeSvsPrime          
#
task Impl {
    input {
        File bnd_bcf
        Int max_distance

        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(bnd_bcf,"GB")) )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_MB=$(( ~{ram_size_gb}*1024 - 500 ))

        # Selecting inter-chromosomal BNDs
        date 1>&2
        bcftools view --threads ${N_THREADS} --output-type v ~{bnd_bcf} | awk 'BEGIN {OFS="\t"} /^#/{print;next} {alt=$5; sub(/^[^\[\]]*[\[\]]/,"",alt); sub(/:.*$/,"",alt); if (alt!=$1) print}' | bcftools view --threads ${N_THREADS} --output-type z --output inter.vcf.gz
        date 1>&2

        # Symmetrizing records
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M BndSymmetrize inter.vcf.gz 0 | bcftools sort --output-type z --output inter_symmetrized.vcf.gz

        # Computing clusters
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M GetCompositeSvsPrime inter_symmetrized.vcf.gz 1 ~{max_distance} 2 1 out.bed
        ${TIME_COMMAND} sort -k4,4nr out.bed > out_sorted.bed
        mv out_sorted.bed out.bed
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
