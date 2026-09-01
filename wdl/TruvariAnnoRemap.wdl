version 1.0

#
workflow TruvariAnnoRemap {
    input {
        File chunk_vcf_gz
        String remote_outdir

        File ref_fa
        File ref_fai

        Float cov_threshold = 0.8

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong_remap:latest"
    }

    call Impl {
        input:
            chunk_vcf_gz = chunk_vcf_gz,
            remote_outdir = remote_outdir,

            ref_fa = ref_fa,
            ref_fai = ref_fai,

            cov_threshold = cov_threshold,

            docker_image = docker_image
    }
}


# TOOL                                                CPU     RAM     TIME
# truvari anno remap                                  300%    20G     15m
#
# Remark: the reference index should be computed only once and provided as 
# input.
# Remark: RAM can reach >64GB in some cases.
#
task Impl {
    input {
        File chunk_vcf_gz
        String remote_outdir

        File ref_fa
        File ref_fai

        Float cov_threshold

        String docker_image
        Int n_cpu = 4
        Int mem_gb = 32
        Int preemptible_number = 2
    }

    Int disk_size_gb = 10*( ceil(size(chunk_vcf_gz,"GB")) )

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        

        # Initializing REF
        mkdir ./ref_files
        mv ~{ref_fa} ./ref_files/
        mv ~{ref_fai} ./ref_files/
        REF_FA_BASENAME=$(basename ~{ref_fa})

        # Mapping
        FILE_NAME=$(basename ~{chunk_vcf_gz} .vcf.gz)
        ${TIME_COMMAND} truvari anno remap --threads ${N_THREADS} --aligner minimap2 --min-length 1 --max-length 3000000000 --cov-threshold ~{cov_threshold} -r ./ref_files/${REF_FA_BASENAME} ~{chunk_vcf_gz} -o ${FILE_NAME}_remap.vcf.gz
        gcloud storage mv ${FILE_NAME}_remap.vcf.gz ~{remote_outdir}/
    >>>
    
    output {
    }

    runtime {
        cpu: n_cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " +  disk_size_gb + " HDD"
        preemptible: preemptible_number
        docker: docker_image
    }
}