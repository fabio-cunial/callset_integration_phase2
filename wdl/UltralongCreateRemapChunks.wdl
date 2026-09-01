version 1.0


# 
#
workflow UltralongCreateRemapChunks {
    input {
        File input_bcf
        Int n_chunks
        String remote_outdir

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            input_bcf = input_bcf,
            n_chunks = n_chunks,
            remote_outdir = remote_outdir,

            docker_image = docker_image
    }
    
    output {
    }
}


# COMMAND                               CPU     RAM     TIME
# bcftools view --drop-genotypes        170%    50M     7s
# chunk creation                                        8s     
#
task Impl {
    input {
        File input_bcf
        Int n_chunks
        String remote_outdir

        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(input_bcf,"GB")) )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_MB=$(( ~{ram_size_gb}*1024 - 500 ))


        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --drop-genotypes --include 'SVTYPE="INS"' --output-type b1 ~{input_bcf} --output ins.bcf
        rm -f ~{input_bcf}
        ${TIME_COMMAND} bcftools index ins.bcf
        N_RECORDS=$(bcftools index --nrecords ins.bcf.csi)
        N_RECORDS_PER_CHUNK=$(( (${N_RECORDS} + ~{n_chunks} - 1) / ~{n_chunks} ))
        bcftools view --header-only ins.bcf > header.vcf
        date 1>&2
        bcftools view --threads ${N_THREADS} --no-header ins.bcf | split -d -a 4 -l ${N_RECORDS_PER_CHUNK} - chunk_
        date 1>&2
        for CHUNK_FILE in chunk_*; do
            cat header.vcf "${CHUNK_FILE}" | bgzip --threads ${N_THREADS} --compress-level 1 > ${CHUNK_FILE}.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${CHUNK_FILE}.vcf.gz
        done
        date 1>&2
        gcloud storage cp 'chunk_*.vcf.gz*' ~{remote_outdir}
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
