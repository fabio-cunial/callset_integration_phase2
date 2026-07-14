version 1.0


# 
#
workflow SameChromPos {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            remote_indir = remote_indir,

            docker_image = docker_image
    }
    
    output {
        File stats_lt_50 = Impl.stats_lt_50
        File stats_ge_50 = Impl.stats_ge_50
    }
}


task Impl {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        
        String docker_image
        Int n_cpu = 1
        Int ram_size_gb = 8
        Int disk_size_gb = 20
        Int preemptible_number = 0
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
        
        touch stats_lt_50.csv stats_ge_50.csv
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}_kanpig.vcf.'gz*' .
            bcftools filter --include 'SVLEN<50' --output-type v ${SAMPLE_ID}_kanpig.vcf.gz --output ${SAMPLE_ID}_lt_50.vcf
            bcftools filter --include 'SVLEN>=50' --output-type v ${SAMPLE_ID}_kanpig.vcf.gz --output ${SAMPLE_ID}_ge_50.vcf
            java -cp ~{docker_dir} SameChromPos ${SAMPLE_ID}_lt_50.vcf >> stats_lt_50.csv
            java -cp ~{docker_dir} SameChromPos ${SAMPLE_ID}_ge_50.vcf >> stats_ge_50.csv

            rm -f ${SAMPLE_ID}_*
        done 3< chunk.csv
    >>>
    
    output {
        File stats_lt_50 = "stats_lt_50.csv"
        File stats_ge_50 = "stats_ge_50.csv"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
