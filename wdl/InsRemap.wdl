version 1.0

#
workflow InsRemap {
    input {
        String chunk_id
        File chunk_csv

        String remote_indir
        String remote_outdir

        File ref_fa
        File ref_fai

        Int max_length = 1000000
        Float cov_threshold = 0.8

        String docker_image = "quay.io/ymostovoy/lr-remap:latest"
    }

    call Impl {
        input:
            chunk_id = chunk_id,
            chunk_csv = chunk_csv,

            remote_indir = remote_indir,
            remote_outdir = remote_outdir,

            ref_fa = ref_fa,
            ref_fai = ref_fai,

            max_length = max_length,
            cov_threshold = cov_threshold,

            docker_image = docker_image
    }
}


#
task Impl {
    input {
        String chunk_id
        File chunk_csv

        String remote_indir
        String remote_outdir

        File ref_fa
        File ref_fai

        Int max_length
        Float cov_threshold

        String docker_image
        Int n_cpu = 8
        Int mem_gb = 32
        Int disk_size_gb = 20
    }

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        

        # Initializing REF
        mkdir ./ref_files
        mv ~{ref_fa} ./ref_files/
        mv ~{ref_fai} ./ref_files/
        REF_FA_BASENAME=$(basename ~{ref_fa})

        # Iterating over samples
        touch ~{chunk_id}_matrix.csv
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            # Skipping the sample if it is not in the input
            TEST=$( gcloud storage ls ~{remote_indir}/${SAMPLE_ID}_ins.vcf.gz || echo "0" )
            if [ ${TEST} = "0" ]; then
                continue
            fi

            # Skipping the sample if it has already been processed
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            # Remapping
            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}_ins.vcf.'gz*' .
            time truvari anno remap --threads ${N_THREADS} --aligner minimap2 --min-length 1 --max-length ~{max_length} --cov-threshold ~{cov_threshold} -r ./ref_files/${REF_FA_BASENAME} ${SAMPLE_ID}_ins.vcf.gz -o ${SAMPLE_ID}_ins_remapped.vcf.gz
            time bcftools index --threads ${N_THREADS} --tbi ${SAMPLE_ID}_ins_remapped.vcf.gz
            bcftools query --format '%INFO/SUPP_SNIFFLES,%INFO/SUPP_PBSV,%INFO/SUPP_PAV,%INFO/SVLEN,%INFO/remap_classification,%INFO/remap_perc\n' ${SAMPLE_ID}_ins_remapped.vcf.gz >> ~{chunk_id}_matrix.csv
            tail -n 10 ~{chunk_id}_matrix.csv 1>&2
            
            # Next iteration
            gcloud storage mv ${SAMPLE_ID}_ins_remapped.vcf.'gz*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            rm -f ${SAMPLE_ID}*
            ls -laht 1>&2
        done 3< ~{chunk_csv}
    >>>
    
    output {
        File matrix_csv = "~{chunk_id}_matrix.csv"
    }

    runtime {
        cpu: n_cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " +  disk_size_gb + " HDD"
        preemptible: 0
        docker: docker_image
    }
}