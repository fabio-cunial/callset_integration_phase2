version 1.0


# For a given chromosome chunk, the procedure performs a CAL_SENS filter and a
# bcftools merge.
#
workflow Workpackage5 {
    input {
        Int chunk_id
        String filter_string
        File sample_ids
        
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 4
        Int ram_size_gb = 128
        Int disk_size_gb = 100
    }
    parameter_meta {
        sample_ids: "Speficies the order of the samples used by bcftools merge."
    }
    
    call Workpackage5Impl {
        input:
            chunk_id = chunk_id,
            filter_string = filter_string,
            sample_ids = sample_ids,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}


# Performance on 10'070 samples, one chunk, 15x, GRCh38:
#
# CAL_SENS  TOOL                      CPU     RAM     TIME
# <=0.7     bcftools merge            150%    26G     1h30m
# <=0.7     bcftools norm             150%    5G      20m
# <=0.9     bcftools merge            
# <=0.9     bcftools norm             
#
task Workpackage5Impl {
    input {
        Int chunk_id
        String filter_string
        File sample_ids
        
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 4
        Int ram_size_gb = 128
        Int disk_size_gb = 100
    }
    parameter_meta {
    }
    
    String docker_dir = "/root"
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        # Localizing all the samples for the given chunk
        N_SAMPLES=$(cat ~{sample_ids} | wc -l)
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir}/'*'_chunk_~{chunk_id}.vcf.'gz*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading some files of chunk ~{chunk_id}. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        N_DOWNLOADED_SAMPLES=$(ls *_chunk_~{chunk_id}.vcf.gz | wc -l)
        if [ ${N_DOWNLOADED_SAMPLES} -ne ${N_SAMPLES} ]; then
            echo "Error: the number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is different from the number of samples specified (${N_SAMPLES})."
            exit 1
        fi
        
        # Filtering
        FILTER_STRING="~{filter_string}"
        if [ ${FILTER_STRING} != none ]; then
            INCLUDE_STR="--include ${FILTER_STRING}"
            while read SAMPLE_ID; do
                bcftools filter --threads ${N_THREADS} ${INCLUDE_STR} --output-type z ${SAMPLE_ID}_chunk_~{chunk_id}.vcf.gz > ${SAMPLE_ID}_filtered.vcf.gz
                tabix -f ${SAMPLE_ID}_filtered.vcf.gz
                echo ${SAMPLE_ID}_filtered.vcf.gz >> list.txt
                rm -f ${SAMPLE_ID}_chunk_~{chunk_id}.vcf.gz*
            done < ~{sample_ids}
        else
            while read SAMPLE_ID; do
                echo ${SAMPLE_ID}_chunk_~{chunk_id}.vcf.gz >> list.txt
            done < ~{sample_ids}
        fi
        df -h
        
        # Bcftools merge
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type z > ~{chunk_id}_merged.vcf.gz
        tabix -f ~{chunk_id}_merged.vcf.gz
        ls -laht ~{chunk_id}_merged.vcf.gz
        df -h
        while read FILE; do
            rm -f ${FILE}
        done < list.txt
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type z ~{chunk_id}_merged.vcf.gz > ~{chunk_id}_normed.vcf.gz
        tabix -f ~{chunk_id}_normed.vcf.gz
        ls -laht ~{chunk_id}_merged.vcf.gz
        
        # Uploading
        gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ~{chunk_id}_normed.vcf.gz ~{remote_outdir}/chunk_~{chunk_id}.vcf.gz
        gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ~{chunk_id}_normed.vcf.gz.tbi ~{remote_outdir}/chunk_~{chunk_id}.vcf.gz.tbi
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
