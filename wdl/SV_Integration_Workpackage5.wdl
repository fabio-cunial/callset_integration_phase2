version 1.0


# For a given chromosome chunk, the procedure performs a CAL_SENS filter and a
# bcftools merge.
#
workflow SV_Integration_Workpackage5 {
    input {
        Int chunk_id
        String filter_string = "FORMAT/CALIBRATION_SENSITIVITY<=0.999"
        File sample_ids
        
        String remote_indir
        String remote_outdir
    }
    parameter_meta {
        sample_ids: "Speficies the order of the samples used by bcftools merge."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            chunk_id = chunk_id,
            filter_string = filter_string,
            sample_ids = sample_ids,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, one 10 MB chunk of chr6,
# CAL_SENS<=0.999:
#
# TOOL               CPU     RAM     TIME
# bcftools merge     120%    23G     30m
# bcftools norm      150%    4G      5m
#
task Impl {
    input {
        Int chunk_id
        String filter_string
        File sample_ids
        
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 4
        Int ram_size_gb = 128
        Int disk_size_gb = 50
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
        
# -----------> FILTERING SHOULD HAVE BEEN DONE BY THE PREVIOUS WORKPACKAGE!!!!!        

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
        
        # Merging
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type z > ~{chunk_id}_merged.vcf.gz
        tabix -f ~{chunk_id}_merged.vcf.gz
        rm -f *_filtered.vcf.gz
        
        # Making sure no multiallelic record is passed downstream
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics - --output-type z ~{chunk_id}_merged.vcf.gz > ~{chunk_id}_normed.vcf.gz
        tabix -f ~{chunk_id}_normed.vcf.gz
        ls -laht
        
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
