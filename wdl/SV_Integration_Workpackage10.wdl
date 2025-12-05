version 1.0


# Merging by ID the re-genotyped personalized BCFs of all samples, in a given 
# chromosome chunk.
#
workflow SV_Integration_Workpackage10 {
    input {
        Int chunk_id
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
            sample_ids = sample_ids,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, one 10 MB chunk of chr6:
#
# TOOL               CPU     RAM     TIME
# ------>bcftools merge     120%    23G     30m
#
# Peak disk usage: ??????
#
task Impl {
    input {
        Int chunk_id
        File sample_ids
        
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 4
        Int ram_size_gb = 32
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        # Localizing all the samples for the given chunk
        N_SAMPLES=$(cat ~{sample_ids} | wc -l)
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir}/'*'_chunk_~{chunk_id}.'bcf*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading some files of chunk ~{chunk_id}. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        N_DOWNLOADED_SAMPLES=$(ls *_chunk_~{chunk_id}.bcf | wc -l)
        if [ ${N_DOWNLOADED_SAMPLES} -ne ${N_SAMPLES} ]; then
            echo "Error: the number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is different from the number of samples specified (${N_SAMPLES})."
            exit 1
        fi
        df -h
        
        # Merging records by ID, since the records in every BCF originate from
        # the same cohort BCF, which had distinct IDs.
        while read SAMPLE_ID; do
            echo ${SAMPLE_ID}_chunk_~{chunk_id}.bcf >> list.txt
        done < ~{sample_ids}
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge id --output-type b --file-list list.txt > chunk_~{chunk_id}.bcf
        bcftools index chunk_~{chunk_id}.bcf        
        ls -laht
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp chunk_~{chunk_id}.'bcf*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading merged BCF. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 3
    }
}
