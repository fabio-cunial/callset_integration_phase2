version 1.0


# Computes truvari collapse chunks from all the bcftools merge chunks of a
# chromosome.
#
workflow SV_Integration_Workpackage6 {
    input {
        String chromosome_id
        String bcftools_chunks
        
        Int n_truvari_chunks
        Int truvari_refdist
        
        String remote_indir
        String remote_outdir
    }
    parameter_meta {
        bcftools_chunks: "Comma-separated and sorted integers. Chunks are assumed to be sorted by POS."
        n_truvari_chunks: "Desired number of chunks. The number of chunks in output might be different from this."
        truvari_refdist: "The actual collapse downstream will run `truvari collapse --refdist X`, where X is this value."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            chromosome_id = chromosome_id,
            bcftools_chunks = bcftools_chunks,
            n_truvari_chunks = n_truvari_chunks,
            truvari_refdist = truvari_refdist,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999:
#
# TOOL               CPU     RAM     TIME
# bcftools concat    12%     20M     1h
# truvari divide     
# 
task Impl {
    input {
        String chromosome_id
        String bcftools_chunks
        
        Int n_truvari_chunks
        Int truvari_refdist
        
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 4
        Int ram_size_gb = 16
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
        
        
        # Localizing all the bcftools merge chunks of the chromosome
        for CHUNK in $(echo ~{bcftools_chunks} | tr ',' ' '); do
            echo ~{remote_indir}/chunk_${CHUNK}.vcf.gz >> uri_list.txt
            echo ~{remote_indir}/chunk_${CHUNK}.vcf.gz.tbi >> uri_list.txt
            echo chunk_${CHUNK}.vcf.gz >> file_list.txt
        done
        while : ; do
            TEST=$(cat uri_list.txt | gsutil -m cp -I . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading bcftools merge chunks. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        # Concatenating all the bcftools merge chunks
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list file_list.txt --output-type z > ~{chromosome_id}.vcf.gz
        tabix -f ~{chromosome_id}.vcf.gz
        rm -f chunk_*.vcf.gz*
        
        # Creating truvari collapse chunks
        N_RECORDS=$(bcftools index --nrecords ~{chromosome_id}.vcf.gz.tbi)
        MIN_RECORDS_PER_CHUNK=$(( ( ${N_RECORDS} / ~{n_truvari_chunks} ) / 10 ))
        ${TIME_COMMAND} truvari divide --threads ${N_THREADS} --min ${MIN_RECORDS_PER_CHUNK} --buffer ~{truvari_refdist} ~{chromosome_id}.vcf.gz ./truvari_chunks/
        i="0"
        N_RECORDS_CHUNKED="0"
        for FILE in $(ls ./truvari_chunks/*.vcf.gz | sort -V); do
            mv ${FILE} ~{chromosome_id}_chunk_${i}.vcf.gz
            mv ${FILE}.tbi ~{chromosome_id}_chunk_${i}.vcf.gz.tbi
            N=$( bcftools index --nrecords ~{chromosome_id}_chunk_${i}.vcf.gz.tbi )
            N_RECORDS_CHUNKED=$(( ${N_RECORDS_CHUNKED} + ${N} ))
            i=$(( ${i} + 1 ))
        done
        if [ ${N_RECORDS_CHUNKED} -ne ${N_RECORDS} ]; then
            echo "ERROR: the truvari collapse chunks contain ${N_RECORDS_CHUNKED} total records, but the chromosome VCF contains ${N_RECORDS} records."
            exit 1
        fi
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{chromosome_id}_chunk_'*'.vcf.'gz*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading truvari collapse chunks. Trying again..."
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
        preemptible: 0
    }
}
