version 1.0


# For a given chromosome, computes truvari collapse chunks from bcftools merge
# chunks.
#
workflow Workpackage6 {
    input {
        String chromosome_id
        String bcftools_chunks
        
        String remote_indir
        String remote_outdir
        
        File reference_fai
        File density_counter_py
        Int max_records_per_chunk = 10000
        Int truvari_chunk_length = 1000000
        Int slack_bp = 2
        
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 100
    }
    parameter_meta {
        bcftools_chunks: "Comma-separated. Chunks must be sorted by POS."
        max_records_per_chunk: "Chunks with more records than this are excluded from truvari collapse."
        truvari_chunk_length: "Min number of bps that a truvari collapse chunk should have. On 10k samples, setting this to 1M gives ~3k chunks for CAL_SENS<=0.7."
        slack_bp: "Extend a chunk [start,end) emitted by `density_counter_py` to [start,end+slack_bp), just for boundary safety."
    }
    
    call Workpackage6Impl {
        input:
            chromosome_id = chromosome_id,
            bcftools_chunks = bcftools_chunks,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            reference_fai = reference_fai,
            density_counter_py = density_counter_py,
            max_records_per_chunk = max_records_per_chunk,
            truvari_chunk_length = truvari_chunk_length,
            slack_bp = slack_bp,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}

   
#
task Workpackage6Impl {
    input {
        String chromosome_id
        String bcftools_chunks
        
        String remote_indir
        String remote_outdir
        
        File reference_fai
        File density_counter_py
        Int max_records_per_chunk
        Int truvari_chunk_length
        Int slack_bp
        
        Int n_cpu
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
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
        
        
        # Localizing all the bcftools merge chunks of the given chromosome
        for CHUNK in $(echo ~{bcftools_chunks} | tr ',' ' '); do
            while : ; do
                TEST=$(gsutil -m cp ~{remote_indir}/chunk_${CHUNK}.vcf.'gz*' . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading chunk ${CHUNK}. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    echo chunk_${CHUNK}.vcf.gz >> list.txt
                    break
                fi
            done
        done
        
        # Concatenating all the bcftools merge chunks
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list list.txt --output-type z > ~{chromosome_id}.vcf.gz
        tabix -f ~{chromosome_id}.vcf.gz
        rm -f chunk_*.vcf.gz*
        
        # Creating truvari collapse chunks
        source activate truvari4
        ${TIME_COMMAND} python ~{density_counter_py} ~{chromosome_id}.vcf.gz > ~{chromosome_id}_tmp.bed
        conda deactivate
        ${TIME_COMMAND} bedtools sort -i ~{chromosome_id}_tmp.bed > ~{chromosome_id}_chunks.bed
        rm -f ~{chromosome_id}_tmp.bed
        java -cp ~{docker_dir} ChunkHistogram ~{chromosome_id}_chunks.bed
        awk '$4 >= ~{max_records_per_chunk}' ~{chromosome_id}_chunks.bed > ~{chromosome_id}_excluded.bed
        cat ~{chromosome_id}_excluded.bed
        ${TIME_COMMAND} bedtools complement -i ~{chromosome_id}_excluded.bed -g ~{reference_fai} > ~{chromosome_id}_tmp.bed
        ${TIME_COMMAND} bedtools sort -i ~{chromosome_id}_tmp.bed > ~{chromosome_id}_included.bed
        rm -f ~{chromosome_id}_tmp.bed
        ${TIME_COMMAND} java -cp ~{docker_dir} SplitForTruvariCollapse ~{chromosome_id} ~{chromosome_id}_chunks.bed ~{reference_fai} ~{truvari_chunk_length} ~{slack_bp} > ~{chromosome_id}_truvari_chunks.csv
        
        # Partitioning the chromosome VCF
        i="0"
        N_RECORDS_BEFORE=$( bcftools index --nrecords ~{chromosome_id}.vcf.gz.tbi )
        N_RECORDS_AFTER="0"
        while read INTERVAL; do
            echo ${INTERVAL} | tr ',' '\t' > tmp.bed
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --regions-file tmp.bed --regions-overlap pos --output-type z ~{chromosome_id}.vcf.gz > ~{chromosome_id}_chunk_${i}.vcf.gz
            tabix -f ~{chromosome_id}_chunk_${i}.vcf.gz
            N=$( bcftools index --nrecords ~{chromosome_id}_chunk_${i}.vcf.gz.tbi )
            N_RECORDS_AFTER=$(( ${N_RECORDS_AFTER} + ${N} ))
            i=$(( ${i} + 1 ))
        done < ~{chromosome_id}_truvari_chunks.csv
        if [ ${N_RECORDS_AFTER} -ne ${N_RECORDS_BEFORE} ]; then
            echo "ERROR: the truvari collapse chunks contain ${N_RECORDS_AFTER} total records, but the chromosome VCF contains ${N_RECORDS_BEFORE} records."
            exit 1
        fi
        
        # Uploading the truvari collapse chunks and the BED files
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{chromosome_id}_chunk_'*'.vcf.'gz*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading truvari collapse chunks. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{chromosome_id}_'*'.bed ~{chromosome_id}_'*'.csv ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading BED/CSV files. Trying again..."
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
