version 1.0


# Computes `truvari collapse` chunks from all the `bcftools merge` chunks of a
# chromosome.
#
workflow SV_Integration_Workpackage6 {
    input {
        String chromosome_id
        String bcftools_chunks
        
        Int truvari_chunk_min_records = 2000
        Int truvari_collapse_refdist = 500
        
        String remote_indir
        String remote_outdir
    }
    parameter_meta {
        bcftools_chunks: "Comma-separated and sorted integers. Chunks are assumed to be sorted by POS."
        truvari_chunk_min_records: "Min number of records per output chunk. It has to be set to a small number (e.g. 2000) for truvari divide not to use too much RAM."
        truvari_collapse_refdist: "The actual collapse downstream will run `truvari collapse --refdist X`, where X is this value."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            chromosome_id = chromosome_id,
            bcftools_chunks = bcftools_chunks,
            truvari_chunk_min_records = truvari_chunk_min_records,
            truvari_collapse_refdist = truvari_collapse_refdist,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999:
#
# TOOL                          CPU     RAM     TIME
# bcftools concat               12%     20M     1h
# truvari divide --min 1769     100%    12G     6h         // Outputs 562 chunks
#
# Output of truvari divide (each chunk is ~5 MB):
# count      562.000000
# mean      3148.354093
# std       2623.478656
# min         68.000000
# 25%       1849.000000
# 50%       2070.500000
# 75%       3062.250000
# max      19218.000000
#
# Remark: truvari divide --min 88468 goes out of RAM with 256GB.
# 
task Impl {
    input {
        String chromosome_id
        String bcftools_chunks
        
        Int truvari_chunk_min_records
        Int truvari_collapse_refdist
        
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
        
        
        # Localizing all the bcftools merge chunks of the chromosome
        rm -f uri_list.txt file_list.txt
        for CHUNK in $(echo ~{bcftools_chunks} | tr ',' ' '); do
            echo ~{remote_indir}/chunk_${CHUNK}.bcf >> uri_list.txt
            echo ~{remote_indir}/chunk_${CHUNK}.bcf.csi >> uri_list.txt
            echo chunk_${CHUNK}.bcf >> file_list.txt
        done
        ${TIME_COMMAND} xargs --arg-file=uri_list.txt --max-lines=1 --max-procs=${N_THREADS} -I {} gcloud storage cp {} .
        df -h
        rm -f uri_list.txt
        
        # Concatenating all the bcftools merge chunks to build a whole-
        # chromosome VCF. This is necessary, since a truvari collapse chunks may
        # straddle different bcftools merge chunks.
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --file-list file_list.txt --output-type z1 > ~{chromosome_id}.vcf.gz
        bcftools index --threads ${N_THREADS} -f ~{chromosome_id}.vcf.gz
        df -h
        rm -f chunk_*.bcf* file_list.txt
        
        # Chunking the chromosome for truvari collapse.
        # Remark: `truvari divide` chunks the VCF directly (i.e. it does not
        # just compute chunk boundaries), so we need to feed it the full
        # bcftools merge VCF with all genotypes.
        N_RECORDS=$(bcftools index --nrecords ~{chromosome_id}.vcf.gz.tbi)
        ${TIME_COMMAND} truvari divide --threads ${N_THREADS} --min ~{truvari_chunk_min_records} --buffer ~{truvari_collapse_refdist} ~{chromosome_id}.vcf.gz ./truvari_chunks/
        i="0"
        N_RECORDS_CHUNKED="0"
        for FILE in $(ls ./truvari_chunks/*.vcf.gz | sort -V); do
            mv ${FILE} chunk_${i}.vcf.gz
            mv ${FILE}.tbi chunk_${i}.vcf.gz.tbi
            N=$( bcftools index --nrecords chunk_${i}.vcf.gz.tbi )
            N_RECORDS_CHUNKED=$(( ${N_RECORDS_CHUNKED} + ${N} ))
            i=$(( ${i} + 1 ))
        done
        if [ ${N_RECORDS_CHUNKED} -ne ${N_RECORDS} ]; then
            echo "ERROR: The truvari collapse chunks contain ${N_RECORDS_CHUNKED} total records, but the chromosome VCF contains ${N_RECORDS} records."
            exit 1
        fi
        
        # Uploading
        ls chunk_*.vcf.gz* > file_list.txt
        xargs --arg-file=file_list.txt --max-lines=1 --max-procs=${N_THREADS} -I {} gcloud storage cp {} ~{remote_outdir}/~{chromosome_id}/
        rm -f file_list.txt
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
