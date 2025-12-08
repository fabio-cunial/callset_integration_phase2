version 1.0


# Combines all the bcftools merge chunks of a re-genotyped chromosome.
#
workflow SV_Integration_Workpackage11 {
    input {
        String chromosome_id
        String bcftools_chunks
        
        String remote_indir
        String remote_outdir
    }
    parameter_meta {
        bcftools_chunks: "Comma-separated and sorted integers. Chunk order is assumed to reflect POS order, and every chunk is assumed to be sorted by POS."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            chromosome_id = chromosome_id,
            bcftools_chunks = bcftools_chunks,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999:
#
# TOOL                          CPU     RAM     TIME
# bcftools concat               12%??     20M??     1h??
# 
task Impl {
    input {
        String chromosome_id
        String bcftools_chunks
        
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 2
        Int ram_size_gb = 4
        Int disk_size_gb = 100
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
            echo ~{remote_indir}/chunk_${CHUNK}.bcf >> uri_list.txt
            echo ~{remote_indir}/chunk_${CHUNK}.bcf.csi >> uri_list.txt
            echo chunk_${CHUNK}.bcf >> file_list.txt
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
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list file_list.txt --output-type b > ~{chromosome_id}.bcf
        bcftools index ~{chromosome_id}.bcf
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{chromosome_id}.'bcf*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading the result of the merge. Trying again..."
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
