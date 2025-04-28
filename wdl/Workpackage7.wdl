version 1.0


# 
#
workflow Workpackage7 {
    input {
        String chromosome_id
        Int chunk_id
        Boolean use_bed
        
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 2
        Int ram_size_gb = 16
        Int disk_size_gb = 100
    }
    parameter_meta {
        remote_indir: "Contains chunks of a bcftools merge VCF that need to be collapsed with truvari."
    }
    
    call Workpackage7Impl {
        input:
            chromosome_id = chromosome_id,
            chunk_id = chunk_id,
            use_bed = use_bed,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}

 
#
task Workpackage7Impl {
    input {
        String chromosome_id
        Int chunk_id
        Boolean use_bed
        
        String remote_indir
        String remote_outdir
        
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
        
        # Localizing the truvari collapse chunk
        while : ; do
            TEST=$(gsutil -m cp ~{chromosome_id}_chunk_~{chunk_id}.vcf.'gz*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading ~{chromosome_id} chunk ~{chunk_id}. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        if ~{use_bed} ; then
            while : ; do
                TEST=$(gsutil -m cp ~{chromosome_id}_included.bed . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading ~{chromosome_id} included BED. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            BED_ARGS="--bed ~{chromosome_id}_included.bed"
        else 
            BED_ARGS=" "
        fi
        
        # Collapsing
        source activate truvari5
        ${TIME_COMMAND} truvari collapse --input ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz --sizemin 0 --sizemax 1000000 --keep common ${BED_ARGS} --gt all --output tmp.vcf
        ${TIME_COMMAND} bcftools sort --max-mem $(( ~{ram_size_gb} - 2 ))G --output-type z tmp.vcf > ~{chromosome_id}_chunk_~{chunk_id}_truvari.vcf.gz
        tabix -f ~{chromosome_id}_chunk_~{chunk_id}_truvari.vcf.gz
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{chromosome_id}_chunk_~{chunk_id}_truvari.vcf.gz ~{remote_outdir}/~{chromosome_id}_chunk_~{chunk_id}.vcf.gz && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading truvari the collapse output. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{chromosome_id}_chunk_~{chunk_id}_truvari.vcf.gz.tbi ~{remote_outdir}/~{chromosome_id}_chunk_~{chunk_id}.vcf.gz.tbi && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading truvari the collapse output. Trying again..."
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
