version 1.0


# 
#
workflow Workpackage7 {
    input {
        String chromosome_id
        Int chunk_id
        Boolean use_bed
        String truvari_flags = "--sizemin 0 --sizemax 1000000 --gt off --keep maxqual"
        Boolean drop_gts = false
        
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 100
    }
    parameter_meta {
        remote_indir: "Contains chunks of a bcftools merge VCF that need to be collapsed with truvari."
        drop_gts: "Remove all the sample GT fields from the VCF before running truvari collapse. This was suggested e.g. in https://github.com/ACEnglish/truvari/issues/220#issuecomment-2830920205"
        truvari_flags: "`--gt all` is very slow on 10k samples. See e.g. https://github.com/ACEnglish/truvari/issues/220#issuecomment-2830920205"
    }
    
    call Workpackage7Impl {
        input:
            chromosome_id = chromosome_id,
            chunk_id = chunk_id,
            use_bed = use_bed,
            truvari_flags = truvari_flags,
            drop_gts = drop_gts,
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
        String truvari_flags
        Boolean drop_gts
        
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
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.21/plugins"
        
        # Localizing the truvari collapse chunk
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir}/~{chromosome_id}_chunk_~{chunk_id}.vcf.'gz*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading ~{chromosome_id} chunk ~{chunk_id}. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        if ~{use_bed} ; then
            while : ; do
                TEST=$(gsutil -m cp ~{remote_indir}/~{chromosome_id}_included.bed . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading ~{chromosome_id} included BED. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            BED_FLAGS="--bed ~{chromosome_id}_included.bed"
        else 
            BED_FLAGS=" "
        fi
        N_RECORDS=$( bcftools index --nrecords ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz.tbi )
        
        # Writing AF in the QUAL field
        ${TIME_COMMAND} bcftools +fill-tags ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz -Oz -o tmp.vcf.gz -- -t AF
        tabix -f tmp.vcf.gz
        bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' tmp.vcf.gz | bgzip -c > annotations.tsv.gz
        tabix -s1 -b2 -e2 annotations.tsv.gz
        rm -f tmp.vcf.gz*
        ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations annotations.tsv.gz --columns CHROM,POS,ID,REF,ALT,QUAL ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz --output-type z > tmp.vcf.gz
        rm -f ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz*
        mv tmp.vcf.gz ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz
        tabix -f ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz*
        
        # Removing GTs
        if ~{drop_gts}; then
            ${TIME_COMMAND} bcftools view --drop-genotypes ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz --output-type z > tmp.vcf.gz
            rm -f ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz*
            mv tmp.vcf.gz ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz
            tabix -f ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz
        fi
        
        # Collapsing
        source activate truvari5
        ${TIME_COMMAND} truvari collapse --input ~{chromosome_id}_chunk_~{chunk_id}.vcf.gz ~{truvari_flags} ${BED_FLAGS} --output tmp.vcf
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
