version 1.0


# Concatenates chunks of the truvari collapse output. Prepares the concatenated
# file for kanpig, by dropping all genotypes and by assigning a unique integer
# ID to every call (moving the original ID to the INFO field): this is
# necessary, otherwise kanpig complains about duplicated IDs, which arise
# naturally from the previous steps of the pipeline.
#
workflow Workpackage8 {
    input {
        String remote_indir
        String remote_outdir
        
        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 1000
    }
    parameter_meta {
        chromosomes: "The order of the chromosomes becomes their order in the output VCF."
    }
    
    call Workpackage8Impl {
        input:
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            chromosomes = chromosomes,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}


#
task Workpackage8Impl {
    input {
        String remote_indir
        String remote_outdir
        
        Array[String] chromosomes
        
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
        
        
        # Localizing all the truvari collapsed chunks of all chromosomes
        CHROMOSOMES=~{sep=',' chromosomes}
        echo ${CHROMOSOMES} | tr ',' '\n' > chr_list.txt
        while read CHROMOSOME; do
            while : ; do
                TEST=$(gsutil -m cp ~{remote_indir}/${CHROMOSOME}_chunk_'*.vcf.gz*' . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading the chunks of ${CHROMOSOME}. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        done < chr_list.txt
        while read CHROMOSOME; do
            ls ${CHROMOSOME}_chunk_*.vcf.gz | sort -V >> chunk_list.txt
        done < chr_list.txt
        rm -f chr_list.txt
        cat chunk_list.txt
        df -h
        
        # Concatenating all the truvari collapsed chunks
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list chunk_list.txt --output-type z > truvari_collapsed.vcf.gz
        tabix -f truvari_collapsed.vcf.gz
        df -h
        
        # Preparing the inter-sample VCF for kanpig
        ${TIME_COMMAND} bcftools view --drop-genotypes truvari_collapsed.vcf.gz --output-type z > truvari_collapsed_nogts.vcf.gz
        tabix -f truvari_collapsed_nogts.vcf.gz
        ${TIME_COMMAND} bcftools view --header-only truvari_collapsed_nogts.vcf.gz > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > truvari_collapsed_for_kanpig.vcf
        echo '##INFO=<ID=ORIGINAL_ID,Number=1,Type=String,Description="Original ID from truvari collapse">' >> truvari_collapsed_for_kanpig.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> truvari_collapsed_for_kanpig.vcf
        bcftools view --no-header truvari_collapsed_no_gts.vcf.gz | awk 'BEGIN { i=0; } { gsub(/;/,"_",$3); printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s;ORIGINAL_ID=%s\n",$1,$2,++i,$4,$5,$6,$7,$8,$3); }' >> truvari_collapsed_for_kanpig.vcf
        rm -f truvari_collapsed_no_gts.vcf.gz*
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} truvari_collapsed_for_kanpig.vcf
        tabix -f truvari_collapsed_for_kanpig.vcf.gz
        df -h
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp 'truvari_collapsed.vcf.gz*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading concatenated VCFs. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp 'truvari_collapsed_for_kanpig.vcf.gz*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading concatenated VCFs. Trying again..."
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
