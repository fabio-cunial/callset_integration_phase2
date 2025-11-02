version 1.0


# Concatenates chunks of the truvari collapse output. Prepares the concatenated
# file for kanpig, by dropping all genotypes and by assigning a unique integer
# ID to every call (the original ID is copied to the INFO field): this is
# necessary, otherwise kanpig may complain about duplicated IDs, which may arise
# naturally from the previous steps of the pipeline.
#
workflow Workpackage8 {
    input {
        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"]
        String remote_indir
        String remote_outdir
    }
    parameter_meta {
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
        chromosomes: "The order of the chromosomes becomes their order in the output VCF."
    }
    
    call Impl {
        input:
            chromosomes = chromosomes,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


# Performance on 10'070 samples, 15x, GRCh38, SSD:
#
# CAL_SENS  TOOL                                CPU     RAM     TIME
# <=0.7     bcftools concat                     50%     300M    2m
# <=0.7     tabix                               100%    90M     20m
# <=0.7     bgzip                               360%    10M     4m
#
# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, SSD:
#
# TOOL                CPU     RAM     TIME
# bcftools concat     70%     90M     10s
# tabix               100%    10M     3m
#
task Impl {
    input {
        Array[String] chromosomes
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 4
        Int ram_size_gb = 4
        Int disk_size_gb = 200
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
        
        
        # Localizing all the chunks of all the chromosomes
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
        cat chunk_list.txt
        df -h
        
        # Concatenating all the chunks
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list chunk_list.txt --output-type z > truvari_collapsed.vcf.gz
        ${TIME_COMMAND} tabix -f truvari_collapsed.vcf.gz
        rm -rf *_chunk_*
        df -h
        
        # Preparing the inter-sample VCF for kanpig
        bcftools view --header-only truvari_collapsed.vcf.gz > header.txt
        N_ROWS=$(wc -l < header.txt)
        (  head -n $(( ${N_ROWS} - 1 )) header.txt ; \
           echo '##INFO=<ID=ORIGINAL_ID,Number=1,Type=String,Description="Original ID from truvari collapse">' ; \
           echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" ; \
           bcftools view --no-header truvari_collapsed.vcf.gz | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { gsub(/;/,"_",$3); printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s;ORIGINAL_ID=%s\tGT\t0/1\n",$1,$2,++i,$4,$5,$6,$7,$8,$3); }' \
        ) | bgzip --compress-level 1 > truvari_collapsed_for_kanpig.vcf.gz
        ${TIME_COMMAND} tabix -f truvari_collapsed_for_kanpig.vcf.gz
        df -h
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp 'truvari_collapsed.vcf.gz*' 'truvari_collapsed_for_kanpig.vcf.gz*' ~{remote_outdir}/ && echo 0 || echo 1)
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
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
