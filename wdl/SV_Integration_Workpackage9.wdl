version 1.0


# Re-genotypes a personalized cohort VCF using the BAM of each sample, and
# partitions the result into ~100 pieces to run bcftools merge over all samples
# on parallel chunks.
#
workflow SV_Integration_Workpackage9 {
    input {
        File sv_integration_chunk_tsv
        File split_for_bcftools_merge_csv
        String remote_indir
        String remote_outdir
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        
        String kanpig_params_cohort = "--neighdist 500 --gpenalty 0.04 --hapsim 0.97"
    }
    parameter_meta {
        sv_integration_chunk_tsv: "A subset of the rows of table `sv_integration_hg38`, without the header."
        split_for_bcftools_merge_csv: "A partition that covers all chromosomes. Every line is a 0-based, half-open, consecutive chunk of a chromosome. Lines are assumed to be sorted."
        remote_indir: "Without final slash. Stores frequent and infrequent BCFs."
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            kanpig_params_cohort = kanpig_params_cohort
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, 8 logical cores, 8GB RAM:
#
# TOOL                      CPU     RAM     TIME
# bcftools view --samples                   3m
# bcftools concat           600%    34M     3s
# kanpig                    400%    500M    3m
# bcftools sort             100%    800M    40s
#
# Performance on 12'680 samples, 30x, GRCh38, 16 logical cores, 16GB RAM:
#
# TOOL                      CPU     RAM     TIME
# bcftools view --samples                   3m                
# bcftools concat           900%    50M     3s
# kanpig                    100%    1G      3m
# bcftools sort             100%    800M    20s
#
task Impl {
    input {
        File sv_integration_chunk_tsv
        File split_for_bcftools_merge_csv
        String remote_indir
        String remote_outdir
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        String kanpig_params_cohort
        
        Int n_cpu = 8
        Int ram_size_gb = 8
        Int disk_size_gb = 50
    }
    parameter_meta {
        disk_size_gb: "50GB is enough for most 15x samples, but sometimes it is not sufficient. 200GB is enough for all 30x samples."
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
        export RUST_BACKTRACE="full"
        INFINITY="1000000000"
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 3)
            ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 4)
            while : ; do
                TEST=$(gsutil -m cp ${ALIGNED_BAM} ./${SAMPLE_ID}_aligned.bam && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${ALIGNED_BAM}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gsutil -m cp ${ALIGNED_BAI} ./${SAMPLE_ID}_aligned.bam.bai && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${ALIGNED_BAI}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            touch ${SAMPLE_ID}_aligned.bam.bai
        }
        
        
        # Deletes all files related to a sample
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}_*
        }
        
        
        function Kanpig() {
            local SAMPLE_ID=$1
            local SEX=$2

            if [ ${SEX} == "M" ]; then
                PLOIDY_BED=$(echo ~{ploidy_bed_male})
            else
                PLOIDY_BED=$(echo ~{ploidy_bed_female})
            fi
            
            # Remark: we set --sizemin 10 instead of zero or one, just because
            # kanpig needs --sizemin >= --kmer . The purpose is still to
            # re-genotype every record in the input VCF.
            ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --sizemin 10 --sizemax ${INFINITY} ~{kanpig_params_cohort} --reference ~{reference_fa} --ploidy-bed ${PLOIDY_BED} --input ${SAMPLE_ID}_personalized.vcf.gz --reads ${SAMPLE_ID}_aligned.bam --out ${SAMPLE_ID}_tmp.vcf --sample ${SAMPLE_ID}
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type b ${SAMPLE_ID}_tmp.vcf > ${SAMPLE_ID}_kanpig.bcf
            rm -f ${SAMPLE_ID}_tmp.vcf
            bcftools index ${SAMPLE_ID}_kanpig.bcf
            
            # Basic statistics on the number of present records
            N_RECORDS=$(bcftools index --nrecords ${SAMPLE_ID}_kanpig.bcf.csi)
            N_PRESENT_RECORDS=$(bcftools view --no-header --include 'GT="alt"' ${SAMPLE_ID}_kanpig.bcf | wc -l)
            echo "${N_RECORDS},${N_PRESENT_RECORDS}" > ${SAMPLE_ID}_kanpig_nrecords.txt
        }
        
        
        # Splits the re-genotyped BCF into chunks for bcftools merge downstream.
        # 
        function ChunkKanpig() {
            local SAMPLE_ID=$1
            
            i="0"
            while read INTERVAL; do
                echo ${INTERVAL} | tr ',' '\t' > ${SAMPLE_ID}.bed
                ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --regions-file ${SAMPLE_ID}.bed --regions-overlap pos --output-type b ${SAMPLE_ID}_kanpig.bcf > ${SAMPLE_ID}_chunk_${i}.bcf
                bcftools index ${SAMPLE_ID}_chunk_${i}.bcf
                i=$(( ${i} + 1 ))
            done < ~{split_for_bcftools_merge_csv}
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        # Localizing frequent and infrequent BCFs
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir}/frequent.'bcf*' ~{remote_indir}/infrequent.'bcf*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading the frequent and infrequent BCFs. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        # Re-genotyping every sample assigned to this task
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            SEX=$(echo ${LINE} | cut -d , -f 2)
            TEST=$(gsutil ls ~{remote_outdir}/${SAMPLE_ID}_chunk_0.bcf && echo "" || echo "ERROR")
            if [ ${TEST} = "ERROR" ]; then
                # Proceeding only if the sample has not already been regenotyped
                LocalizeSample ${SAMPLE_ID} ${LINE}
                date
                bcftools view --samples ${SAMPLE_ID} --output-type u infrequent.bcf | bcftools view --drop-genotypes --include 'COUNT(GT="alt")>0' --output-type b > ${SAMPLE_ID}_infrequent.bcf
                date
                bcftools index ${SAMPLE_ID}_infrequent.bcf
                ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --rm-dups exact --output-type z frequent.bcf ${SAMPLE_ID}_infrequent.bcf > ${SAMPLE_ID}_personalized.vcf.gz
                ${TIME_COMMAND} tabix -f ${SAMPLE_ID}_personalized.vcf.gz
                Kanpig ${SAMPLE_ID} ${SEX}
                ChunkKanpig ${SAMPLE_ID}
                while : ; do
                    TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${SAMPLE_ID}_chunk_'*.bcf*' ${SAMPLE_ID}_kanpig_nrecords.txt ~{remote_outdir}/ && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error uploading chunks for sample ${SAMPLE_ID}. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                DelocalizeSample ${SAMPLE_ID}
                ls -laht
            fi
        done < chunk.csv
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
