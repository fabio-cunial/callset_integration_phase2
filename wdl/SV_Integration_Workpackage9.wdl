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
        File autosomes_bed
        
        String kanpig_params_cohort = "--neighdist 500 --gpenalty 0.04 --hapsim 0.97"
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
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
            autosomes_bed = autosomes_bed,
            
            kanpig_params_cohort = kanpig_params_cohort,
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, 8 logical cores, 8GB RAM,
# chr6,X,Y:
#
# TOOL                      CPU     RAM     TIME
# BAM download                               3m
# bcftools view --samples                    3m
# bcftools concat           600%     34M     5s
# kanpig                    400%    1.5G     3m
# bcftools sort             100%    1.5G    40s
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
        File autosomes_bed
        
        String kanpig_params_cohort
        String docker_image
        
        Int n_cpu = 6
        Int ram_size_gb = 8
        Int disk_size_gb = 50
        Int preemptible_number = 4
    }
    parameter_meta {
        disk_size_gb: "50GB is enough for many 15x samples, but for ~1500 samples it is not sufficient. 200GB is enough for all 15x and 30x samples."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        GSUTIL_DELAY_S="600"
        export RUST_BACKTRACE="full"
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 3)
            ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 4)
            
            # Failing immediately if the BAM is too large. Otherwise the VM
            # may get stuck forever, and this is even worse with preemption.
            AVAILABLE_GB=$(df -h | grep "cromwell_root" | tr -s ' ' | cut -d ' ' -f 4)
            AVAILABLE_GB=${AVAILABLE_GB%G}
            AVAILABLE_GB=${AVAILABLE_GB%.*}
            BAM_GB=$(gsutil ls -lh ${ALIGNED_BAM} | head -n 1 | sed 's/^[ ]*//' | tr -s ' ' | cut -d ' ' -f 1)
            BAM_GB=${BAM_GB%.*}
            SLACK_GB="5"
            BAM_GB=$(( ${BAM_GB} + ${SLACK_GB} ))
            if [ ${BAM_GB} -gt ${AVAILABLE_GB} ]; then
                echo "ERROR: the BAM is larger than the available disk space. BAM size + slack: ${BAM_GB}GB. Available disk: ${AVAILABLE_GB}GB."
                exit 1
            fi
            
            # Downloading
            date 1>&2
            gcloud storage cp ${ALIGNED_BAM} ./${SAMPLE_ID}_aligned.bam
            date 1>&2
            gcloud storage cp ${ALIGNED_BAI} ./${SAMPLE_ID}_aligned.bam.bai
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
            
            # Remark: kanpig needs --sizemin >= --kmer
            ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --sizemin 10 --sizemax ${INFINITY} ~{kanpig_params_cohort} --reference ~{reference_fa} --ploidy-bed ${PLOIDY_BED} --input ${SAMPLE_ID}_personalized.vcf.gz --reads ${SAMPLE_ID}_aligned.bam --out ${SAMPLE_ID}_out.vcf --sample ${SAMPLE_ID}
            rm -f ${SAMPLE_ID}_personalized.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # Sorting
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_kanpig.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_kanpig.vcf.gz.tbi
            
            # Printing debug information
            N_RECORDS=$(bcftools index --nrecords ${SAMPLE_ID}_kanpig.vcf.gz)
            N_PRESENT_RECORDS=$( bcftools query --format '%ID' --include 'GT="alt"' ${SAMPLE_ID}_kanpig.vcf.gz | wc -l )
            PERCENT=$( echo "scale=2; 100 * ${N_PRESENT_RECORDS} / ${N_RECORDS}" | bc )
            echo "${N_PRESENT_RECORDS},${N_RECORDS},${PERCENT},Number of records that are marked as ALT by kanpig" >> ${SAMPLE_ID}_kanpig.csv
            N_HETS_IN_AUTOSOMES=$( bcftools query --format '%ID' --include 'GT="het"' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l )
            N_PRESENT_RECORDS_IN_AUTOSOMES=$( bcftools query --format '%ID' --include 'GT="alt"' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l )
            PERCENT=$( echo "scale=2; 100 * ${N_HETS_IN_AUTOSOMES} / ${N_PRESENT_RECORDS_IN_AUTOSOMES}" | bc )
            echo "${N_HETS_IN_AUTOSOMES},${N_PRESENT_RECORDS_IN_AUTOSOMES},${PERCENT},Number of records in autosomes that are marked as HET by kanpig" >> ${SAMPLE_ID}_kanpig.csv
            ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigWindows ${SAMPLE_ID}_kanpig.vcf.gz | bgzip > ${SAMPLE_ID}_kanpig.bed.gz
        }
        
        
        # Splits the re-genotyped VCF into chunks for bcftools merge downstream.
        # 
        function ChunkAndUpload() {
            local SAMPLE_ID=$1
            
            i="0"
            while read INTERVAL; do
                echo ${INTERVAL} | tr ',' '\t' > ${SAMPLE_ID}.bed
                ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --regions-file ${SAMPLE_ID}.bed --regions-overlap pos --output-type b ${SAMPLE_ID}_kanpig.vcf.gz --output ${SAMPLE_ID}_chunk_${i}.bcf
                bcftools index --threads ${N_THREADS} ${SAMPLE_ID}_chunk_${i}.bcf
                gcloud storage cp ${SAMPLE_ID}_chunk_${i}.bcf ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf
                gcloud storage cp ${SAMPLE_ID}_chunk_${i}.bcf.csi ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf.csi
                i=$(( ${i} + 1 ))
            done < ~{split_for_bcftools_merge_csv}
            gcloud storage cp ${SAMPLE_ID}_kanpig.bed.gz ${SAMPLE_ID}_kanpig.csv ~{remote_outdir}/
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        ~{docker_dir}/kanpig --version 1>&2
        
        # Localizing frequent and infrequent BCFs
        date 1>&2
        gcloud storage cp ~{remote_indir}/frequent.'bcf*' ~{remote_indir}/infrequent.'bcf*' .
        date 1>&2
        
        # Re-genotyping every sample assigned to this task
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            SEX=$(echo ${LINE} | cut -d , -f 2)
            
            # Skipping the sample if it has already been processed
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            # Re-genotyping
            LocalizeSample ${SAMPLE_ID} ${LINE}
            date 1>&2
            bcftools view --samples ${SAMPLE_ID} --output-type u infrequent.bcf | bcftools view --include 'GT="alt"' --drop-genotypes --output-type b --output ${SAMPLE_ID}_infrequent.bcf
            date 1>&2
            bcftools index --threads ${N_THREADS} ${SAMPLE_ID}_infrequent.bcf
            ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --rm-dups exact --output-type z frequent.bcf ${SAMPLE_ID}_infrequent.bcf --output ${SAMPLE_ID}_personalized.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_personalized.vcf.gz
            Kanpig ${SAMPLE_ID} ${SEX}
            ChunkAndUpload ${SAMPLE_ID}
            
            # Marking the sample as completed
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done < chunk.csv
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
