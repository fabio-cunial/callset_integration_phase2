version 1.0


# Annotates ultralong records with Sniffles and cuteFC features, and extracts TP
# records for XGBoost.
#
workflow UltralongAnnotate {
    input {
        File chunk_csv
        String remote_outdir
        
        File reference_fa
        File reference_fai
        File ultralong_training_resource_bed
        File ultralong_training_resource_vcf_gz
        File ultralong_training_resource_tbi
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong"
        Int preemptible_number = 4
    }
    parameter_meta {
        chunk_csv: "Format: ID,bai,bam,csi,bcf"
        ultralong_training_resource_vcf_gz: "This should be a training resource specifically built for ultralong records. The standard training resource excludes such records."
    }
    
    call Impl {
        input:
            chunk_csv = chunk_csv,
            remote_outdir = remote_outdir,
    
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ultralong_training_resource_bed = ultralong_training_resource_bed,
            ultralong_training_resource_vcf_gz = ultralong_training_resource_vcf_gz,
            ultralong_training_resource_tbi = ultralong_training_resource_tbi,
    
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    
    output {
    }
}


#
task Impl {
    input {
        File chunk_csv
        String remote_outdir
        
        File reference_fa
        File reference_fai
        File ultralong_training_resource_bed
        File ultralong_training_resource_vcf_gz
        File ultralong_training_resource_tbi
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong"
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 128
        Int preemptible_number
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))


        
        # ----------------------- Steps of the pipeline ------------------------
        
        # @param 
        # $2 A row of `chunk_csv`.
        #
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 2)
            ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 3)
            ULTRALONG_CSI=$(echo ${LINE} | cut -d , -f 4)
            ULTRALONG_BCF=$(echo ${LINE} | cut -d , -f 5)
            
            date 1>&2
            gcloud storage cp ${ALIGNED_BAM} ./${SAMPLE_ID}.bam
            date 1>&2
            gcloud storage cp ${ALIGNED_BAI} ./${SAMPLE_ID}.bam.bai
            gcloud storage cp ${ULTRALONG_BCF} ./${SAMPLE_ID}.bcf
            gcloud storage cp ${ULTRALONG_CSI} ./${SAMPLE_ID}.csi
            
            # Converting to .vcf.gz for the genotypers
            bcftools view --threads ${N_THREADS} --output-type z ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}.vcf.gz
            bcftools index --threads ${N_THREADS} -t ${SAMPLE_ID}.vcf.gz
            rm -f ${SAMPLE_ID}.bcf*
            
            # Checking the integrity of the BAM
            ${TIME_COMMAND} samtools quickcheck -v ${SAMPLE_ID}.bam && echo "" || echo "ERROR: the BAM is corrupted."
        }
        
        
        # Deletes all and only the files downloaded by `LocalizeSample()`.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}*.bam* ${SAMPLE_ID}*.bcf* ${SAMPLE_ID}*.vcf.gz*
        }
        
        
        # Ensures the VCF is correctly formatted for the genotypers
        #
        function CanonizeVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
    
            #java -cp ~{docker_dir} FixUltralongRecords ${INPUT_VCF_GZ} ~{reference_fai} | bgzip --compress-level 1 > ${SAMPLE_ID}_canonized.vcf.gz
            #rm -f ${INPUT_VCF_GZ}* ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_canonized.vcf.gz
            
            # 1. Fixing symbolic records
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G FixSymbolicRecords ${INPUT_VCF_GZ} ~{reference_fa} > ${SAMPLE_ID}_out.vcf
            rm -f ${INPUT_VCF_GZ}.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 2. Fixing REF
            ${TIME_COMMAND} bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 3. Cleaning REF, ALT, QUAL, FILTER.
            # - REF and ALT must be uppercase for XGBoost scoring downstream to
            #   work.
            # - We force every record to PASS, to rule out any filter-dependent
            #   effect in downstream tools.
            DEFAULT_QUAL="60"   # Arbitrary
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanRefAltQual ${SAMPLE_ID}_in.vcf ${DEFAULT_QUAL} > ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 4. Removing END, since its values may be inconsistent and make
            # GATK crash downstream.
            ${TIME_COMMAND} bcftools annotate --remove INFO/END --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -t ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_canonized.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_canonized.vcf.gz.tbi
        }
        
        
        # Remark: the procedure stores in a TSV just the features created by the
        # genotyper (the re-genotyped VCF is not saved). In this way we do not
        # care if the genotyper removes fields from the input VCF.
        #
        function Sniffles() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local ALIGNMENTS_BAM=$3
            
            source /opt/sniffles_env/bin/activate
            sniffles --version 1>&2
            ${TIME_COMMAND} sniffles --threads ${N_THREADS} --input ${ALIGNMENTS_BAM} --reference ~{reference_fa} --genotype-vcf ${INPUT_VCF_GZ} --vcf ${SAMPLE_ID}_out.vcf
            deactivate
            mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            bcftools query --format '%CHROM\t%POS\t%ID\t[%GT]\t[%GQ]\t[%DR]\t[%DV]\n' ${SAMPLE_ID}_in.vcf.gz | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                GT_COUNT=-1; \
                if ($4=="0/0" || $4=="0|0" || $4=="./."  || $4==".|." || $4=="./0" || $4==".|0" || $4=="0/." || $4=="0|." || $4=="0" || $4==".") GT_COUNT=0; \
                else if ($4=="0/1" || $4=="0|1" || $4=="1/0" || $4=="1|0" || $4=="./1" || $4==".|1" || $4=="1/." || $4=="1|." || $4=="1") GT_COUNT=1; \
                else if ($4=="1/1" || $4=="1|1") GT_COUNT=2; \
                \
                if ($5==".") GQ=-1; \
                else GQ=$5; \
                \
                if ($6==".") DR=-1; \
                else DR=$6; \
                \
                if ($7==".") DV=-1; \
                else DV=$7; \
                \
                printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\n",$1,$2,$3,GT_COUNT,GQ,DR,DV); \
            }' | bgzip -c > sniffles_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz*
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 sniffles_annotations.tsv.gz
            echo '##INFO=<ID=SNIFFLES_GT_COUNT,Number=1,Type=Integer,Description="Sniffles GT converted to an integer in {0,1,2}.">' > sniffles_header.txt
            echo '##INFO=<ID=SNIFFLES_GQ,Number=1,Type=Integer,Description="Genotype quality according to sniffles">' >> sniffles_header.txt
            echo '##INFO=<ID=SNIFFLES_DR,Number=1,Type=Integer,Description="Number of reference reads according to sniffles">' >> sniffles_header.txt
            echo '##INFO=<ID=SNIFFLES_DV,Number=1,Type=Integer,Description="Number of variant reads according to sniffles">' >> sniffles_header.txt
            SNIFFLES_COLUMNS='CHROM,POS,~ID,INFO/SNIFFLES_GT_COUNT,INFO/SNIFFLES_GQ,INFO/SNIFFLES_DR,INFO/SNIFFLES_DV'
        }
        
        
        # Remark: the procedure stores in a TSV just the features created by the
        # genotyper (the re-genotyped VCF is not saved). In this way we do not
        # care if the genotyper removes fields from the input VCF (cuteFC does,
        # for example).
        #
        function Cutefc() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local ALIGNMENTS_BAM=$3
                
            mkdir ./cutefc_dir/
            source /opt/cutefc_env/bin/activate
            cuteFC --version 1>&2
            ${TIME_COMMAND} cuteFC --threads ${N_THREADS} --genotype --max_size -1 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -Ivcf ${INPUT_VCF_GZ} ${ALIGNMENTS_BAM} ~{reference_fa} ${SAMPLE_ID}_out.vcf ./cutefc_dir
            deactivate
            rm -rf ./cutefc_dir ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            bcftools query --format '%CHROM\t%POS\t%ID\t[%GT]\t[%GQ]\t[%DR]\t[%DV]\t[%PL]\t%INFO/CIPOS\t%INFO/CILEN\t%INFO/RE\t%INFO/STRAND\n' ${SAMPLE_ID}_in.vcf.gz | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                GT_COUNT=-1; \
                if ($4=="0/0" || $4=="0|0" || $4=="./."  || $4==".|." || $4=="./0" || $4==".|0" || $4=="0/." || $4=="0|." || $4=="0" || $4==".") GT_COUNT=0; \
                else if ($4=="0/1" || $4=="0|1" || $4=="1/0" || $4=="1|0" || $4=="./1" || $4==".|1" || $4=="1/." || $4=="1|." || $4=="1") GT_COUNT=1; \
                else if ($4=="1/1" || $4=="1|1") GT_COUNT=2; \
                \
                if ($5==".") GQ=-1; \
                else GQ=$5; \
                \
                if ($6==".") DR=-1; \
                else DR=$6; \
                \
                if ($7==".") DV=-1; \
                else DV=$7; \
                \
                PL_1=-1; PL_2=-1; PL_3=-1; \
                p=0; \
                for (i=1; i<=length($8); i++) { \
                    if (substr($8,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    PL_1=substr($8,1,p-1); \
                    q=0; \
                    for (i=p+1; i<=length($8); i++) { \
                        if (substr($8,i,1)==",") { q=i; break; } \
                    } \
                    if (q>0) { \
                        PL_2=substr($8,p+1,q-1-p); \
                        PL_3=substr($8,q+1); \
                    } \
                    else { PL_2=substr($8,p+1); } \
                } \
                else { PL_1=$8; }
                \
                CIPOS_1=-1; CIPOS_2=-1; \
                p=0; \
                for (i=1; i<=length($9); i++) { \
                    if (substr($9,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    CIPOS_1=substr($9,1,p-1); \
                    CIPOS_2=substr($9,p+1); \
                } \
                else { CIPOS_1=$9; } \
                \
                CILEN_1=-1; CILEN_2=-1; \
                p=0; \
                for (i=1; i<=length($10); i++) { \
                    if (substr($10,i,1)==",") { p=i; break; } \
                } \
                if (p>0) { \
                    CILEN_1=substr($10,1,p-1); \
                    CILEN_2=substr($10,p+1); \
                } \
                else { CILEN_1=$10; } \
                \
                if ($11==".") RE=-1; \
                else RE=$11; \
                \
                STRAND=-1; \
                if ($12=="--") STRAND=0; \
                else if ($12=="-+") STRAND=1; \
                else if ($12=="+-") STRAND=2; \
                else if ($12=="++") STRAND=3; \
                \
                printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n",$1,$2,$3,GT_COUNT,GQ,DR,DV,PL_1,PL_2,PL_3,CIPOS_1,CIPOS_2,CILEN_1,CILEN_2,RE,STRAND); \
            }' | bgzip -c > cutefc_annotations.tsv.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz*
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 cutefc_annotations.tsv.gz
            echo '##INFO=<ID=CUTEFC_GT_COUNT,Number=1,Type=Integer,Description="Cutefc GT converted to an integer in {0,1,2}.">' > cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_GQ,Number=1,Type=Integer,Description="Genotype quality according to cutefc">' >> cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_DR,Number=1,Type=Integer,Description="High-quality reference reads according to cutefc">' >> cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_DV,Number=1,Type=Integer,Description="High-quality variant reads according to cutefc">' >> cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_PL_1,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_PL_2,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_PL_3,Number=1,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer according to cutefc">' >> cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_CIPOS_1,Number=1,Type=Integer,Description="Confidence interval around POS for imprecise variants according to cutefc">' >> cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_CIPOS_2,Number=1,Type=Integer,Description="Confidence interval around POS for imprecise variants according to cutefc">' >> cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_CILEN_1,Number=1,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends according to cutefc">' >> cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_CILEN_2,Number=1,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends according to cutefc">' >> cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_RE,Number=1,Type=Integer,Description="Number of read support this record according to cutefc">' >> cutefc_header.txt
            echo '##INFO=<ID=CUTEFC_STRAND,Number=1,Type=Integer,Description="Cutefc strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--) converted to an integer in {0,1,2,3}.">' >> cutefc_header.txt
            CUTEFC_COLUMNS='CHROM,POS,~ID,INFO/CUTEFC_GT_COUNT,INFO/CUTEFC_GQ,INFO/CUTEFC_DR,INFO/CUTEFC_DV,INFO/CUTEFC_PL_1,INFO/CUTEFC_PL_2,INFO/CUTEFC_PL_3,INFO/CUTEFC_CIPOS_1,INFO/CUTEFC_CIPOS_2,INFO/CUTEFC_CILEN_1,INFO/CUTEFC_CILEN_2,INFO/CUTEFC_RE,INFO/CUTEFC_STRAND'
        }
        
        
        # Copies the annotations of all genotypers to a single VCF.
        #
        # Remark: SUPP_* fields from intra-sample truvari are already in INFO,
        # thanks to the pipeline steps upstream.
        #
        function Annotate() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local OUTPUT_VCF_GZ=$3
            
            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_in.vcf.gz
            mv ${INPUT_VCF_GZ}.tbi ${SAMPLE_ID}_in.vcf.gz.tbi
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations sniffles_annotations.tsv.gz --header-lines sniffles_header.txt --columns ${SNIFFLES_COLUMNS} --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations cutefc_annotations.tsv.gz --header-lines cutefc_header.txt --columns ${CUTEFC_COLUMNS} --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            rm -f ${INPUT_VCF_GZ}*
            
            mv ${SAMPLE_ID}_in.vcf.gz ${OUTPUT_VCF_GZ}
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${OUTPUT_VCF_GZ}.tbi
        }
        
        
        # Remark: we keep sequence similarity on, even if the records are
        # ultralong, to be as specific as possible.
        #
        cat << 'END' > truvari_bench.sh
#!/bin/bash
TRAINING_RESOURCE_VCF_GZ=$1
INFINITY=$2
INPUT_VCF=$3

CHUNK_ID=${INPUT_VCF#chunk_}
CHUNK_ID=${CHUNK_ID%.vcf}
${TIME_COMMAND} truvari bench -b ${TRAINING_RESOURCE_VCF_GZ} -c ${INPUT_VCF_GZ} --dup-to-ins --max-resolve ${INFINITY} --no-roll --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick single -o truvari_${CHUNK_ID}/
END
        chmod +x truvari_bench.sh
                
        
        # Extracts every record that has a stringent `truvari bench` match with
        # some records in the resource. This is parallelized by record.
        #
        function GetTrainingRecords() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            # Storing every record in a separate VCF
            rm -f tasks.wsv
            bcftools view --header-only ${INPUT_VCF_GZ} > header.txt
            bcftools view --no-header ${INPUT_VCF_GZ} | split -d -a 4 -l 1 chunk_
            for FILE in $( ls chunk_* ); do
                cat header.txt ${FILE} > ${FILE}.vcf
                rm -f ${FILE}
                echo "${FILE}.vcf" >> tasks.wsv
            done
            rm -f header.txt
            
            # Running every record in parallel
            ${TIME_COMMAND} xargs --arg-file=tasks.wsv --max-lines=1 --max-procs=${N_THREADS} ./truvari_bench.sh ~{ultralong_training_resource_vcf_gz} ${INFINITY}
            
            # Concatenating TPs
            rm -f file_list.txt
            while read FILE; do
                CHUNK_ID=${FILE#chunk_}
                CHUNK_ID=${CHUNK_ID%.vcf}
                echo truvari_${CHUNK_ID}/tp-comp.vcf.gz >> file_list.txt
            done < tasks.wsv
            ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list file_list.txt --output-type z --output ${SAMPLE_ID}_training.vcf.gz
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_training.vcf.gz
            
            # Removing temporary files
            rm -rf ./truvari_*/ ./chunk_*.vcf file_list.txt
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        truvari --help 1>&2
        df -h 1>&2
        
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            # Skipping the sample if it has already been processed
            TEST=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
        
            # Annotating and marking training records
            LocalizeSample ${SAMPLE_ID} ${LINE}
            CanonizeVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz
            Sniffles ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}.bam
            Cutefc ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}.bam
            Annotate ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}_annotated.vcf.gz
            GetTrainingRecords ${SAMPLE_ID} ${SAMPLE_ID}_annotated.vcf.gz
        
            # Uploading
            gcloud storage mv ${SAMPLE_ID}_annotated.vcf.'gz*' ${SAMPLE_ID}_training.vcf.'gz*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done < ~{chunk_csv}
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
