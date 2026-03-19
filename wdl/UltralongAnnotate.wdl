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
        
        File feature_extraction_py
        
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
    
            feature_extraction_py = feature_extraction_py,
    
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    
    output {
    }
}


# TOOL                      CPU     RAM     TIME
# sniffles           
# cutefc       
# bcftools query
# bcftools annotate
# truvari bench      
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
        
        File feature_extraction_py
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong"
        Int n_cpu = 4
        Int ram_size_gb = 32
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
        
        
        # Ensures that the VCF is correctly formatted for the genotypers
        #
        function CanonizeVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            # 1. Fixing END
            java -cp ~{docker_dir} FixUltralongRecords ${INPUT_VCF_GZ} ~{reference_fai} | bgzip --compress-level 1 > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${INPUT_VCF_GZ}* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            # 2. Adding sequence to symbolic records
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G FixSymbolicRecords ${SAMPLE_ID}_in.vcf.gz ~{reference_fa} > ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 3. Fixing REF
            ${TIME_COMMAND} bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # 4. Cleaning REF, ALT, QUAL, FILTER.
            # - REF and ALT must be uppercase for XGBoost scoring downstream to
            #   work.
            # - We force every record to PASS, to rule out any filter-dependent
            #   effect in downstream tools.
            DEFAULT_QUAL="60"   # Arbitrary
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanRefAltQual ${SAMPLE_ID}_in.vcf ${DEFAULT_QUAL} > ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            bcftools view --threads ${N_THREADS} --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_canonized.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_canonized.vcf.gz
        }
        
        
        # Remark: the procedure stores in a TSV just the features created by the
        # genotyper (the re-genotyped VCF is not saved). In this way we do not
        # care if the genotyper removes fields from the input VCF.
        #
        # TOOL                     CPU%         RAM        TIME     NOTES
        # sniffles 2.6.3 --vcf     100%         3G         20s      Empty output
        # sniffles 2.7.3 --vcf                                      Crashes
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
        # TOOL                  CPU%      RAM        TIME       NOTES
        # cuteFC 1.0.2 -Ivcf    120%      3G         44s        All zeros in out
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
        
        
        cat << 'END' > lrcaller.sh
#!/bin/bash
set -euxo pipefail

SAMPLE_ID=$1
INPUT_VCF_GZ=$2
ALIGNMENTS_BAM=$3
BREAKPOINT=$4
REFERENCE_FA=$5
N_THREADS=$6

if [ ${BREAKPOINT} -eq 0 ]; then
    BREAKPOINT_FLAG=" "
else
    BREAKPOINT_FLAG="--right_breakpoint"
fi
lrcaller --number_of_threads ${N_THREADS} ${BREAKPOINT_FLAG} --dyn-w-size --fa ${REFERENCE_FA} ${ALIGNMENTS_BAM} ${INPUT_VCF_GZ} ${SAMPLE_ID}_out.vcf 2> /dev/null
END
        chmod +x lrcaller.sh
        
        
        # Remark: the procedure stores in a TSV just the features created by the
        # genotyper (the re-genotyped VCF is not saved). In this way we do not
        # care if the genotyper removes fields from the input VCF.
        #
        # TOOL                           CPU%         RAM        TIME    NOTES
        # lrcaller                       300%       30.5G         10m
        # lrcaller --right_breakpoint                                    Crashes
        #
        function Lrcaller() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local ALIGNMENTS_BAM=$3
            local BREAKPOINT=$4
            
            bcftools view --threads ${N_THREADS} --drop-genotypes --output-type z ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_in.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            #lrcaller --version 1>&2
            ${TIME_COMMAND} ./lrcaller.sh ${SAMPLE_ID} ${SAMPLE_ID}_in.vcf.gz ${ALIGNMENTS_BAM} ${BREAKPOINT} ~{reference_fa} ${N_THREADS}
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            if [ ${BREAKPOINT} -eq 0 ]; then
                SUFFIX="left"
            else
                SUFFIX="right"
            fi
            grep '^[^#]' ${SAMPLE_ID}_in.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                printf("%s",$1); \
                for (i=2; i<=3; i++) printf("\t%s",$i); \
                for (i=10; i<=14; i++) { \
                    gsub(/[:,]/,"\t",$i); \
                    printf("\t%s",$i); \
                } \
                printf("\n"); \
            }' | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                GT_COUNT=-1; \
                if ($4=="0/0" || $4=="0|0" || $4=="./."  || $4==".|." || $4=="./0" || $4==".|0" || $4=="0/." || $4=="0|." || $4=="0" || $4==".") GT_COUNT=0; \
                else if ($4=="0/1" || $4=="0|1" || $4=="1/0" || $4=="1|0" || $4=="./1" || $4==".|1" || $4=="1/." || $4=="1|." || $4=="1") GT_COUNT=1; \
                else if ($4=="1/1" || $4=="1|1") GT_COUNT=2; \
                $4=GT_COUNT; \
                \
                GT_COUNT=-1; \
                if ($14=="0/0" || $14=="0|0" || $14=="./."  || $14==".|." || $14=="./0" || $14==".|0" || $14=="0/." || $14=="0|." || $14=="0" || $14==".") GT_COUNT=0; \
                else if ($14=="0/1" || $14=="0|1" || $14=="1/0" || $14=="1|0" || $14=="./1" || $14==".|1" || $14=="1/." || $14=="1|." || $14=="1") GT_COUNT=1; \
                else if ($14=="1/1" || $14=="1|1") GT_COUNT=2; \
                $14=GT_COUNT; \
                \
                GT_COUNT=-1; \
                if ($24=="0/0" || $24=="0|0" || $24=="./."  || $24==".|." || $24=="./0" || $24==".|0" || $24=="0/." || $24=="0|." || $24=="0" || $24==".") GT_COUNT=0; \
                else if ($24=="0/1" || $24=="0|1" || $24=="1/0" || $24=="1|0" || $24=="./1" || $24==".|1" || $24=="1/." || $24=="1|." || $24=="1") GT_COUNT=1; \
                else if ($24=="1/1" || $24=="1|1") GT_COUNT=2; \
                $24=GT_COUNT; \
                \
                GT_COUNT=-1; \
                if ($34=="0/0" || $34=="0|0" || $34=="./."  || $34==".|." || $34=="./0" || $34==".|0" || $34=="0/." || $34=="0|." || $34=="0" || $34==".") GT_COUNT=0; \
                else if ($34=="0/1" || $34=="0|1" || $34=="1/0" || $34=="1|0" || $34=="./1" || $34==".|1" || $34=="1/." || $34=="1|." || $34=="1") GT_COUNT=1; \
                else if ($34=="1/1" || $34=="1|1") GT_COUNT=2; \
                $34=GT_COUNT; \
                \
                GT_COUNT=-1; \
                if ($44=="0/0" || $44=="0|0" || $44=="./."  || $44==".|." || $44=="./0" || $44==".|0" || $44=="0/." || $44=="0|." || $44=="0" || $44==".") GT_COUNT=0; \
                else if ($44=="0/1" || $44=="0|1" || $44=="1/0" || $44=="1|0" || $44=="./1" || $44==".|1" || $44=="1/." || $44=="1|." || $44=="1") GT_COUNT=1; \
                else if ($44=="1/1" || $44=="1|1") GT_COUNT=2; \
                $44=GT_COUNT; \
                \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' | bgzip -c > lrcaller_annotations_${SUFFIX}.tsv.gz
            zcat lrcaller_annotations_${SUFFIX}.tsv.gz | head -n 10 1>&2 || echo "0"
            rm -f ${SAMPLE_ID}_in.vcf
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 lrcaller_annotations_${SUFFIX}.tsv.gz
            echo '##INFO=<ID=GTCOUNT1_'${SUFFIX}',Number=1,Type=Integer,Description="Genotype count">' > lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=GTCOUNT2_'${SUFFIX}',Number=1,Type=Integer,Description="Genotype count">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=GTCOUNT3_'${SUFFIX}',Number=1,Type=Integer,Description="Genotype count">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=GTCOUNT4_'${SUFFIX}',Number=1,Type=Integer,Description="Genotype count">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=GTCOUNT5_'${SUFFIX}',Number=1,Type=Integer,Description="Genotype count">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=AD11_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=AD12_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=AD13_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=AD21_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=AD22_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=AD23_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=AD31_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=AD32_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=AD33_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=AD41_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=AD42_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=AD43_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=AD51_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=AD52_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=AD53_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=VA11_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=VA12_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=VA13_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=VA21_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=VA22_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=VA23_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=VA31_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=VA32_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=VA33_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=VA41_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=VA42_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=VA43_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=VA51_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=VA52_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=VA53_'${SUFFIX}',Number=1,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=PL11_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=PL12_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=PL13_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=PL21_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=PL22_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=PL23_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=PL31_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=PL32_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=PL33_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=PL41_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=PL42_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=PL43_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            
            echo '##INFO=<ID=PL51_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=PL52_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            echo '##INFO=<ID=PL53_'${SUFFIX}',Number=1,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> lrcaller_header_${SUFFIX}.txt
            
            if [ ${BREAKPOINT} -eq 0 ]; then
                LRCALLER_COLUMNS_LEFT='CHROM,POS,~ID,INFO/GTCOUNT1_'${SUFFIX}',INFO/AD11_'${SUFFIX}',INFO/AD12_'${SUFFIX}',INFO/AD13_'${SUFFIX}',INFO/VA11_'${SUFFIX}',INFO/VA12_'${SUFFIX}',INFO/VA13_'${SUFFIX}',INFO/PL11_'${SUFFIX}',INFO/PL12_'${SUFFIX}',INFO/PL13_'${SUFFIX}',INFO/GTCOUNT2_'${SUFFIX}',INFO/AD21_'${SUFFIX}',INFO/AD22_'${SUFFIX}',INFO/AD23_'${SUFFIX}',INFO/VA21_'${SUFFIX}',INFO/VA22_'${SUFFIX}',INFO/VA23_'${SUFFIX}',INFO/PL21_'${SUFFIX}',INFO/PL22_'${SUFFIX}',INFO/PL23_'${SUFFIX}',INFO/GTCOUNT3_'${SUFFIX}',INFO/AD31_'${SUFFIX}',INFO/AD32_'${SUFFIX}',INFO/AD33_'${SUFFIX}',INFO/VA31_'${SUFFIX}',INFO/VA32_'${SUFFIX}',INFO/VA33_'${SUFFIX}',INFO/PL31_'${SUFFIX}',INFO/PL32_'${SUFFIX}',INFO/PL33_'${SUFFIX}',INFO/GTCOUNT4_'${SUFFIX}',INFO/AD41_'${SUFFIX}',INFO/AD42_'${SUFFIX}',INFO/AD43_'${SUFFIX}',INFO/VA41_'${SUFFIX}',INFO/VA42_'${SUFFIX}',INFO/VA43_'${SUFFIX}',INFO/PL41_'${SUFFIX}',INFO/PL42_'${SUFFIX}',INFO/PL43_'${SUFFIX}',INFO/GTCOUNT5_'${SUFFIX}',INFO/AD51_'${SUFFIX}',INFO/AD52_'${SUFFIX}',INFO/AD53_'${SUFFIX}',INFO/VA51_'${SUFFIX}',INFO/VA52_'${SUFFIX}',INFO/VA53_'${SUFFIX}',INFO/PL51_'${SUFFIX}',INFO/PL52_'${SUFFIX}',INFO/PL53_'${SUFFIX}
            else
                LRCALLER_COLUMNS_RIGHT='CHROM,POS,~ID,INFO/GTCOUNT1_'${SUFFIX}',INFO/AD11_'${SUFFIX}',INFO/AD12_'${SUFFIX}',INFO/AD13_'${SUFFIX}',INFO/VA11_'${SUFFIX}',INFO/VA12_'${SUFFIX}',INFO/VA13_'${SUFFIX}',INFO/PL11_'${SUFFIX}',INFO/PL12_'${SUFFIX}',INFO/PL13_'${SUFFIX}',INFO/GTCOUNT2_'${SUFFIX}',INFO/AD21_'${SUFFIX}',INFO/AD22_'${SUFFIX}',INFO/AD23_'${SUFFIX}',INFO/VA21_'${SUFFIX}',INFO/VA22_'${SUFFIX}',INFO/VA23_'${SUFFIX}',INFO/PL21_'${SUFFIX}',INFO/PL22_'${SUFFIX}',INFO/PL23_'${SUFFIX}',INFO/GTCOUNT3_'${SUFFIX}',INFO/AD31_'${SUFFIX}',INFO/AD32_'${SUFFIX}',INFO/AD33_'${SUFFIX}',INFO/VA31_'${SUFFIX}',INFO/VA32_'${SUFFIX}',INFO/VA33_'${SUFFIX}',INFO/PL31_'${SUFFIX}',INFO/PL32_'${SUFFIX}',INFO/PL33_'${SUFFIX}',INFO/GTCOUNT4_'${SUFFIX}',INFO/AD41_'${SUFFIX}',INFO/AD42_'${SUFFIX}',INFO/AD43_'${SUFFIX}',INFO/VA41_'${SUFFIX}',INFO/VA42_'${SUFFIX}',INFO/VA43_'${SUFFIX}',INFO/PL41_'${SUFFIX}',INFO/PL42_'${SUFFIX}',INFO/PL43_'${SUFFIX}',INFO/GTCOUNT5_'${SUFFIX}',INFO/AD51_'${SUFFIX}',INFO/AD52_'${SUFFIX}',INFO/AD53_'${SUFFIX}',INFO/VA51_'${SUFFIX}',INFO/VA52_'${SUFFIX}',INFO/VA53_'${SUFFIX}',INFO/PL51_'${SUFFIX}',INFO/PL52_'${SUFFIX}',INFO/PL53_'${SUFFIX}
            fi
        }
        
        
        # This script suffers from the same problem as sniffles and cutesv, i.e.
        # it crashes on reading the BAM and does not write the corresponding
        # annotations:
        #
        # [E::bgzf_read] Read block operation failed with error 3 after 0 of 4
        # bytes
        #
        # So it is likely an issue with pysam. The only nontrivial annotations
        # written are gc_frac, homopolymer_max.
        #
        # TOOL                           CPU%         RAM        TIME    NOTES
        # feature_extraction.py          100%        3.5G         10s
        #
        function FeatureExtraction() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local ALIGNMENTS_BAM=$3
            
            ${TIME_COMMAND} python ~{feature_extraction_py} ${INPUT_VCF_GZ} ${ALIGNMENTS_BAM} ~{reference_fa} 1>&2
            head -n 20 features.csv 1>&2
            N_DEL=$(bcftools query --format '%ID\n' --include 'SVTYPE="DEL"' ${INPUT_VCF_GZ} | wc -l)
            N_INS=$(bcftools query --format '%ID\n' --include 'SVTYPE="INS"' ${INPUT_VCF_GZ} | wc -l)
            N_INV=$(bcftools query --format '%ID\n' --include 'SVTYPE="INV"' ${INPUT_VCF_GZ} | wc -l)
            N_DUP=$(bcftools query --format '%ID\n' --include 'SVTYPE="DUP"' ${INPUT_VCF_GZ} | wc -l)
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
            
            #${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations sniffles_annotations.tsv.gz --header-lines sniffles_header.txt --columns ${SNIFFLES_COLUMNS} --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            #rm -f ${SAMPLE_ID}_in.vcf.gz ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            #${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations cutefc_annotations.tsv.gz --header-lines cutefc_header.txt --columns ${CUTEFC_COLUMNS} --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            #rm -f ${SAMPLE_ID}_in.vcf.gz ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations lrcaller_annotations_left.tsv.gz --header-lines lrcaller_header_left.txt --columns ${LRCALLER_COLUMNS_LEFT} --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            #${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations lrcaller_annotations_right.tsv.gz --header-lines lrcaller_header_right.txt --columns ${LRCALLER_COLUMNS_RIGHT} --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            #rm -f ${SAMPLE_ID}_in.vcf.gz ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${OUTPUT_VCF_GZ}
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${OUTPUT_VCF_GZ}.tbi
        }
        
        
        # Remark: we keep sequence similarity on, even if the records are
        # ultralong, to be as specific as possible. However, this is too slow in
        # practice, even after parallelizing by record.
        #
        cat << 'END' > truvari_bench.sh
#!/bin/bash
set -euxo pipefail

TRAINING_RESOURCE_VCF_GZ=$1
INFINITY=$2
INPUT_VCF_GZ=$3

CHUNK_ID=${INPUT_VCF_GZ#chunk_}
CHUNK_ID=${CHUNK_ID%.vcf.gz}
SVTYPE=$(bcftools query --format '%INFO/SVTYPE\n' ${INPUT_VCF_GZ})
if [ ${SVTYPE} = "INS" ]; then
    PCTSEQ_FLAG="--pctseq 0.9 --no-roll"
else
    PCTSEQ_FLAG="--pctseq 0"
fi
${TIME_COMMAND} truvari bench -b ${TRAINING_RESOURCE_VCF_GZ} -c ${INPUT_VCF_GZ} --dup-to-ins --max-resolve ${INFINITY} --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 ${PCTSEQ_FLAG} --pick single -o truvari_${CHUNK_ID}/
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
            bcftools view --no-header ${INPUT_VCF_GZ} | split -d -a 4 -l 1 - chunk_
            for FILE in chunk_* ; do
                cat header.txt ${FILE} | bgzip > ${FILE}.vcf.gz
                bcftools index -t -f ${FILE}.vcf.gz
                rm -f ${FILE}
                echo "${FILE}.vcf.gz" >> tasks.wsv
            done
            rm -f header.txt
            
            # Running every record in parallel
            ${TIME_COMMAND} xargs --arg-file=tasks.wsv --max-lines=1 --max-procs=${N_THREADS} ./truvari_bench.sh ~{ultralong_training_resource_vcf_gz} ${INFINITY}
            
            # Concatenating TPs
            rm -f file_list.txt
            while read FILE; do
                CHUNK_ID=${FILE#chunk_}
                CHUNK_ID=${CHUNK_ID%.vcf.gz}
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
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
        
            # Annotating and marking training records
            LocalizeSample ${SAMPLE_ID} ${LINE}
            CanonizeVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz
            #Sniffles ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}.bam
            #Cutefc ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}.bam
            Lrcaller ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}.bam 0
            #Lrcaller ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}.bam 1
            Annotate ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}_annotated.vcf.gz
            gcloud storage cp ${SAMPLE_ID}_annotated.vcf.'gz*' ~{remote_outdir}/
            #FeatureExtraction ${SAMPLE_ID} ${SAMPLE_ID}_canonized.vcf.gz ${SAMPLE_ID}.bam
            GetTrainingRecords ${SAMPLE_ID} ${SAMPLE_ID}_annotated.vcf.gz
            gcloud storage cp ${SAMPLE_ID}_training.vcf.'gz*' ~{remote_outdir}/
            
            # Next iteration
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
