version 1.0


# Builds a training resource for XGBoost.
#
workflow SV_Integration_BuildTrainingResource {
    input {
        File dipcall_tsv
        
        Int min_sv_length = 20
        Int max_sv_length = 10000
        
        File reference_fai
        File reference_agp
        
        String remote_outdir
    }
    parameter_meta {
        dipcall_tsv: "Format of each row: ID, DIPCALL_BED, DIPCALL_VCF. We assume that every VCF is sorted and sequence-resolved."
    }
    
    call Impl {
        input:
            dipcall_tsv = dipcall_tsv,
            min_sv_length = min_sv_length,
            max_sv_length = max_sv_length,
            reference_fai = reference_fai,
            reference_agp = reference_agp,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


#
task Impl {
    input {
        File dipcall_tsv
        
        Int min_sv_length
        Int max_sv_length
        
        File reference_fai
        File reference_agp
        
        String remote_outdir
        
        Int n_cpu = 8
        Int ram_size_gb = 32
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
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            DIPCALL_BED=$(echo ${LINE} | cut -d , -f 2)
            DIPCALL_VCF_GZ=$(echo ${LINE} | cut -d , -f 3)
            while : ; do
                TEST=$(gsutil -m cp ${DIPCALL_BED} ./${SAMPLE_ID}.bed && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${DIPCALL_BED}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gsutil -m cp ${DIPCALL_VCF_GZ} ./${SAMPLE_ID}.vcf.gz && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${DIPCALL_VCF_GZ}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            tabix -f ${SAMPLE_ID}.vcf.gz
        }
        
        
        # Deletes all and only the files downloaded by `LocalizeSample()`.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}.bed ${SAMPLE_ID}.vcf.gz*
        }
        
        
        # Returns a BED file that excludes every gap from the AGP file of the
        # reference.
        #
        function GetReferenceGaps() {
            local INPUT_AGP=$1
            local OUTPUT_BED=$2
            
            awk 'BEGIN { FS="\t"; OFS="\t"; } { if ($1=="chr1" || $1=="chr2" || $1=="chr3" || $1=="chr4" || $1=="chr5" || $1=="chr6" || $1=="chr7" || $1=="chr8" || $1=="chr9" || $1=="chr10" || $1=="chr11" || $1=="chr12" || $1=="chr13" || $1=="chr14" || $1=="chr15" || $1=="chr16" || $1=="chr17" || $1=="chr18" || $1=="chr19" || $1=="chr20" || $1=="chr21" || $1=="chr22" || $1=="chrX" || $1=="chrY" || $1=="chrM") print $0 }' ${INPUT_AGP} > in.bed
            awk 'BEGIN { FS="\t"; OFS="\t"; } { if ($5=="N") print $0 }' in.bed > out.bed
            mv out.bed in.bed
            bedtools sort -i in.bed -faidx ~{reference_fai} > out.bed
            mv out.bed in.bed
            bedtools complement -i in.bed -g ~{reference_fai} > out.bed
            mv out.bed in.bed
            
            mv in.bed ${OUTPUT_BED}
        }
        
        
        # Puts in canonical form a raw VCF from dipcall.
        #
        function CanonizeVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local INPUT_TBI=$3
            local MIN_SV_LENGTH=$4
            local MAX_SV_LENGTH=$5
            local NOT_GAPS_BED=$6
            
            
            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_in.vcf.gz
            mv ${INPUT_TBI} ${SAMPLE_ID}_in.vcf.gz.tbi
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Removing SNVs, records that are not marked as present, records
            # with a FILTER, and records with unresolved REF/ALT.
            ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || COUNT(GT="alt")=0 || (FILTER!="PASS" && FILTER!=".") || REF="*" || ALT="*"' --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the dipcall BED
            ${TIME_COMMAND} bcftools filter --regions-file ${SAMPLE_ID}.bed --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Removing records in reference gaps
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Making sure SVLEN and SVTYPE are consistently annotated
            truvari anno svinfo --minsize 1 ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the given length range
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_canonized.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_canonized.vcf.gz.tbi
        }
        
        
        # A simple `bcftools merge` without any additional collapse.
        #
        function Merge() {
            local FILE_LIST_TXT=$1
            
            # Merging
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list ${FILE_LIST_TXT} --output-type v > out.vcf
            mv out.vcf in.vcf
            
            # Enforcing a single sample column with artificial GTs
            bcftools view --header-only in.vcf > header.txt
            N_ROWS=$(wc -l < header.txt)
            head -n $(( ${N_ROWS} - 1 )) header.txt > out.vcf
            tail -n 1 header.txt | awk '{ \
                printf("%s",$1); \
                for (i=2; i<=9; i++) printf("\t%s",$i); \
                printf("\tSAMPLE\n"); \
            }' >> out.vcf
            rm -f header.txt
            bcftools view --no-header in.vcf | awk '{ \
                printf("%s",$1); \
                for (i=2; i<=8; i++) printf("\t%s",$i); \
                printf("\tGT\t0/1\n"); \
            }' >> out.vcf
            rm -f in.vcf
            bgzip --threads ${N_THREADS} out.vcf
            tabix -f out.vcf.gz
            
            mv out.vcf.gz training_resource.vcf.gz
            mv out.vcf.gz.tbi training_resource.vcf.gz.tbi
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        GetReferenceGaps ~{reference_agp} not_gaps.bed
        
        # Downloading and canonizing the single-sample dipcall VCFs
        touch list.txt
        cat ~{dipcall_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            LocalizeSample ${SAMPLE_ID} ${LINE}
            CanonizeVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz ${SAMPLE_ID}.vcf.gz.tbi ~{min_sv_length} ~{max_sv_length} not_gaps.bed
            echo ${SAMPLE_ID}_canonized.vcf.gz >> list.txt
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done < chunk.csv
        
        # Merging
        Merge list.txt
        
        # Uploading
        gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv training_resource.vcf.'gz*' ~{remote_outdir}/
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
