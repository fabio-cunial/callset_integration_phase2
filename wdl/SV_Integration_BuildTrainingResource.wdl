version 1.0


# Builds a resource of high-confidence true calls for XGBoost training, by
# canonizing and trivially merging multiple dipcall VCFs.
#
workflow SV_Integration_BuildTrainingResource {
    input {
        File dipcall_tsv
        
        Int min_sv_length = 20
        Int max_sv_length = 10000
        
        File reference_fai
        File standard_chromosomes_bed
        File reference_agp
        
        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        dipcall_tsv: "Format of each row: ID, DIPCALL_BED, DIPCALL_VCF. We assume that every VCF is sorted and sequence-resolved."
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            dipcall_tsv = dipcall_tsv,
            
            min_sv_length = min_sv_length,
            max_sv_length = max_sv_length,
            
            reference_fai = reference_fai,
            standard_chromosomes_bed = standard_chromosomes_bed,
            reference_agp = reference_agp,
            
            remote_outdir = remote_outdir,
            
            docker_image = docker_image
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
        File standard_chromosomes_bed
        File reference_agp
        
        String remote_outdir
        
        String docker_image
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
        GSUTIL_DELAY_S="600"
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            DIPCALL_BED=$(echo ${LINE} | cut -d , -f 2)
            DIPCALL_VCF_GZ=$(echo ${LINE} | cut -d , -f 3)
            while : ; do
                TEST=$(gcloud storage cp ${DIPCALL_BED} ./${SAMPLE_ID}.bed && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${DIPCALL_BED}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gcloud storage cp ${DIPCALL_VCF_GZ} ./${SAMPLE_ID}.vcf.gz && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${DIPCALL_VCF_GZ}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            tabix -@ ${N_THREADS} -f ${SAMPLE_ID}.vcf.gz
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
        # Remark: we only keep INS and DEL records, since dipcall's replacement
        # records have no SVTYPE, they would be assigned SVTYPE=SUB by our
        # annotation script, and they would not be matched to any record in an
        # intra-sample VCF since `truvari bench` uses type as a matching
        # criterion.
        #
        # Alternatively, we could keep all records and run `truvari bench`
        # without considering SVTYPE, but dipcall's replacement records can be
        # strange, e.g.:
        #
        # REF: CAAAAAAAAAAAAAAAAAAAAAA
        # ALT: C
        # ALT: CAAAAAAAAAAAAAAAAAAAAAAAAAA
        #
        # REF: CGGTGGTCCTCCTTGCCGGTGGTCCTCCTTCCTGGTGGTTCTCCTTCCTGGTGGTCCTCCTTCCT
        # ALT: C
        # ALT: TGGTGGTCCTCCTTGCCGGTGGTCCTCCTTCCTGGTGGTTCTCCTTCCTGGTGGTCCTCCTTCCT
        #
        # This means that e.g. INVs are never part of the training resource:
        # this is probably fine, since the resource does not need to be
        # comprehensive.
        #
        function CanonizeDipcallVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local INPUT_TBI=$3
            local MIN_SV_LENGTH=$4
            local MAX_SV_LENGTH=$5
            local STANDARD_CHROMOSOMES_BED=$6
            local NOT_GAPS_BED=$7
            
            
            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_in.vcf.gz
            mv ${INPUT_TBI} ${SAMPLE_ID}_in.vcf.gz.tbi
            
            # Keeping only records in the standard chromosomes
            ${TIME_COMMAND} bcftools filter --regions-file ${STANDARD_CHROMOSOMES_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Removing records in reference gaps
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the dipcall BED
            ${TIME_COMMAND} bcftools filter --regions-file ${SAMPLE_ID}.bed --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics -any --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Removing SNVs, replacement records, records that are not marked
            # as ALT, records with a FILTER, and records with unresolved
            # REF/ALT.
            ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || (STRLEN(REF)>1 && STRLEN(ALT)>1) || GT!="alt" || (FILTER!="PASS" && FILTER!=".") || REF="*" || ALT="*"' --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Making sure SVLEN and SVTYPE are consistently annotated        
            ${TIME_COMMAND} java -cp ~{docker_dir} AddSvtypeSvlen ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only INS and DEL in the given length range
            ${TIME_COMMAND} bcftools filter --include '(SVTYPE="INS" || SVTYPE="DEL") && ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_canonized.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_canonized.vcf.gz.tbi
        }
        
        
        function CanonizeThread() {
            local THREAD_ID=$1
            local CHUNK_CSV=$2
            
            while read LINE; do
                SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
                LocalizeSample ${SAMPLE_ID} ${LINE}
                CanonizeDipcallVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz ${SAMPLE_ID}.vcf.gz.tbi ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
                echo ${SAMPLE_ID}_canonized.vcf.gz >> ${THREAD_ID}_list.txt
                DelocalizeSample ${SAMPLE_ID}
            done < ${CHUNK_CSV}
        }
        
        
        # A simple `bcftools merge` without any additional collapse.
        #
        function Merge() {
            local FILE_LIST_TXT=$1
            
            # Merging
            date
            bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list ${FILE_LIST_TXT} --output-type v | cut -f 1-10 > out.vcf
            date
            mv out.vcf in.vcf
            
            # Setting sample name
            echo "MERGED" > samples.txt
            bcftools reheader --samples samples.txt in.vcf > out.vcf
            mv out.vcf in.vcf
            
            # Enforcing artificial GTs on the single sample column
            ${TIME_COMMAND} bcftools +setGT --output-type z in.vcf -- --target-gt a --new-gt c:0/1 > out.vcf.gz
            rm -f in.vcf ; mv out.vcf.gz in.vcf.gz ; tabix -@ ${N_THREADS} in.vcf.gz
            
            mv in.vcf.gz training_resource.vcf.gz
            mv in.vcf.gz.tbi training_resource.vcf.gz.tbi
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        GetReferenceGaps ~{reference_agp} not_gaps.bed
        
        # Downloading and canonizing the single-sample dipcall VCFs in parallel
        touch list.txt
        cat ~{dipcall_tsv} | tr '\t' ',' > samples.csv
        N_ROWS=$(wc -l < samples.csv)
        if [ ${N_ROWS} -gt ${N_THREADS} ]; then
            N_ROWS_PER_THREAD=$(( ${N_ROWS} / ${N_THREADS} ))
            split -l ${N_ROWS_PER_THREAD} -d -a 4 samples.csv chunk_
        else
            mv samples.csv chunk_0
        fi
        for FILE in $(ls chunk_*); do
            CanonizeThread ${FILE#chunk_*} ${FILE} &
        done
        wait
        cat *_list.txt > list.txt
        
        # Merging
        Merge list.txt
        
        # Uploading
        gcloud storage mv training_resource.vcf.'gz*' ~{remote_outdir}/
    >>>

    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
