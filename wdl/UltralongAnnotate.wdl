version 1.0


# 
#
workflow UltralongAnnotate {
    input {
        File chunk_csv
        String remote_outdir
        
        File reference_fai
        File reference_agp
        File training_resource_bed
        File training_resource_vcf_gz
        File training_resource_tbi
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong"
        Int preemptible_number = 4
    }
    parameter_meta {
        chunk_csv: "Format: ID,bai,bam,csi,bcf"
    }
    
    call Impl {
        input:
            chunk_csv = chunk_csv,
            remote_outdir = remote_outdir,
    
            reference_fai = reference_fai,
            reference_agp = reference_agp,
            training_resource_bed = training_resource_bed,
            training_resource_vcf_gz = training_resource_vcf_gz,
            training_resource_tbi = training_resource_tbi,
    
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
        
        File reference_fai
        File reference_agp
        File training_resource_bed
        File training_resource_vcf_gz
        File training_resource_tbi
        
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



        
        # ----------------------- Steps of the pipeline ------------------------
        
        # Builds a BED file that excludes every gap from the AGP file of
        # the reference.
        #
        function GetReferenceGaps() {
            # Computing non-gap regions
            awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    if ( ( $1=="chr1" || $1=="chr2" || $1=="chr3" || $1=="chr4" || $1=="chr5" || $1=="chr6" || $1=="chr7" || $1=="chr8" || $1=="chr9" || $1=="chr10" || \
                           $1=="chr11" || $1=="chr12" || $1=="chr13" || $1=="chr14" || $1=="chr15" || $1=="chr16" || $1=="chr17" || $1=="chr18" || $1=="chr19" || $1=="chr20" || \
                           $1=="chr21" || $1=="chr22" || $1=="chrX" || $1=="chrY" || $1=="chrM" \
                         ) && $5=="N" \
                       ) print $0 \
                 }' ~{reference_agp} > gaps_unsorted.bed
            bedtools sort -i gaps_unsorted.bed -faidx ~{reference_fai} > gaps.bed
            bedtools complement -L -i gaps.bed -g ~{reference_fai} > not_gaps.bed
            
            # Intersecting non-gap regions with the training BED
            bedtools sort -i ~{training_resource_bed} -faidx ~{reference_fai} > training_resource_sorted.bed
            rm -f training_not_gaps_beds.wsv
            ID="0"
            while read ROW; do
                ID=$(( ${ID} + 1 ))
                echo "${ROW}" > ${ID}.bed
                bedtools intersect -a ${ID}.bed -b training_resource_sorted.bed -sorted -g ~{reference_fai} > training_not_gaps_${ID}.bed
                if [ -s training_not_gaps_${ID}.bed ]; then
                    echo "${ID} training_not_gaps_${ID}.bed" >> training_not_gaps_beds.wsv
                else
                    rm -f training_not_gaps_${ID}.bed
                fi
                rm -f ${ID}.bed
            done < not_gaps.bed
            ls -lht *.bed 1>&2
            
            # Removing temporary files
            rm -f gaps_unsorted.bed training_resource_sorted.bed
        }
        
        
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
            
            # Converting to .vcf.gz for sniffles
            bcftools view --threads ${N_THREADS} --output-type z ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}.vcf.gz
            bcftools index --threads ${N_THREADS} -t ${SAMPLE_ID}.vcf.gz
            rm -f ${SAMPLE_ID}.bcf*
        }
        
        
        # Deletes all and only the files downloaded by `LocalizeSample()`.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}*.bam* ${SAMPLE_ID}*.bcf* ${SAMPLE_ID}*.vcf.gz*
        }
        
        
        # Copies truvari's SUPP field from SAMPLE to three tags in INFO. This is
        # necessary, since sniffles might overwrite the SAMPLE column.
        #
        # Remark: the funtion requires an indexed `.vcf.gz` in input, and it
        # outputs an indexed `.vcf.gz` of `bcf`, depending on `OUTPUT_FORMAT`
        # (`z` or `b`).
        #
        function CopySuppToInfo() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local OUTPUT_FORMAT=$3
            local OUTPUT_VCF_GZ=$4
            
            bcftools query --format '%CHROM\t%POS\t%ID\t[%SUPP]\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                printf("%s",$1); \
                for (i=2; i<=NF-1; i++) printf("\t%s",$i); \
                if ($4=="0") printf("\t0\t0\t0");
                else if ($4=="1") printf("\t0\t0\t1");
                else if ($4=="2") printf("\t0\t1\t0");
                else if ($4=="3") printf("\t0\t1\t1");
                else if ($4=="4") printf("\t1\t0\t0");
                else if ($4=="5") printf("\t1\t0\t1");
                else if ($4=="6") printf("\t1\t1\t0");
                else if ($4=="7") printf("\t1\t1\t1");
                printf("\n"); \
            }' | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by pav">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by sniffles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> ${SAMPLE_ID}_header.txt
            # Remark: the order of the callers is now the reverse of the one in
            # which they were bcftools-merged.
            ${TIME_COMMAND} bcftools annotate --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns CHROM,POS,~ID,INFO/SUPP_SNIFFLES,INFO/SUPP_PBSV,INFO/SUPP_PAV --output-type ${OUTPUT_FORMAT} ${INPUT_VCF_GZ} --output ${OUTPUT_VCF_GZ}
            if [ ${OUTPUT_FORMAT} = z ]; then
                bcftools index --threads ${N_THREADS} -f -t ${OUTPUT_VCF_GZ}
            elif [ ${OUTPUT_FORMAT} = b ]; then
                bcftools index --threads ${N_THREADS} -f -c ${OUTPUT_VCF_GZ}
            fi
            
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt ${INPUT_VCF_GZ}*
        }
        
        
        # Remark: we consider sniffles just as a feature annotator, so we keep
        # all recors, even those that are not genotyped as present by sniffles.
        #
        # Remark: the function outputs a `.vcf.gz`, since it's needed by the 
        # following steps.
        #
        function Sniffles() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local ALIGNMENTS_BAM=$3
            
            ${TIME_COMMAND} sniffles --input ${ALIGNMENTS_BAM} --genotype-vcf ${INPUT_VCF} --vcf ${SAMPLE_ID}_out.vcf
            rm -f ${INPUT_VCF} ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # Sorting
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_sniffles.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_sniffles.vcf.gz.tbi
        }
        
        
        cat << 'END' > truvari_bench.sh
#!/bin/bash
SAMPLE_ID=$1
INPUT_VCF_GZ=$2
TRAINING_RESOURCE_VCF_GZ=$3
INFINITY=$4
CHUNK_ID=$5
INCLUDE_BED=$6
${TIME_COMMAND} truvari bench -b ${TRAINING_RESOURCE_VCF_GZ} -c ${INPUT_VCF_GZ} --includebed ${INCLUDE_BED} --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick single -o ${SAMPLE_ID}_truvari_${CHUNK_ID}/
END
        chmod +x truvari_bench.sh
                
        
        # Extracts every record that has a stringent `truvari bench` match with
        # some records in the resource.
        #
        # Remark: we use `--pick single` to force every resource record to be
        # matched with at most one sample record, which is hopefully the
        # most similar to it. This is because we assume that using a
        # contaminated training set in XGBoost downstream is worse than using a
        # slightly smaller training set. With `--pick multi` e.g. two records in
        # the sample VCF might be matched to the same record in the resource 
        # VCF (probably not good) and vice versa (good).
        #
        # Remark: multiple instances of `truvari bench` are run in parallel
        # using `not_gaps.bed`.
        #
        # Remark: in few anecdotal tests, `--pick multi` seems a bit faster than
        # `--pick single` (4m vs 5m with 6 hyperthreading cores).
        #
        # Remark: both the inputs and the output of the function are indexed
        # `.vcf.gz`, since they are needed by `truvari bench`.
        #
        function GetTrainingRecords() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            # Running in parallel
            ${TIME_COMMAND} xargs --arg-file=training_not_gaps_beds.wsv --max-lines=1 --max-procs=${N_THREADS} ./truvari_bench.sh ${SAMPLE_ID} ${INPUT_VCF_GZ} ~{training_resource_vcf_gz} ${INFINITY}
            
            # Concatenating outputs
            rm -f ${SAMPLE_ID}_outputs.txt
            while read ROW; do
                ID=$(echo ${ROW} | cut -d ' ' -f 1)
                echo ${SAMPLE_ID}_truvari_${ID}/tp-comp.vcf.gz >> ${SAMPLE_ID}_outputs.txt
            done < training_not_gaps_beds.wsv
            ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list ${SAMPLE_ID}_outputs.txt --output-type z --output ${SAMPLE_ID}_training.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_training.vcf.gz
            
            # Removing temporary files
            rm -rf ${SAMPLE_ID}_script.sh ${SAMPLE_ID}_outputs.txt ./${SAMPLE_ID}_truvari_*/
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        truvari --help 1>&2
        sniffles --version 1>&2
        
        GetReferenceGaps ~{reference_agp} not_gaps.bed
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            # Skipping the sample if it has already been processed
            TEST=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
        
            # Annotating and marking training records
            LocalizeSample ${SAMPLE_ID} ${LINE}
            CopySuppToInfo ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz z ${SAMPLE_ID}_supp.vcf.gz
            
            echo "Before sniffles:" 1>&2
            bcftools view --no-header ${SAMPLE_ID}_supp.vcf.gz | head -n 5 1>&2
            Sniffles ${SAMPLE_ID} ${SAMPLE_ID}_supp.vcf.gz ${SAMPLE_ID}.bam
            echo "After sniffles:" 1>&2
            bcftools view --header-only ${SAMPLE_ID}_sniffles.vcf.gz 1>&2
            bcftools view --no-header ${SAMPLE_ID}_sniffles.vcf.gz | head -n 5 1>&2
            
            GetTrainingRecords ${SAMPLE_ID} ${SAMPLE_ID}_sniffles.vcf.gz            
            #-------->CopyKanpigFieldsToInfo ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz
        
            # Uploading
            gcloud storage mv ${SAMPLE_ID}_sniffles.vcf.'gz*' ${SAMPLE_ID}_training.vcf.'gz*' ~{remote_outdir}/
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
