version 1.0


# --------->Studies precision, recall, Mendelian error, de novo rate of the final cohort
# VCF, on the whole genome and for a given set of samples.
#
# Structure of `remote_outdir`:
#
# ├── samples/                for each sample, all the records in the final VCF;
# ├── precision_recall/                 for each sample, precision/recall stats;
# └── mendelian/                         for each sample, mendelian error stats.
#
workflow SV_Integration_UltralongAnalysis {
    input {
        String remote_workpackage_11_dir
        String remote_outdir
        
        String precision_recall_samples_csv
        String truvari_bench_flags
        Int use_dipcall_bed
        
        File reference_fa
        File reference_fai
        File reference_agp
        File tandem_bed
        File ploidy_bed_male
        File ploidy_bed_female
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages:latest"
        Int preemptible_number = 4
    }
    parameter_meta {
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai,
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    
    
    call SplitBcfBySample as split {
        input:
            samples_csv = precision_recall_samples_csv,
            
            remote_workpackage_11_dir = remote_workpackage_11_dir,
            remote_outdir = remote_outdir+"/samples",
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    
    # 1. Precision/recall analysis
    call PrecisionRecallAnalysis {
        input:
            samples_csv = precision_recall_samples_csv,
            remote_outdir = remote_outdir,
        
            truvari_bench_flags = truvari_bench_flags,
            use_dipcall_bed = use_dipcall_bed,
        
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            reference_agp = reference_agp,
            tandem_bed = ComplementBed.sorted_bed,
            not_tandem_bed = ComplementBed.complement_bed,
            
            in_flag = split.out_flag,
            
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    
    output {
    }
}


# Writes to a separate file every sample column.
#
# Remark: we keep every record, not just those genotyped as present, to support
# all types of analysis downstream.
#
# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, SSD:
#
# TOOL                              CPU     RAM     TIME
# bcftools +split                   100%    500M    3m
#
task SplitBcfBySample {
    input {
        File samples_csv
        
        String remote_workpackage_11_dir
        String remote_outdir
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 50
        Int preemptible_number
    }
    parameter_meta {
        remote_outdir: "The result of the split is stored in this bucket location."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        
        gcloud storage cp ~{remote_workpackage_11_dir}/merged.'bcf*' .
        cut -d , -f 1 ~{samples_csv} | sort | uniq > samples.txt
        ${TIME_COMMAND} bcftools +split --samples-file samples.txt --output-type b --output . merged.bcf
        rm -f merged.bcf*
        for FILE in $(ls *.bcf); do
            bcftools index --threads ${N_THREADS} -f ${FILE}
        done
        ls -laht
        gcloud storage cp '*.bcf*' ~{remote_outdir}/        
        
        # Fake output
        echo "done" > out.txt
    >>>
    
    output {
        File out_flag = "out.txt"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}


#
task ComplementBed {
    input {
        File tandem_bed
        File reference_fai
        
        String docker_image
        Int n_cpu = 1
        Int ram_size_gb = 4
        Int preemptible_number
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(tandem_bed,"GB"))
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"

        
        ${TIME_COMMAND} bedtools sort -i ~{tandem_bed} -faidx ~{reference_fai} > sorted.bed
        ${TIME_COMMAND} bedtools complement -i sorted.bed -L -g ~{reference_fai} > complement.bed
    >>>
    
    output {
        File sorted_bed = "sorted.bed"
        File complement_bed = "complement.bed"
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


# Performance with 4 cores and 32GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# truvari bench             
# vcfdist
#
task PrecisionRecallAnalysis {
    input {
        File samples_csv
        String remote_outdir
        
        Int min_sv_length = 10000
        Int max_sv_length = 999999999
        String truvari_bench_flags
        Int use_dipcall_bed
        
        File reference_fa
        File reference_fai
        File reference_agp
        File tandem_bed
        File not_tandem_bed
        
        File in_flag
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 4
        Int disk_size_gb = 20
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
        
        
        # Puts in canonical form a raw VCF from dipcall. This is similar to
        # `SV_Integration_BuildTrainingResource.wdl`.
        #
        function CanonizeDipcallVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local INPUT_TBI=$3
            local DIPCALL_BED=$4
            local MIN_SV_LENGTH=$5
            local MAX_SV_LENGTH=$6
            local NOT_GAPS_BED=$7
            
            
            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_in.vcf.gz
            mv ${INPUT_TBI} ${SAMPLE_ID}_in.vcf.gz.tbi
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type b ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.bcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.bcf ${SAMPLE_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_in.bcf
            
            # Removing SNVs, records with unresolved REF/ALT, records that are
            # not marked as present, and records with a FILTER. 
            ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || (GT!="alt" && GT!=".|1" && GT!="1|." && GT!="./1" && GT!="1/.") || (FILTER!="PASS" && FILTER!=".") || REF="*" || ALT="*"' --output-type b ${SAMPLE_ID}_in.bcf --output ${SAMPLE_ID}_out.bcf
            rm -f ${SAMPLE_ID}_in.bcf* ; mv ${SAMPLE_ID}_out.bcf ${SAMPLE_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_in.bcf
            
            # Removing records in reference gaps
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.bcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.bcf* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the dipcall BED, if needed.
            if [ ~{use_dipcall_bed} -eq 1 ]; then
                ${TIME_COMMAND} bcftools filter --regions-file ${DIPCALL_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            fi
            
            # Making sure SVLEN and SVTYPE are consistently annotated
            truvari anno svinfo --minsize 1 ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the given length range
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_truth.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_truth.vcf.gz.tbi
        }

        
        #
        function Benchmark() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local DIPCALL_BED=$3
            local OUTPUT_PREFIX=$4
            
            if [ ~{use_dipcall_bed} -eq 1 ]; then
                BED_FLAG_TRUVARI="--includebed ${DIPCALL_BED}"
                BED_FLAG_VCFDIST="--bed ${DIPCALL_BED}"
            else
                BED_FLAG_TRUVARI=" "
                BED_FLAG_VCFDIST=" "
            fi
            if [ ${INPUT_VCF#*.} = vcf.gz ]; then
                cp ${INPUT_VCF} ${OUTPUT_PREFIX}_input.vcf.gz
                cp ${INPUT_VCF}.tbi ${OUTPUT_PREFIX}_input.vcf.gz.tbi
            else
                bcftools view --output-type z ${INPUT_VCF} --output ${OUTPUT_PREFIX}_input.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${OUTPUT_PREFIX}_input.vcf.gz
            fi
            
            # Extracting calls with POS inside and outside TRs
            ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz --output ${OUTPUT_PREFIX}_tr.vcf.gz &
            ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz --output ${OUTPUT_PREFIX}_not_tr.vcf.gz &
            wait
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${OUTPUT_PREFIX}_tr.vcf.gz &
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${OUTPUT_PREFIX}_not_tr.vcf.gz &
            wait
        
            # Benchmarking
            rm -rf ./${OUTPUT_PREFIX}_truvari_*
            ${TIME_COMMAND} truvari bench ${BED_FLAG_TRUVARI} -b ${SAMPLE_ID}_truth.vcf.gz        -c ${OUTPUT_PREFIX}_input.vcf.gz  --sizemin ~{min_sv_length} --sizefilt ~{min_sv_length} --sizemax ~{max_sv_length} ~{truvari_bench_flags} -o ./${OUTPUT_PREFIX}_truvari_all/ &
            ${TIME_COMMAND} truvari bench ${BED_FLAG_TRUVARI} -b ${SAMPLE_ID}_truth_tr.vcf.gz     -c ${OUTPUT_PREFIX}_tr.vcf.gz     --sizemin ~{min_sv_length} --sizefilt ~{min_sv_length} --sizemax ~{max_sv_length} ~{truvari_bench_flags} -o ./${OUTPUT_PREFIX}_truvari_tr/ &
            ${TIME_COMMAND} truvari bench ${BED_FLAG_TRUVARI} -b ${SAMPLE_ID}_truth_not_tr.vcf.gz -c ${OUTPUT_PREFIX}_not_tr.vcf.gz --sizemin ~{min_sv_length} --sizefilt ~{min_sv_length} --sizemax ~{max_sv_length} ~{truvari_bench_flags} -o ./${OUTPUT_PREFIX}_truvari_not_tr/ &
            wait
            mv ./${OUTPUT_PREFIX}_truvari_all/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_all.txt
            mv ./${OUTPUT_PREFIX}_truvari_tr/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_tr.txt
            mv ./${OUTPUT_PREFIX}_truvari_not_tr/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_not_tr.txt
            
            # Removing temporary files
            rm -f ${OUTPUT_PREFIX}_input.vcf.gz*
        }




        # --------------------------- Main program -----------------------------
        
        GetReferenceGaps ~{reference_agp} not_gaps.bed
        
        while read ROW; do
            SAMPLE_ID=$(echo ${ROW} | cut -d , -f 1)
            DIPCALL_BED_URI=$(echo ${ROW} | cut -d , -f 2)
            DIPCALL_VCF_URI=$(echo ${ROW} | cut -d , -f 3)
            
            # Downloading
            gcloud storage cp ${DIPCALL_BED_URI} ${SAMPLE_ID}_dipcall.bed
            gcloud storage cp ${DIPCALL_VCF_URI} ${SAMPLE_ID}_dipcall.vcf.gz
            
            # Canonizing the dipcall VCF
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_dipcall.vcf.gz
            CanonizeDipcallVcf ${SAMPLE_ID} ${SAMPLE_ID}_dipcall.vcf.gz ${SAMPLE_ID}_dipcall.vcf.gz.tbi ${SAMPLE_ID}_dipcall.bed ~{min_sv_length} ~{max_sv_length} not_gaps.bed
            rm -f ${SAMPLE_ID}_dipcall.vcf.gz*
            ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${SAMPLE_ID}_truth.vcf.gz --output ${SAMPLE_ID}_truth_tr.vcf.gz &
            ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${SAMPLE_ID}_truth.vcf.gz --output ${SAMPLE_ID}_truth_not_tr.vcf.gz &
            wait
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_truth_tr.vcf.gz &
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_truth_not_tr.vcf.gz &
            wait
        
            # Localizing the VCF to benchmark
            gcloud storage cp ~{remote_outdir}/samples/${SAMPLE_ID}.'bcf*' .
        
            # Keeping only records in the given length range
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type b ${SAMPLE_ID}.bcf --output out.bcf
            rm -f ${SAMPLE_ID}.bcf* ; mv out.bcf ${SAMPLE_ID}.bcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}.bcf
        
            # Keeping only records that are genotyped as present. This is
            # important, since truvari bench does not consider GTs when
            # matching.
            ${TIME_COMMAND} bcftools filter --include 'GT="alt" || GT=".|1" || GT="1|." || GT="./1" || GT="1/."' --output-type z ${SAMPLE_ID}.bcf --output out.vcf.gz
            rm -f ${SAMPLE_ID}.bcf* ; mv out.vcf.gz ${SAMPLE_ID}.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}.vcf.gz
        
            # Benchmarking
            Benchmark ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz ${SAMPLE_ID}_dipcall.bed ultralong
        
            # Uploading
            gcloud storage cp '*.txt' ~{remote_outdir}/precision_recall/ 
        done < ~{samples_csv}
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
