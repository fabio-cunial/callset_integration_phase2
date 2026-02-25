version 1.0


#
workflow InvestigateMaleSamples2 {
    input {
        String chromosome_id
        String remote_dir
        
        Int precision_recall_bench_method
        Int max_sv_length = 10000
        
        Array[String] precision_recall_samples
        Array[File] precision_recall_samples_dipcall_vcf_gz
        Array[File] precision_recall_samples_dipcall_bed
        
        File reference_fa
        File reference_fai
        File reference_agp
        File tandem_bed
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages:latest"
        Int preemptible_number = 4
    }
    parameter_meta {
        precision_recall_bench_method: "0=truvari bench, 1=vcfdist."
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai,
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    scatter (j in range(length(precision_recall_samples))) {
        call PrecisionRecallAnalysis as pr_analysis_20 {
            input:
                sample_id = precision_recall_samples[j],
                sample_dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[j],
                sample_dipcall_bed = precision_recall_samples_dipcall_bed[j],
                chromosome = chromosome_id,
                
                remote_dir = remote_dir,
            
                bench_method = precision_recall_bench_method,
                min_sv_length = 20,
                max_sv_length = max_sv_length,    
            
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_agp = reference_agp,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                
                docker_image = docker_image,
                preemptible_number = preemptible_number
        }
        call PrecisionRecallAnalysis as pr_analysis_50 {
            input:
                sample_id = precision_recall_samples[j],
                sample_dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[j],
                sample_dipcall_bed = precision_recall_samples_dipcall_bed[j],
                chromosome = chromosome_id,
                
                remote_dir = remote_dir,
            
                bench_method = precision_recall_bench_method,
                min_sv_length = 50,
                max_sv_length = max_sv_length,    
            
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_agp = reference_agp,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                
                docker_image = docker_image,
                preemptible_number = preemptible_number
        }
    }
    
    output {
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
        String sample_id
        File sample_dipcall_vcf_gz
        File sample_dipcall_bed
        String chromosome
        
        String remote_dir
        
        Int bench_method
        Int min_sv_length
        Int max_sv_length
        
        File reference_fa
        File reference_fai
        File reference_agp
        File tandem_bed
        File not_tandem_bed
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 4
        Int disk_size_gb = 20
        Int preemptible_number
    }
    parameter_meta {  
        min_sv_length: "The input VCFs (truvari, kanpig and dipcall) are first hard-filtered based on SVLEN, and fed to the chosen benchmarking tool."
        bench_method: "0=truvari bench with default parameters; 1=vcfdist."
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
            local MIN_SV_LENGTH=$4
            local MAX_SV_LENGTH=$5
            local NOT_GAPS_BED=$6
            
            
            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_in.vcf.gz
            mv ${INPUT_TBI} ${SAMPLE_ID}_in.vcf.gz.tbi
            
            # Keeping only records in the given chromosome
            N_RECORDS_BEFORE=$(bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz)
            ${TIME_COMMAND} bcftools view --output-type b ${SAMPLE_ID}_in.vcf.gz ~{chromosome} --output ${SAMPLE_ID}_out.bcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.bcf ${SAMPLE_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_in.bcf
            N_RECORDS_AFTER=$(bcftools index --nrecords ${SAMPLE_ID}_in.bcf)
            
            # Splitting multiallelic records into biallelic records
            N_RECORDS_BEFORE=$(bcftools index --nrecords ${SAMPLE_ID}_in.bcf)
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type b ${SAMPLE_ID}_in.bcf --output ${SAMPLE_ID}_out.bcf
            rm -f ${SAMPLE_ID}_in.bcf* ; mv ${SAMPLE_ID}_out.bcf ${SAMPLE_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_in.bcf
            N_RECORDS_AFTER=$(bcftools index --nrecords ${SAMPLE_ID}_in.bcf)
            
            # Removing SNVs, records with unresolved REF/ALT, records that are
            # not marked as present, and records with a FILTER. 
            # Remark: in chrY we keep calls with any FILTER and any GT,
            # otherwise the number of calls becomes very small.
            N_RECORDS_BEFORE=$(bcftools index --nrecords ${SAMPLE_ID}_in.bcf)
            if [ ~{chromosome} = "chrY" ]; then
                ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1)                                                 || REF="*" || ALT="*"' --output-type b ${SAMPLE_ID}_in.bcf --output ${SAMPLE_ID}_out.bcf
            else
                ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || GT!="alt" || (FILTER!="PASS" && FILTER!=".") || REF="*" || ALT="*"' --output-type b ${SAMPLE_ID}_in.bcf --output ${SAMPLE_ID}_out.bcf
            fi
            rm -f ${SAMPLE_ID}_in.bcf* ; mv ${SAMPLE_ID}_out.bcf ${SAMPLE_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_in.bcf
            N_RECORDS_AFTER=$(bcftools index --nrecords ${SAMPLE_ID}_in.bcf)
            
            # Removing records in reference gaps
            N_RECORDS_BEFORE=$(bcftools index --nrecords ${SAMPLE_ID}_in.bcf)
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.bcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.bcf* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            N_RECORDS_AFTER=$(bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz)
            
            # Keeping only records in the dipcall BED.
            # Remark: we do not do this in chrY, otherwise the number of calls
            # becomes very small.
            N_RECORDS_BEFORE=$(bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz)
            if [ ~{chromosome} != "chrY" ]; then
                ${TIME_COMMAND} bcftools filter --regions-file ~{sample_dipcall_bed} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            fi
            N_RECORDS_AFTER=$(bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz)
            
            # Making sure SVLEN and SVTYPE are consistently annotated
            truvari anno svinfo --minsize 1 ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the given length range
            N_RECORDS_BEFORE=$(bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz)
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            N_RECORDS_AFTER=$(bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz)
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_truth.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_truth.vcf.gz.tbi
        }

        
        #
        function Benchmark() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local OUTPUT_PREFIX=$3
            
            if [ ~{chromosome} == "chrY" ]; then
                BED_FLAG_TRUVARI=" "
                BED_FLAG_VCFDIST=" "
            else
                BED_FLAG_TRUVARI="--includebed ~{sample_dipcall_bed}"
                BED_FLAG_VCFDIST="--bed ~{sample_dipcall_bed}"
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
            if [ ~{bench_method} -eq 0 ]; then
                rm -rf ./${OUTPUT_PREFIX}_truvari_*
                ${TIME_COMMAND} truvari bench ${BED_FLAG_TRUVARI} -b ${SAMPLE_ID}_truth.vcf.gz -c ${OUTPUT_PREFIX}_input.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_all/ &
                ${TIME_COMMAND} truvari bench ${BED_FLAG_TRUVARI} -b ${SAMPLE_ID}_truth_tr.vcf.gz -c ${OUTPUT_PREFIX}_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_tr/ &
                ${TIME_COMMAND} truvari bench ${BED_FLAG_TRUVARI} -b ${SAMPLE_ID}_truth_not_tr.vcf.gz -c ${OUTPUT_PREFIX}_not_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_not_tr/ &
                wait
                mv ./${OUTPUT_PREFIX}_truvari_all/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_all.txt
                mv ./${OUTPUT_PREFIX}_truvari_tr/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_tr.txt
                mv ./${OUTPUT_PREFIX}_truvari_not_tr/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_not_tr.txt
            else
                # Assumed to be the only job running, so using all cores.
                rm -f ./${OUTPUT_PREFIX}_vcfdist_*
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_input.vcf.gz ${SAMPLE_ID}_truth.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} ${BED_FLAG_VCFDIST} --prefix ./${OUTPUT_PREFIX}_vcfdist_all/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_tr.vcf.gz ${SAMPLE_ID}_truth_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} ${BED_FLAG_VCFDIST} --prefix ./${OUTPUT_PREFIX}_vcfdist_tr/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_not_tr.vcf.gz ${SAMPLE_ID}_truth_not_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} ${BED_FLAG_VCFDIST} --prefix ./${OUTPUT_PREFIX}_vcfdist_not_tr/
                mv ./${OUTPUT_PREFIX}_vcfdist_all/precision-recall-summary.tsv ./${SAMPLE_ID}_${OUTPUT_PREFIX}_all.txt
                mv ./${OUTPUT_PREFIX}_vcfdist_tr/precision-recall-summary.tsv ./${SAMPLE_ID}_${OUTPUT_PREFIX}_tr.txt
                mv ./${OUTPUT_PREFIX}_vcfdist_not_tr/precision-recall-summary.tsv ./${SAMPLE_ID}_${OUTPUT_PREFIX}_not_tr.txt
            fi
            
            # Removing temporary files
            rm -f ${OUTPUT_PREFIX}_input.vcf.gz*
        }




        # --------------------------- Main program -----------------------------
        
        FILTER_STRING_TRUVARI="--sizemin ~{min_sv_length} --sizefilt ~{min_sv_length} --sizemax ~{max_sv_length}"
        FILTER_STRING_VCFDIST="--sv-threshold ~{min_sv_length} --largest-variant ~{max_sv_length}"
        # See https://github.com/TimD1/vcfdist/wiki/02-Parameters-and-Usage
        # Remark: `--max-supercluster-size` has to be >= `--largest-variant + 2`
        #         and we set it to 10002 to mimic kanpig's `--sizemax`.
        # Remark: we choose `--cluster gap` since it is faster. We choose 500
        #         to mimic kanpig inter-sample's `--neighdist` (the intra-
        #         sample value would be 1000, which might be too big).
        SV_STRING_VCFDIST="--cluster gap 500 --max-supercluster-size $((~{max_sv_length}+2)) --realign-query --realign-truth"
        
        
        GetReferenceGaps ~{reference_agp} not_gaps.bed
        
        # Canonizing the dipcall VCF
        bcftools index --threads ${N_THREADS} -f -t ~{sample_dipcall_vcf_gz}
        CanonizeDipcallVcf ~{sample_id} ~{sample_dipcall_vcf_gz} ~{sample_dipcall_vcf_gz}.tbi ~{min_sv_length} ~{max_sv_length} not_gaps.bed
        rm -f ~{sample_dipcall_vcf_gz}
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ~{sample_id}_truth.vcf.gz --output ~{sample_id}_truth_tr.vcf.gz &
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ~{sample_id}_truth.vcf.gz --output ~{sample_id}_truth_not_tr.vcf.gz &
        wait
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ~{sample_id}_truth_tr.vcf.gz &
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ~{sample_id}_truth_not_tr.vcf.gz &
        wait
        
        # Localizing the VCFs to benchmark
        gcloud storage cp ~{remote_dir}/~{sample_id}_raw_truvari.vcf.'gz*' .
        gcloud storage cp ~{remote_dir}/~{sample_id}_kanpig.vcf.'gz*' .
        
        # Keeping only records in the given chromosome
        N_RECORDS_BEFORE=$(bcftools index --nrecords ~{sample_id}_raw_truvari.vcf.gz)
        ${TIME_COMMAND} bcftools view --output-type b ~{sample_id}_raw_truvari.vcf.gz ~{chromosome} --output out.bcf
        rm -f ~{sample_id}_raw_truvari.vcf.gz* ; mv out.bcf ~{sample_id}_raw_truvari.bcf ; bcftools index --threads ${N_THREADS} -f ~{sample_id}_raw_truvari.bcf
        N_RECORDS_AFTER=$(bcftools index --nrecords ~{sample_id}_raw_truvari.bcf)
        N_RECORDS_BEFORE=$(bcftools index --nrecords ~{sample_id}_kanpig.vcf.gz)
        ${TIME_COMMAND} bcftools view --output-type b ~{sample_id}_kanpig.vcf.gz ~{chromosome} --output out.bcf
        rm -f ~{sample_id}_kanpig.vcf.gz* ; mv out.bcf ~{sample_id}_kanpig.bcf ; bcftools index --threads ${N_THREADS} -f ~{sample_id}_kanpig.bcf
        N_RECORDS_AFTER=$(bcftools index --nrecords ~{sample_id}_kanpig.bcf)
        
        # Keeping only records in the given length range
        N_RECORDS_BEFORE=$(bcftools index --nrecords ~{sample_id}_raw_truvari.bcf)
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type b ~{sample_id}_raw_truvari.bcf --output out.bcf
        rm -f ~{sample_id}_raw_truvari.bcf* ; mv out.bcf ~{sample_id}_truvari.bcf ; bcftools index --threads ${N_THREADS} -f ~{sample_id}_truvari.bcf
        N_RECORDS_AFTER=$(bcftools index --nrecords ~{sample_id}_truvari.bcf)
        N_RECORDS_BEFORE=$(bcftools index --nrecords ~{sample_id}_kanpig.bcf)
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ~{sample_id}_kanpig.bcf --output out.bcf
        rm -f ~{sample_id}_kanpig.bcf* ; mv out.bcf ~{sample_id}_kanpig.bcf ; bcftools index --threads ${N_THREADS} -f ~{sample_id}_kanpig.bcf
        N_RECORDS_AFTER=$(bcftools index --nrecords ~{sample_id}_kanpig.bcf)
        
        # Keeping only records that are genotyped as present. This is important,
        # since truvari bench does not consider GTs when matching.
        bcftools view --output-type z ~{sample_id}_truvari.bcf --output ~{sample_id}_truvari_all.vcf.gz
        bcftools index --threads ${N_THREADS} -f -t ~{sample_id}_truvari_all.vcf.gz
        
        N_RECORDS_BEFORE=$(bcftools index --nrecords ~{sample_id}_truvari.bcf)
        ${TIME_COMMAND} bcftools filter --include 'GT="alt"' --output-type z ~{sample_id}_truvari.bcf --output out.vcf.gz
        rm -f ~{sample_id}_truvari.bcf* ; mv out.vcf.gz ~{sample_id}_truvari_present.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ~{sample_id}_truvari_present.vcf.gz
        N_RECORDS_AFTER=$(bcftools index --nrecords ~{sample_id}_truvari_present.vcf.gz)
        
        N_RECORDS_BEFORE=$(bcftools index --nrecords ~{sample_id}_kanpig.bcf)
        ${TIME_COMMAND} bcftools filter --include 'GT="alt"' --output-type z ~{sample_id}_kanpig.bcf --output out.vcf.gz
        rm -f ~{sample_id}_kanpig.bcf* ; mv out.vcf.gz ~{sample_id}_kanpig.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ~{sample_id}_kanpig.vcf.gz
        N_RECORDS_AFTER=$(bcftools index --nrecords ~{sample_id}_kanpig.vcf.gz)
        
        # Benchmarking
        Benchmark ~{sample_id} ~{sample_id}_truvari_all.vcf.gz truvari_all_~{min_sv_length}bp
        Benchmark ~{sample_id} ~{sample_id}_truvari_present.vcf.gz truvari_present_~{min_sv_length}bp
        Benchmark ~{sample_id} ~{sample_id}_kanpig.vcf.gz kanpig_~{min_sv_length}bp
        
        # Uploading
        gcloud storage cp '*_truvari_*.txt' ~{remote_dir}/precision_recall/~{chromosome}/
        gcloud storage cp '*_kanpig_*.txt' ~{remote_dir}/precision_recall/~{chromosome}/
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
