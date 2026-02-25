version 1.0


# Studies precision, recall, Mendelian error, de novo rate of the final cohort
# VCF, on the whole genome and for a given set of samples.
#
# Structure of `remote_outdir`:
#
# ├── samples/                for each sample, all the records in the final VCF;
# ├── precision_recall/                 for each sample, precision/recall stats;
# └── mendelian/                         for each sample, mendelian error stats.
#
workflow SV_Integration_PlotPrme {
    input {
        Int max_sv_length = 10000
        
        String remote_workpackage_11_dir
        String remote_outdir
        
        Int precision_recall_bench_method
        Array[String] precision_recall_samples
        Array[String] precision_recall_sex
        Array[File] precision_recall_samples_dipcall_vcf_gz
        Array[File] precision_recall_samples_dipcall_bed
        
        Array[String] mendelian_error_samples
        Array[String] mendelian_error_sex
        File mendelian_error_ped
        Int mendelian_error_n_trios
        
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
        precision_recall_bench_method: "0=truvari bench, 1=vcfdist."
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
            samples = flatten([precision_recall_samples, mendelian_error_samples]),
            
            remote_workpackage_11_dir = remote_workpackage_11_dir,
            remote_outdir = remote_outdir+"/samples",
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    
    # 1. Mendelian error analysis
    scatter (i in range(mendelian_error_n_trios)) {
        call BenchTrio as me_bench_20 {
            input:
                ped_tsv = mendelian_error_ped,
                ped_tsv_row = i+1,
                remote_indir = remote_outdir+"/samples",
                remote_outdir = remote_outdir+"/mendelian",
                min_sv_length = 20,
                max_sv_length = max_sv_length,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                in_flag = [split.out_flag],
                docker_image = docker_image,
                preemptible_number = preemptible_number
        }
        call BenchTrio as me_bench_50 {
            input:
                ped_tsv = mendelian_error_ped,
                ped_tsv_row = i+1,
                remote_indir = remote_outdir+"/samples",
                remote_outdir = remote_outdir+"/mendelian",
                min_sv_length = 50,
                max_sv_length = max_sv_length,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                in_flag = [split.out_flag],
                docker_image = docker_image,
                preemptible_number = preemptible_number
        }
    }
    
    # 2. Precision/recall analysis
    scatter (j in range(length(precision_recall_samples))) {
        call PrecisionRecallAnalysis as pr_analysis_20 {
            input:
                sample_id = precision_recall_samples[j],
                sample_dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[j],
                sample_dipcall_bed = precision_recall_samples_dipcall_bed[j],
                
                remote_outdir = remote_outdir,
            
                bench_method = precision_recall_bench_method,
                min_sv_length = 20,
                max_sv_length = max_sv_length,    
            
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_agp = reference_agp,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                
                in_flag = split.out_flag,
                
                docker_image = docker_image,
                preemptible_number = preemptible_number
        }
        call PrecisionRecallAnalysis as pr_analysis_50 {
            input:
                sample_id = precision_recall_samples[j],
                sample_dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[j],
                sample_dipcall_bed = precision_recall_samples_dipcall_bed[j],
                
                remote_outdir = remote_outdir,
            
                bench_method = precision_recall_bench_method,
                min_sv_length = 50,
                max_sv_length = max_sv_length,    
            
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_agp = reference_agp,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                
                in_flag = split.out_flag,
                
                docker_image = docker_image,
                preemptible_number = preemptible_number
        }
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
        Array[String] samples
        
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
        echo ~{sep="," samples} | tr ',' '\n' > samples.txt
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
        String sample_id
        File sample_dipcall_vcf_gz
        File sample_dipcall_bed
        
        String remote_outdir
        
        Int bench_method
        Int min_sv_length
        Int max_sv_length
        
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
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type b ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.bcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.bcf ${SAMPLE_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_in.bcf
            
            # Removing SNVs, records with unresolved REF/ALT, records that are
            # not marked as present, and records with a FILTER. 
            ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || GT!="alt" || (FILTER!="PASS" && FILTER!=".") || REF="*" || ALT="*"' --output-type b ${SAMPLE_ID}_in.bcf --output ${SAMPLE_ID}_out.bcf
            rm -f ${SAMPLE_ID}_in.bcf* ; mv ${SAMPLE_ID}_out.bcf ${SAMPLE_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_in.bcf
            
            # Removing records in reference gaps
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.bcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.bcf* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the dipcall BED
            ${TIME_COMMAND} bcftools filter --regions-file ~{sample_dipcall_bed} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
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
            local OUTPUT_PREFIX=$3
            
            BED_FLAG_TRUVARI="--includebed ~{sample_dipcall_bed}"
            BED_FLAG_VCFDIST="--bed ~{sample_dipcall_bed}"
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
        
        # Localizing the VCF to benchmark
        gcloud storage cp ~{remote_outdir}/samples/~{sample_id}.'bcf*' .
        
        # Keeping only records in the given length range
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type b ~{sample_id}.bcf --output out.bcf
        rm -f ~{sample_id}.bcf* ; mv out.bcf ~{sample_id}.bcf ; bcftools index --threads ${N_THREADS} -f ~{sample_id}.bcf
        
        # Keeping only records that are genotyped as present. This is important,
        # since truvari bench does not consider GTs when matching.
        ${TIME_COMMAND} bcftools filter --include 'GT="alt"' --output-type z ~{sample_id}.bcf --output out.vcf.gz
        rm -f ~{sample_id}.bcf* ; mv out.vcf.gz ~{sample_id}.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ~{sample_id}.vcf.gz
        
        # Benchmarking
        Benchmark ~{sample_id} ~{sample_id}.vcf.gz final
        
        # Uploading
        gcloud storage cp '*_final_*.txt' ~{remote_outdir}/precision_recall/ 
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


# Remark: `bcftools +trio-dnm2` marks the following triplets (child, father,
# mother) as de novos:
#
# 0|0, 0|0, 1|1
# 0|0, 0|1, 1|1
# 1|1, 1|1, 0|0
# 1|1, 0|0, 1|0
#
# I.e. the disappearance of a record is considered a de novo, and the
# appearance of a record on a haplotype is considered a de novo.
#
# Performance with 8 cores and 16GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# bcftools merge            700%        50M     3s
# bcftools +mendelian2      100%        20M     2s
# bcftools +trio-dnm2       100%        20M     10s
#
task BenchTrio {
    input {
        File ped_tsv
        Int ped_tsv_row
        
        String remote_indir
        String remote_outdir
        
        Int min_sv_length
        Int max_sv_length
        
        File tandem_bed
        File not_tandem_bed
        
        Array[File] in_flag
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 20
        Int preemptible_number
    }
    parameter_meta {
        ped_tsv_row: "The row (one-based) in `ped_tsv` that corresponds to this trio."
        min_sv_length: "The input VCFs are first hard-filtered based on SVLEN, and then fed to the chosen benchmarking tool."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        # @param 1: a VCF that contains only 3 samples (child, parents).
        # 
        function Benchmark() {
            local INPUT_VCF_GZ=$1
            local SAMPLE_ID=$2
            local SUFFIX=$3
            
            # 1. Mendelian error
            ${TIME_COMMAND} bcftools +mendelian2 ${INPUT_VCF_GZ} --ped ped.tsv > ${SAMPLE_ID}_mendelian_${SUFFIX}.txt
            
            # 2. De novo rate
            
            # 2.1 +trio-dnm2
            ${TIME_COMMAND} bcftools +trio-dnm2 --use-NAIVE --chrX GRCh38 --ped ped.tsv --output-type z ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_annotated.vcf.gz
            NUMERATOR=$( bcftools view --no-header --include 'FORMAT/DNM[0]=1' ${SAMPLE_ID}_annotated.vcf.gz | wc -l )
            DENOMINATOR=$( bcftools view --no-header --include 'GT[0]="alt" && COUNT(GT="mis")=0' ${SAMPLE_ID}_annotated.vcf.gz | wc -l )
            echo -e "${NUMERATOR},${DENOMINATOR}" > ${SAMPLE_ID}_dnm1_${SUFFIX}.txt
            rm -f ${SAMPLE_ID}_annotated.vcf.gz*
            
            # 2.2 Simple count (and saving the whole matrix for future analysis)
            ${TIME_COMMAND} truvari anno numneigh --sizemin 1 --refdist 1000 ${INPUT_VCF_GZ} | bgzip > ${SAMPLE_ID}_annotated.vcf.gz
            ${TIME_COMMAND} bcftools query --format '[%GT,][%SQ,][%GQ,][%DP,][%AD,][%KS,]%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors\n' ${SAMPLE_ID}_annotated.vcf.gz > ${SAMPLE_ID}_matrix.txt
            ${TIME_COMMAND} java -cp ~{docker_dir} CountDeNovoSimple ${SAMPLE_ID}_matrix.txt > ${SAMPLE_ID}_dnm2_${SUFFIX}.txt
            rm -f ${SAMPLE_ID}_annotated.vcf.gz*
        }
        
        
        
        
        # --------------------------- Main program -----------------------------
        
        head -n ~{ped_tsv_row} ~{ped_tsv} | tail -n 1 > ped.tsv
        PROBAND_ID=$(cut -f 2 ped.tsv)
        FATHER_ID=$(cut -f 3 ped.tsv)
        MOTHER_ID=$(cut -f 4 ped.tsv)
        
        # Localizing. These files may contain 0/0 records.
        TEST=$(gcloud storage cp ~{remote_indir}/${PROBAND_ID}'*.vcf.gz*' ~{remote_indir}/${FATHER_ID}'*.vcf.gz*' ~{remote_indir}/${MOTHER_ID}'*.vcf.gz*' . && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            gcloud storage cp ~{remote_indir}/${PROBAND_ID}'*.bcf*' ~{remote_indir}/${FATHER_ID}'*.bcf*' ~{remote_indir}/${MOTHER_ID}'*.bcf*' .
        fi
        
        # Ensuring a consistent format
        TEST=$(ls *.vcf.gz && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} bcftools view --output-type z $(ls ${PROBAND_ID}_*.bcf) --output ${PROBAND_ID}_in.vcf.gz &
            ${TIME_COMMAND} bcftools view --output-type z $(ls ${FATHER_ID}_*.bcf) --output ${FATHER_ID}_in.vcf.gz &
            ${TIME_COMMAND} bcftools view --output-type z $(ls ${MOTHER_ID}_*.bcf) --output ${MOTHER_ID}_in.vcf.gz &
            wait
            bcftools index --threads ${N_THREADS} -f -t ${PROBAND_ID}_in.vcf.gz &
            bcftools index --threads ${N_THREADS} -f -t ${FATHER_ID}_in.vcf.gz &
            bcftools index --threads ${N_THREADS} -f -t ${MOTHER_ID}_in.vcf.gz &
            wait
        else
            mv $(ls ${PROBAND_ID}_*.vcf.gz) ${PROBAND_ID}_in.vcf.gz
            mv $(ls ${FATHER_ID}_*.vcf.gz) ${FATHER_ID}_in.vcf.gz
            mv $(ls ${MOTHER_ID}_*.vcf.gz) ${MOTHER_ID}_in.vcf.gz
        fi
        echo ${PROBAND_ID}_in.vcf.gz > list.txt
        echo ${FATHER_ID}_in.vcf.gz >> list.txt
        echo ${MOTHER_ID}_in.vcf.gz >> list.txt
        ls -laht
        
        # Remark: we keep all records from the three files, since otherwise
        # bcftools merge might transform a 0/0 (that would have disappeared)
        # into a ./., affecting the number of records over which Mendelian
        # error is computed.
        
        # Keeping only records in the given length range
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ${PROBAND_ID}_in.vcf.gz --output out.vcf.gz
        rm -f ${PROBAND_ID}_in.vcf.gz* ; mv out.vcf.gz ${PROBAND_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${PROBAND_ID}_in.vcf.gz
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ${FATHER_ID}_in.vcf.gz --output out.vcf.gz
        rm -f ${FATHER_ID}_in.vcf.gz* ; mv out.vcf.gz ${FATHER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${FATHER_ID}_in.vcf.gz
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ${MOTHER_ID}_in.vcf.gz --output out.vcf.gz
        rm -f ${MOTHER_ID}_in.vcf.gz* ; mv out.vcf.gz ${MOTHER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${MOTHER_ID}_in.vcf.gz
        
        # Merging records by ID, since the records in every VCF originate from
        # the same cohort VCF, which had distinct IDs.
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge id --output-type z --file-list list.txt > trio.vcf.gz
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t trio.vcf.gz
        rm -f ${PROBAND_ID}_* ${FATHER_ID}_* ${MOTHER_ID}_*
        ls -laht
        
        # Benchmarking: original merged VCF.
        Benchmark trio.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_all
        
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z trio.vcf.gz --output tr.vcf.gz
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t tr.vcf.gz
        Benchmark tr.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_tr
        
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z trio.vcf.gz --output not_tr.vcf.gz
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t not_tr.vcf.gz
        Benchmark not_tr.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_not_tr
        
        # Benchmarking: VCF with missing->ref.
        ${TIME_COMMAND} bcftools +setGT trio.vcf.gz --output-type z -- --target-gt . --new-gt 0 > trio_no_missing.vcf.gz
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t trio_no_missing.vcf.gz
        Benchmark trio_no_missing.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_all_no_missing
        
        ${TIME_COMMAND} bcftools +setGT tr.vcf.gz --output-type z -- --target-gt . --new-gt 0 > tr_no_missing.vcf.gz
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t tr_no_missing.vcf.gz
        Benchmark tr_no_missing.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_tr_no_missing
        
        ${TIME_COMMAND} bcftools +setGT not_tr.vcf.gz --output-type z -- --target-gt . --new-gt 0 > not_tr_no_missing.vcf.gz
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t not_tr_no_missing.vcf.gz
        Benchmark not_tr_no_missing.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_not_tr_no_missing
        
        rm -f tr*.vcf.gz*
        rm -f not_tr*.vcf.gz*
        
        # Uploading
        gcloud storage cp '*.txt' ~{remote_outdir}/
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
