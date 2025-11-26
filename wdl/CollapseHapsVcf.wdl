version 1.0


# Structure of `remote_output_dir`:
#
# ├── kanpig/            for each sample, its full personalized and re-GT'd VCF;
# ├── precision_recall/                 for each sample, precision/recall stats;
# └── mendelian/                         for each sample, mendelian error stats;
#
workflow CollapseHapsVcf {
    input {
        File kanpig_haps_vcf_gz
        File kanpig_haps_tbi
        
        Float similarity_fraction = 0.9
        String remote_dir
        
        Array[String] precision_recall_samples
        Array[String] precision_recall_sex
        Array[File] precision_recall_bam
        Array[File] precision_recall_bai
        Array[File] precision_recall_samples_dipcall_vcf_gz
        Array[File] precision_recall_samples_dipcall_bed
        
        Array[String] mendelian_error_samples
        Array[String] mendelian_error_sex
        Array[File] mendelian_error_bam
        Array[File] mendelian_error_bai
        File mendelian_error_ped
        Int mendelian_error_n_trios
        
        File reference_fa
        File reference_fai
        File reference_agp
        File standard_chromosomes_bed
        File autosomes_bed
        File tandem_bed
        File ploidy_bed_male
        File ploidy_bed_female
    }
    parameter_meta {
        similarity_fraction: "Applied to both size and sequence similarity."
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    call CollapseHapsVcf {
        input:
            kanpig_haps_vcf_gz = kanpig_haps_vcf_gz,
            kanpig_haps_tbi = kanpig_haps_tbi,
            similarity_fraction = similarity_fraction
    }
    scatter (j in range(length(precision_recall_samples))) {
        call Kanpig as pr_kanpig {
            input:
                sample_id = precision_recall_samples[j],
                sex = precision_recall_sex[j],
                personalized_vcf_gz = CollapseHapsVcf.out_vcf_gz,
                personalized_tbi = CollapseHapsVcf.out_tbi,
                alignments_bam = precision_recall_bam[j],
                alignments_bai = precision_recall_bai[j],
                remote_indir = remote_dir+"/chunks",
                remote_outdir = remote_dir+"/kanpig",
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                ploidy_bed_male = ploidy_bed_male,
                ploidy_bed_female = ploidy_bed_female
        }
        call PrecisionRecallAnalysis as pr_analysis {
            input:
                sample_id = precision_recall_samples[j],
                sample_dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[j],
                sample_dipcall_bed = precision_recall_samples_dipcall_bed[j],
                remote_dir = remote_dir,
            
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_agp = reference_agp,
                standard_chromosomes_bed = standard_chromosomes_bed,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                
                in_flag_kanpig = pr_kanpig.out_flag
        }
    }
    scatter (j in range(length(mendelian_error_samples))) {
        call Kanpig as me_kanpig {
            input:
                sample_id = mendelian_error_samples[j],
                sex = mendelian_error_sex[j],
                personalized_vcf_gz = CollapseHapsVcf.out_vcf_gz,
                personalized_tbi = CollapseHapsVcf.out_tbi,
                alignments_bam = mendelian_error_bam[j],
                alignments_bai = mendelian_error_bai[j],
                remote_indir = remote_dir+"/chunks",
                remote_outdir = remote_dir+"/kanpig",
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                ploidy_bed_male = ploidy_bed_male,
                ploidy_bed_female = ploidy_bed_female
        }
    }
    scatter (j in range(mendelian_error_n_trios)) {
        call BenchTrio as me_bench {
            input:
                ped_tsv = mendelian_error_ped,
                ped_tsv_row = j+1,
                remote_indir = remote_dir+"/kanpig",
                remote_outdir = remote_dir+"/mendelian",
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                autosomes_bed = autosomes_bed,
                in_flag = me_kanpig.out_flag
        }
    }
    
    output {
    }
}




#----------------------- Personalized VCF construction -------------------------

#
task ComplementBed {
    input {
        File tandem_bed
        File reference_fai
        
        Int n_cpu = 1
        Int ram_size_gb = 4
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
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}


#
task CollapseHapsVcf {
    input {
        File kanpig_haps_vcf_gz
        File kanpig_haps_tbi
        
        Float similarity_fraction = 0.9
        
        Int n_cpu = 1
        Int ram_size_gb = 4
    }
    parameter_meta {
        similarity_fraction: "Applied to both size and sequence similarity."
    }
    
    Int disk_size_gb = 10*ceil(size(kanpig_haps_vcf_gz,"GB"))
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
        INFINITY="1000000000"
        
        mv ~{kanpig_haps_vcf_gz} in.vcf.gz
        mv ~{kanpig_haps_tbi} in.vcf.gz.tbi
        
        ${TIME_COMMAND} truvari collapse --input in.vcf.gz --keep maxqual --refdist 0 --pctseq ~{similarity_fraction} --pctsize ~{similarity_fraction} --sizemin 0 --sizemax ${INFINITY} --output out.vcf
        rm -f in.vcf.gz* ; mv out.vcf in.vcf

        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z in.vcf > in.vcf.gz
        rm -f in.vcf ; tabix -f in.vcf.gz
        
        mv in.vcf.gz out.vcf.gz
        mv in.vcf.gz.tbi out.vcf.gz.tbi
    >>>
    
    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999:
#
# TOOL                CPU     RAM     TIME
# kanpig              700%    2.5G    6m
#
task Kanpig {
    input {
        String sample_id
        String sex
        File personalized_vcf_gz
        File personalized_tbi
        File alignments_bam
        File alignments_bai
        
        String remote_indir
        String remote_outdir
        
        String kanpig_params_cohort = "--neighdist 500 --gpenalty 0.04 --hapsim 0.97"
        File reference_fa
        File reference_fai
        File ploidy_bed_male
        File ploidy_bed_female
        
        Int n_cpu = 8
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 4*ceil( size(alignments_bam,"GB") )
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
        export RUST_BACKTRACE="full"
        INFINITY="1000000000"
        
        if [ ~{sex} == "M" ]; then
            PLOIDY_BED=$(echo ~{ploidy_bed_male})
        else
            PLOIDY_BED=$(echo ~{ploidy_bed_female})
        fi
        
        mv ~{personalized_vcf_gz} in.vcf.gz
        mv ~{personalized_tbi} in.vcf.gz.tbi
        
        # Remark: we set --sizemin 10 instead of zero or one, just because
        # kanpig needs --sizemin >= --kmer . The purpose is still to
        # re-genotype every record in the input VCF.
        ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --sizemin 10 --sizemax ${INFINITY} ~{kanpig_params_cohort} --reference ~{reference_fa} --ploidy-bed ${PLOIDY_BED} --input in.vcf.gz --reads ~{alignments_bam} --out out.vcf --sample ~{sample_id}
        rm -f in.vcf.gz* ; mv out.vcf in.vcf
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z in.vcf > ~{sample_id}_kanpig_haps.vcf.gz
        rm -f in.vcf
        tabix -f ~{sample_id}_kanpig_haps.vcf.gz
        N_RECORDS=$(bcftools index --nrecords ~{sample_id}_kanpig_haps.vcf.gz.tbi)
        N_PRESENT_RECORDS=$(bcftools view --no-header --include 'GT="alt"' ~{sample_id}_kanpig_haps.vcf.gz | wc -l)
        echo "${N_RECORDS},${N_PRESENT_RECORDS}" > ~{sample_id}_kanpig_haps_nrecords.txt
        rm -f ~{alignments_bam} ~{reference_fa}
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{sample_id}_kanpig_'*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading ~{sample_id}_kanpig_haps.vcf.gz. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        echo "done" > out.txt
    >>>
    
    output {
        File out_flag = "out.txt"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}




#------------------------------- Benchmarking ----------------------------------

# Performance with 4 cores and 32GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# truvari bench             
# vcfdist                   200%        11G     6m
#
task PrecisionRecallAnalysis {
    input {
        String sample_id
        File sample_dipcall_vcf_gz
        File sample_dipcall_bed
        String chromosome = "chr6"
        
        String remote_dir
        
        Int bench_method = 1
        Int min_sv_length = 20
        Int max_sv_length = 10000
        
        File reference_fa
        File reference_fai
        File reference_agp
        File standard_chromosomes_bed
        File tandem_bed
        File not_tandem_bed
        
        File in_flag_kanpig
        
        Int n_cpu = 4
        Int ram_size_gb = 16
        Int disk_size_gb = 20
    }
    parameter_meta {
        chromosome: "all=do not restrict to a chromosome."        
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
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
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
        
        
        # Puts in canonical form a raw VCF from dipcall. This is identical to
        # `SV_Integration_BuildTrainingResource.wdl`.
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
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Removing SNVs, records that are not marked as present, records
            # with a FILTER, and records with unresolved REF/ALT.
            ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || COUNT(GT="alt")=0 || (FILTER!="PASS" && FILTER!=".") || REF="*" || ALT="*"' --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the standard chromosomes
            ${TIME_COMMAND} bcftools filter --regions-file ${STANDARD_CHROMOSOMES_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Removing records in reference gaps
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the dipcall BED
            ${TIME_COMMAND} bcftools filter --regions-file ~{sample_dipcall_bed} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Making sure SVLEN and SVTYPE are consistently annotated
            truvari anno svinfo --minsize 1 ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the given length range
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_truth.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_truth.vcf.gz.tbi
        }

        
        #
        function Benchmark() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local OUTPUT_PREFIX=$3
                
            if [ ${INPUT_VCF#*.} -eq vcf.gz ]; then
                cp ${INPUT_VCF} ${OUTPUT_PREFIX}_input.vcf.gz
                cp ${INPUT_VCF}.tbi ${OUTPUT_PREFIX}_input.vcf.gz.tbi
            else
                bcftools view --output-type z ${INPUT_VCF} > ${OUTPUT_PREFIX}_input.vcf.gz
                tabix ${OUTPUT_PREFIX}_input.vcf.gz
            fi
            
            # Extracting calls with POS inside and outside TRs
            ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz > ${OUTPUT_PREFIX}_tr.vcf.gz &
            ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz > ${OUTPUT_PREFIX}_not_tr.vcf.gz &
            wait
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_tr.vcf.gz
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_not_tr.vcf.gz
        
            # Benchmarking
            if [ ~{bench_method} -eq 0 ]; then
                rm -rf ./${OUTPUT_PREFIX}_truvari_*
                ${TIME_COMMAND} truvari bench --includebed ~{sample_dipcall_bed} -b ${SAMPLE_ID}_truth.vcf.gz -c ${OUTPUT_PREFIX}_input.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_all/ &
                ${TIME_COMMAND} truvari bench --includebed ~{sample_dipcall_bed} -b ${SAMPLE_ID}_truth_tr.vcf.gz -c ${OUTPUT_PREFIX}_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_tr/ &
                ${TIME_COMMAND} truvari bench --includebed ~{sample_dipcall_bed} -b ${SAMPLE_ID}_truth_not_tr.vcf.gz -c ${OUTPUT_PREFIX}_not_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_not_tr/ &
                wait
                mv ./${OUTPUT_PREFIX}_truvari_all/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_all.txt
                mv ./${OUTPUT_PREFIX}_truvari_tr/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_tr.txt
                mv ./${OUTPUT_PREFIX}_truvari_not_tr/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_not_tr.txt
            else
                # Assumed to be the only job running, so using all cores.
                rm -f ./${OUTPUT_PREFIX}_vcfdist_*
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_input.vcf.gz ${SAMPLE_ID}_truth.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_all/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_tr.vcf.gz ${SAMPLE_ID}_truth_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_tr/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_not_tr.vcf.gz ${SAMPLE_ID}_truth_not_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_not_tr/
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
        tabix -f ~{sample_dipcall_vcf_gz}
        CanonizeDipcallVcf ~{sample_id} ~{sample_dipcall_vcf_gz} ~{sample_dipcall_vcf_gz}.tbi ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
        rm -f ~{sample_dipcall_vcf_gz}
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ~{sample_id}_truth.vcf.gz > ~{sample_id}_truth_tr.vcf.gz &
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ~{sample_id}_truth.vcf.gz > ~{sample_id}_truth_not_tr.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f ~{sample_id}_truth_tr.vcf.gz &
        ${TIME_COMMAND} tabix -f ~{sample_id}_truth_not_tr.vcf.gz &
        wait
        
        # Subsetting the dipcall VCF to the given chromosome, if specified.
        if [ ~{chromosome} != "all" ]; then
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth.vcf.gz > ~{sample_id}_truth_prime.vcf.gz
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth_tr.vcf.gz > ~{sample_id}_truth_tr_prime.vcf.gz
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth_not_tr.vcf.gz > ~{sample_id}_truth_not_tr_prime.vcf.gz
            
            rm -f ~{sample_id}_truth.vcf.gz* ; mv ~{sample_id}_truth_prime.vcf.gz ~{sample_id}_truth.vcf.gz ; bcftools index ~{sample_id}_truth.vcf.gz
            rm -f ~{sample_id}_truth_tr.vcf.gz* ; mv ~{sample_id}_truth_tr_prime.vcf.gz ~{sample_id}_truth_tr.vcf.gz ; bcftools index ~{sample_id}_truth_tr.vcf.gz
            rm -f ~{sample_id}_truth_not_tr.vcf.gz* ; mv ~{sample_id}_truth_not_tr_prime.vcf.gz ~{sample_id}_truth_not_tr.vcf.gz ; bcftools index ~{sample_id}_truth_not_tr.vcf.gz
        fi
        
        # Localizing the re-genotyped haps and records VCFs
        gsutil -m cp ~{remote_dir}/kanpig/~{sample_id}_kanpig_'*.vcf.gz*' .
        
        # Keeping only records that are genotyped as present
        ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type z ~{sample_id}_kanpig_haps.vcf.gz > out.vcf.gz
        rm -f ~{sample_id}_kanpig_haps.vcf.gz* ; mv out.vcf.gz ~{sample_id}_kanpig_haps.vcf.gz ; tabix -f ~{sample_id}_kanpig_haps.vcf.gz
        
        # Benchmarking
        Benchmark ~{sample_id} ~{sample_id}_kanpig_haps.vcf.gz kanpig_haps_~{min_sv_length}bp
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp '*_kanpig_*.txt' ~{remote_dir}/precision_recall/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading kanpig benchmarks. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
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
        
        File tandem_bed
        File not_tandem_bed
        File autosomes_bed
        
        Array[File] in_flag
        
        Int n_cpu = 8
        Int ram_size_gb = 8
        Int disk_size_gb = 20
    }
    parameter_meta {
        ped_tsv_row: "The row (one-based) in `ped_tsv` that corresponds to this trio."
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
            ${TIME_COMMAND} bcftools +trio-dnm2 --use-NAIVE --chrX GRCh38 --ped ped.tsv --output-type z ${INPUT_VCF_GZ} > ${SAMPLE_ID}_annotated.vcf.gz
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
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir}/${PROBAND_ID}_'*.vcf.gz*' ~{remote_indir}/${FATHER_ID}_'*.vcf.gz*' ~{remote_indir}/${MOTHER_ID}_'*.vcf.gz*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading VCF files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        # Restricting to autosomes. This is just for simplicity and is not
        # strictly necessary.
        ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${PROBAND_ID}_kanpig_haps.vcf.gz) > ${PROBAND_ID}_autosomes_haps.vcf.gz &
        ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${FATHER_ID}_kanpig_haps.vcf.gz) > ${FATHER_ID}_autosomes_haps.vcf.gz &
        ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${MOTHER_ID}_kanpig_haps.vcf.gz) > ${MOTHER_ID}_autosomes_haps.vcf.gz &
        wait
        tabix -f ${PROBAND_ID}_autosomes_haps.vcf.gz &
        tabix -f ${FATHER_ID}_autosomes_haps.vcf.gz &
        tabix -f ${MOTHER_ID}_autosomes_haps.vcf.gz &
        wait
        echo ${PROBAND_ID}_autosomes_haps.vcf.gz > list_haps.txt
        echo ${FATHER_ID}_autosomes_haps.vcf.gz >> list_haps.txt
        echo ${MOTHER_ID}_autosomes_haps.vcf.gz >> list_haps.txt
        ls -laht
        
        # Remark: we keep all records from the three files, since otherwise
        # bcftools merge might transform a 0/0 (that would have disappeared)
        # into a ./., affecting the number of records over which Mendelian
        # error is computed.
        
        # Merging records by ID, since the records in every VCF originate from
        # the same cohort VCF, which had distinct IDs.
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge id --output-type z --file-list list_haps.txt > trio_haps.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_haps.vcf.gz
        rm -f ${PROBAND_ID}_* ${FATHER_ID}_* ${MOTHER_ID}_*
        ls -laht
        
        # Benchmarking
        Benchmark trio_haps.vcf.gz ${PROBAND_ID} all_haps
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp '*.txt' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading stats. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
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
