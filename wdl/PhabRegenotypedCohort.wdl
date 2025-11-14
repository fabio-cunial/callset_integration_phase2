version 1.0


# 
#
workflow PhabRegenotypedCohort {
    input {
        File cohort_regenotyped_vcf_gz
        File cohort_regenotyped_tbi
        String limit_to_chromosome = "chr6"
        
        Int min_sv_length = 20
        Int max_sv_length = 10000
        Int n_phab_tasks = 1
        String phab_align = "wfa"
        String remote_output_dir

        File reference_fa
        File reference_fai
    }
    parameter_meta {
        cohort_regenotyped_vcf_gz: "The result of re-genotyping a truvari-collapsed VCF with kanpig."
        phab_align: "Possible values: mafft, wfa, poa."
        min_sv_length: "Keep only calls >=X from the input regenotyped cohort, and from the output of truvari phab."
    }
    
    call PrepareCohortVcf {
        input:
            cohort_regenotyped_vcf_gz = cohort_regenotyped_vcf_gz,
            cohort_regenotyped_tbi = cohort_regenotyped_tbi,
            limit_to_chromosome = limit_to_chromosome,
            min_sv_length = min_sv_length,
            max_sv_length = max_sv_length
    }
    call GetRegions {
        input:
            cohort_regenotyped_vcf_gz = PrepareCohortVcf.out_vcf_gz,
            cohort_regenotyped_tbi = PrepareCohortVcf.out_tbi,
            n_phab_tasks = n_phab_tasks
    }
    scatter (chunk_bed in GetRegions.region_chunks) {
        call PhabChunk {
            input:
                regions_bed = chunk_bed,
                cohort_regenotyped_vcf_gz = PrepareCohortVcf.out_vcf_gz,
                cohort_regenotyped_tbi = PrepareCohortVcf.out_tbi,
                remote_output_dir = remote_output_dir,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                phab_align = phab_align,
                min_sv_length = min_sv_length
        }
    }
    
    output {
    }
}


# Performance on a VM with 8 virtual cores and 8GB of RAM, chr6:
#
# TOOL              CPU%    RAM     TIME
# bcftools view     300%    300M    20m
#
task PrepareCohortVcf {
    input {
        File cohort_regenotyped_vcf_gz
        File cohort_regenotyped_tbi
        
        String limit_to_chromosome = "chr6"
        Int min_sv_length
        Int max_sv_length
        
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
        cohort_regenotyped_vcf_gz: "The result of re-genotyping a truvari-collapsed VCF with kanpig."
    }
    
    Int disk_size_gb = 4*ceil(size(cohort_regenotyped_vcf_gz,"GB"))
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        INFINITY="1000000000"
        
        mv ~{cohort_regenotyped_vcf_gz} in.vcf.gz
        mv ~{cohort_regenotyped_tbi} in.vcf.gz.tbi
        
        # Subsetting to chromosome
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --regions ~{limit_to_chromosome} --output-type z in.vcf.gz > out.vcf.gz
        rm -f in.vcf.gz* ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        # Preparing for delocalization
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
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Performance on a VM with 8 virtual cores and 8GB of RAM:
#
# TOOL              CPU%    RAM     TIME
# GetKanpigRegions  100%    400M    10m
#
task GetRegions {
    input {
        File cohort_regenotyped_vcf_gz
        File cohort_regenotyped_tbi

        Int n_phab_tasks
        
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
        cohort_regenotyped_vcf_gz: "The result of re-genotyping a truvari-collapsed VCF with kanpig."
    }
    
    Int disk_size_gb = 4*ceil(size(cohort_regenotyped_vcf_gz,"GB"))
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        INFINITY="1000000000"
        
        mv ~{cohort_regenotyped_vcf_gz} in.vcf.gz
        mv ~{cohort_regenotyped_tbi} in.vcf.gz.tbi
        
        # Computing regions, using WLOG the first sample. Each sample might have
        # different region IDs, since kanpig was run in isolation on each
        # sample, but regions should be identical across all samples.
        ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigRegions in.vcf.gz ${INFINITY} > regions.bed
        N_ROWS=$(wc -l < regions.bed)
        if [ ~{n_phab_tasks} -gt 1 -a ${N_ROWS} -gt ~{n_phab_tasks} ]; then
            N_ROWS_PER_TASK=$(( ${N_ROWS} / ~{n_phab_tasks} ))
            split -l ${N_ROWS_PER_TASK} -d -a 4 regions.bed chunk_ && echo 0 || echo "Error"
        else
            mv regions.bed chunk_0000
        fi
        ls -laht
    >>>
    
    output {
        Array[File] region_chunks = glob("chunk_*")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Remark: the output of `truvari phab` can include also short calls and does
# not preserve any input field, e.g.:
#
# chr1    996401  .       A       C       .       .       .       GT      1/0
# chr1    996410  .       T       G       .       .       .       GT      1/0
# chr1    996414  .       G       A       .       .       .       GT      1/0
# chr1    996427  .       G       A       .       .       .       GT      1/0
# chr1    996458  .       A       G       .       .       .       GT      1/0
#
# Performance on a VM with 8 virtual cores and 8GB of RAM, one 90bp window:
#
# METHOD    CPU%    RAM     TIME
# wfa       
# poa       800%    150M    15m
# mafft     
#
task PhabChunk {
    input {
        File regions_bed
        File cohort_regenotyped_vcf_gz
        File cohort_regenotyped_tbi

        File reference_fa
        File reference_fai
        
        String phab_align
        Int min_sv_length
        
        String remote_output_dir
        
        Int n_cpu = 8
        Int ram_size_gb = 8
    }
    parameter_meta {
        remote_output_dir: "Without final slash"
    }
    
    Int disk_size_gb = 4*ceil(size(cohort_regenotyped_vcf_gz,"GB"))
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
        
        CHUNK_ID=$(basename ~{regions_bed})
        CHUNK_ID=${CHUNK_ID#chunk_*}
        
        mv ~{cohort_regenotyped_vcf_gz} in.vcf.gz
        mv ~{cohort_regenotyped_tbi} in.vcf.gz.tbi
        
        head -n 1 ~{regions_bed} > regions.bed
        
        ${TIME_COMMAND} truvari phab --debug -t ${N_THREADS} --align ~{phab_align} -b in.vcf.gz -r regions.bed -f ~{reference_fa} -o out.vcf.gz
        rm -f in.vcf.gz* ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        ${TIME_COMMAND} truvari anno svinfo -o out.vcf.gz in.vcf.gz
        rm -f in.vcf.gz* ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length} --output-type z in.vcf.gz > out.vcf.gz
        rm -f in.vcf.gz* ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        mv in.vcf.gz chunk_${CHUNK_ID}.vcf.gz
        mv in.vcf.gz.tbi chunk_${CHUNK_ID}.vcf.gz.tbi
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp chunk_${CHUNK_ID}.vcf.'gz*' ~{remote_output_dir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading VCF. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
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




#------------------------------- Benchmarking ----------------------------------

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
        String min_n_samples
        
        Int bench_method
        Int min_sv_length
        Int max_sv_length
        
        File reference_fa
        File reference_fai
        File reference_agp
        File standard_chromosomes_bed
        File tandem_bed
        File not_tandem_bed
        
        File in_flag_truvari
        File in_flag_kanpig
        
        Int n_cpu = 2
        Int ram_size_gb = 4
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
        
        # Localizing the VCFs to benchmark
        gsutil -m cp ~{remote_dir}/truvari/~{sample_id}_truvari.'bcf*' .
        gsutil -m cp ~{remote_dir}/~{min_n_samples}_samples/kanpig/~{sample_id}_kanpig.vcf.'gz*' .
        
        # Subsetting every VCF to the given chromosome, if specified.
        if [ ~{chromosome} != "all" ]; then
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth.vcf.gz > ~{sample_id}_truth_prime.vcf.gz
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth_tr.vcf.gz > ~{sample_id}_truth_tr_prime.vcf.gz
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth_not_tr.vcf.gz > ~{sample_id}_truth_not_tr_prime.vcf.gz
            
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type b ~{sample_id}_truvari.bcf > ~{sample_id}_truvari_prime.bcf
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_kanpig.vcf.gz > ~{sample_id}_kanpig_prime.vcf.gz
            
            rm -f ~{sample_id}_truth.vcf.gz* ; mv ~{sample_id}_truth_prime.vcf.gz ~{sample_id}_truth.vcf.gz ; bcftools index ~{sample_id}_truth.vcf.gz
            rm -f ~{sample_id}_truth_tr.vcf.gz* ; mv ~{sample_id}_truth_tr_prime.vcf.gz ~{sample_id}_truth_tr.vcf.gz ; bcftools index ~{sample_id}_truth_tr.vcf.gz
            rm -f ~{sample_id}_truth_not_tr.vcf.gz* ; mv ~{sample_id}_truth_not_tr_prime.vcf.gz ~{sample_id}_truth_not_tr.vcf.gz ; bcftools index ~{sample_id}_truth_not_tr.vcf.gz
            
            rm -f ~{sample_id}_truvari.bcf* ; mv ~{sample_id}_truvari_prime.bcf ~{sample_id}_truvari.bcf ; bcftools index ~{sample_id}_truvari.bcf
            rm -f ~{sample_id}_kanpig.vcf.gz* ; mv ~{sample_id}_kanpig_prime.vcf.gz ~{sample_id}_kanpig.vcf.gz ; bcftools index ~{sample_id}_kanpig.vcf.gz
        fi
        
        # Keeping only records in the given length range
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ~{sample_id}_truvari.bcf > out.vcf.gz
        rm -f ~{sample_id}_truvari.bcf* ; mv out.vcf.gz ~{sample_id}_truvari.vcf.gz ; tabix -f ~{sample_id}_truvari.vcf.gz
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ~{sample_id}_kanpig.vcf.gz > out.vcf.gz
        rm -f ~{sample_id}_kanpig.vcf.gz* ; mv out.vcf.gz ~{sample_id}_kanpig.vcf.gz ; tabix -f ~{sample_id}_kanpig.vcf.gz
        
        # Keeping only records that are genotyped as present. This is important,
        # since truvari bench does not consider GTs when matching.
        ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type z ~{sample_id}_truvari.vcf.gz > out.vcf.gz
        rm -f ~{sample_id}_truvari.vcf.gz* ; mv out.vcf.gz ~{sample_id}_truvari.vcf.gz ; tabix -f ~{sample_id}_truvari.vcf.gz
        ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type z ~{sample_id}_kanpig.vcf.gz > out.vcf.gz
        rm -f ~{sample_id}_kanpig.vcf.gz* ; mv out.vcf.gz ~{sample_id}_kanpig.vcf.gz ; tabix -f ~{sample_id}_kanpig.vcf.gz
        
        # Benchmarking
        Benchmark ~{sample_id} ~{sample_id}_truvari.vcf.gz truvari_~{min_sv_length}bp
        Benchmark ~{sample_id} ~{sample_id}_kanpig.vcf.gz kanpig_~{min_sv_length}bp
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp '*_truvari_*.txt' ~{remote_dir}/truvari/precision_recall/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading truvari benchmarks. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp '*_kanpig_*.txt' ~{remote_dir}/~{min_n_samples}_samples/precision_recall/ && echo 0 || echo 1)
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
        
        Int min_sv_length
        Int max_sv_length
        
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
                TEST=$(gsutil -m cp ~{remote_indir}/${PROBAND_ID}_'*.bcf*' ~{remote_indir}/${FATHER_ID}_'*.bcf*' ~{remote_indir}/${MOTHER_ID}_'*.bcf*' . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading VCF files. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            else 
                break
            fi
        done
        
        # Restricting to autosomes. This is just for simplicity and is not
        # strictly necessary.
        TEST=$(ls *.vcf.gz && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${PROBAND_ID}_*.bcf) > ${PROBAND_ID}_autosomes.vcf.gz &
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${FATHER_ID}_*.bcf) > ${FATHER_ID}_autosomes.vcf.gz &
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${MOTHER_ID}_*.bcf) > ${MOTHER_ID}_autosomes.vcf.gz &
        else
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${PROBAND_ID}_*.vcf.gz) > ${PROBAND_ID}_autosomes.vcf.gz &
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${FATHER_ID}_*.vcf.gz) > ${FATHER_ID}_autosomes.vcf.gz &
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${MOTHER_ID}_*.vcf.gz) > ${MOTHER_ID}_autosomes.vcf.gz &
        fi
        wait
        tabix -f ${PROBAND_ID}_autosomes.vcf.gz &
        tabix -f ${FATHER_ID}_autosomes.vcf.gz &
        tabix -f ${MOTHER_ID}_autosomes.vcf.gz &
        wait
        echo ${PROBAND_ID}_autosomes.vcf.gz > list.txt
        echo ${FATHER_ID}_autosomes.vcf.gz >> list.txt
        echo ${MOTHER_ID}_autosomes.vcf.gz >> list.txt
        ls -laht
        
        # Remark: we keep all records from the three files, since otherwise
        # bcftools merge might transform a 0/0 (that would have disappeared)
        # into a ./., affecting the number of records over which Mendelian
        # error is computed.
        
        # Keeping only records in the given length range
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ${PROBAND_ID}_autosomes.vcf.gz > out.vcf.gz
        rm -f ${PROBAND_ID}_autosomes.vcf.gz* ; mv out.vcf.gz ${PROBAND_ID}_autosomes.vcf.gz ; tabix -f ${PROBAND_ID}_autosomes.vcf.gz
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ${FATHER_ID}_autosomes.vcf.gz > out.vcf.gz
        rm -f ${FATHER_ID}_autosomes.vcf.gz* ; mv out.vcf.gz ${FATHER_ID}_autosomes.vcf.gz ; tabix -f ${FATHER_ID}_autosomes.vcf.gz
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ${MOTHER_ID}_autosomes.vcf.gz > out.vcf.gz
        rm -f ${MOTHER_ID}_autosomes.vcf.gz* ; mv out.vcf.gz ${MOTHER_ID}_autosomes.vcf.gz ; tabix -f ${MOTHER_ID}_autosomes.vcf.gz
        
        # Merging records by ID, since the records in every VCF originate from
        # the same cohort VCF, which had distinct IDs.
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge id --output-type z --file-list list.txt > trio.vcf.gz
        ${TIME_COMMAND} tabix -f trio.vcf.gz
        rm -f ${PROBAND_ID}_* ${FATHER_ID}_* ${MOTHER_ID}_*
        ls -laht
        
        # Benchmarking
        Benchmark trio.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_all
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z trio.vcf.gz > tr.vcf.gz
        ${TIME_COMMAND} tabix -f tr.vcf.gz
        Benchmark tr.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_tr
        rm -f tr.vcf.gz*
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z trio.vcf.gz > not_tr.vcf.gz
        ${TIME_COMMAND} tabix -f not_tr.vcf.gz
        Benchmark not_tr.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_not_tr
        rm -f not_tr.vcf.gz*
        
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
