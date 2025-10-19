version 1.0


# Benchmarks personalized VCFs and the V1 cohort VCF, against per-sample
# dipcall truth VCFs. This is similar to `BenchCohortSamples.wdl`.
#
workflow BenchCohortSamples_PersonalizedCohortVcf {
    input {
        Array[String] sample_ids = ["HG03579", "NA18906", "NA19240", "NA20129", "NA21309"]
        Array[File] single_sample_dipcall_vcf_gz
        Array[File] single_sample_dipcall_bed
        
        String remote_input_dir
        Array[Int] min_n_samples = [2, 3, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
        
        File v1_07_cohort_vcf_gz
        File v1_07_cohort_tbi
        
        Int min_sv_length = 50
        Int max_sv_length = 10000
        Int bench_method

        File tandem_bed
        File reference_fa
        File reference_fai
    }
    parameter_meta {
        bench_method: "0=truvari bench, 1=vcfdist."
        single_sample_dipcall_vcf_gz: "In the same order as `sample_ids`."
        single_sample_dipcall_bed: "In the same order as `sample_ids`."
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    call SubsetToSamples as v1 {
        input:
            cohort_vcf_gz = v1_07_cohort_vcf_gz,
            cohort_tbi = v1_07_cohort_tbi,
            sample_ids = sample_ids
    }
    scatter (i in range(length(sample_ids))) {
        call BenchSample {
            input:
                sample_id = sample_ids[i],
                single_sample_dipcall_vcf_gz = single_sample_dipcall_vcf_gz[i],
                single_sample_dipcall_bed = single_sample_dipcall_bed[i],
                
                remote_input_dir = remote_input_dir,
                min_n_samples = min_n_samples,
                
                v1_07_cohort_vcf_gz = v1.out_vcf_gz,
                v1_07_cohort_tbi = v1.out_tbi,
                
                min_sv_length = min_sv_length,
                max_sv_length = max_sv_length,
                bench_method = bench_method,
                
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                reference_fa = reference_fa,
                reference_fai = reference_fai
        }
    }
    
    output {
        Array[Array[File]] out_jsons = BenchSample.out_jsons
    }
}


# Restricts a cohort VCF to records that occur in a given set of samples, to
# speed up the following steps.
#
task SubsetToSamples {
    input {
        File cohort_vcf_gz
        File cohort_tbi
        Array[String] sample_ids
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(cohort_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        
        INPUT_FILES=~{sep=',' sample_ids}
        echo ${INPUT_FILES} | tr ',' '\n' > list.txt
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file list.txt --output-type z ~{cohort_vcf_gz} > tmp1.vcf.gz
        ${TIME_COMMAND} tabix -f tmp1.vcf.gz
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1.vcf.gz > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


#
task ComplementBed {
    input {
        File tandem_bed
        File reference_fai
        
        Int n_cpu = 1
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(tandem_bed,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        
        ${TIME_COMMAND} bedtools sort -i ~{tandem_bed} -faidx ~{reference_fai} > sorted.bed
        ${TIME_COMMAND} bedtools complement -i sorted.bed -L -g ~{reference_fai} > complement.bed
    >>>
    
    output {
        File sorted_bed = "sorted.bed"
        File complement_bed = "complement.bed"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Remark: this uses `truvari bench` with default parameters.
#
# Performance with 16 cores and 32GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# truvari bench             
# vcfdist
#
task BenchSample {
    input {
        String sample_id
        File single_sample_dipcall_vcf_gz
        File single_sample_dipcall_bed
        
        String remote_input_dir
        Array[Int] min_n_samples
        
        File v1_07_cohort_vcf_gz
        File v1_07_cohort_tbi
        
        Int min_sv_length
        Int max_sv_length
        Int bench_method
        
        File tandem_bed
        File not_tandem_bed
        File reference_fa
        File reference_fai
        
        Int ram_size_gb = 32
    }
    parameter_meta {
    }
    
    Int n_personalized_vcfs = length(min_n_samples)
    Int disk_size_gb = ( 1+n_personalized_vcfs )*( ceil(size(v1_07_cohort_vcf_gz,"GB")) )
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        if [ ~{min_sv_length} -ne 0 ]; then
            FILTER_STRING_TRUVARI="--sizemin ~{min_sv_length} --sizefilt ~{min_sv_length} --sizemax ~{max_sv_length}"
            FILTER_STRING_VCFDIST="--sv-threshold ~{min_sv_length} --largest-variant ~{max_sv_length}"
        else
            FILTER_STRING_TRUVARI="--sizemin 0 --sizefilt 0 --sizemax ~{max_sv_length}"
            FILTER_STRING_VCFDIST="--sv-threshold 0 --largest-variant ~{max_sv_length}"
        fi
        # See https://github.com/TimD1/vcfdist/wiki/02-Parameters-and-Usage
        # Remark: `--max-supercluster-size` has to be >= `--largest-variant + 2`
        #         and we set it to 10002 to mimic kanpig's `--sizemax`.
        # Remark: we choose `--cluster gap` since it is faster. We choose 500
        #         to mimic kanpig inter-sample's `--neighdist` (the intra-
        #         sample value would be 1000, which might be too big).
        SV_STRING_VCFDIST="--cluster gap 500 --max-supercluster-size $((~{max_sv_length}+2)) --realign-query --realign-truth"
        
        
        function bench_thread() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_PREFIX=$2
            
            # Fltering by length, if needed.
            if [ ~{min_sv_length} -ne 0 ]; then
                truvari anno svinfo -m 1 ${INPUT_VCF_GZ} | bcftools view --include "(SVLEN>=~{min_sv_length} && SVLEN<=~{max_sv_length}) || (SVLEN<=-~{min_sv_length} && SVLEN>=-~{max_sv_length})" --output-type z > ${OUTPUT_PREFIX}_input.vcf.gz
                tabix -f ${OUTPUT_PREFIX}_input.vcf.gz
            else
                cp ${INPUT_VCF_GZ} ${OUTPUT_PREFIX}_input.vcf.gz
                cp ${INPUT_VCF_GZ}.tbi ${OUTPUT_PREFIX}_input.vcf.gz.tbi
            fi
            
            # Extracting calls with POS inside and outside TRs
            ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz > ${OUTPUT_PREFIX}_tr.vcf.gz
            ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz > ${OUTPUT_PREFIX}_not_tr.vcf.gz
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_tr.vcf.gz
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_not_tr.vcf.gz
        
            # Benchmarking
            if [ ~{bench_method} -eq 0 ]; then
                # Assumed to be one of many jobs running in parallel, so using
                # only one core.
                rm -rf ./${OUTPUT_PREFIX}_truvari_*
                ${TIME_COMMAND} truvari bench --includebed ~{single_sample_dipcall_bed} -b truth.vcf.gz -c ${OUTPUT_PREFIX}_input.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_all/
                ${TIME_COMMAND} truvari bench --includebed ~{single_sample_dipcall_bed} -b truth_tr.vcf.gz -c ${OUTPUT_PREFIX}_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_tr/
                ${TIME_COMMAND} truvari bench --includebed ~{single_sample_dipcall_bed} -b truth_not_tr.vcf.gz -c ${OUTPUT_PREFIX}_not_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_not_tr/
                mv ./${OUTPUT_PREFIX}_truvari_all/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_all.txt
                mv ./${OUTPUT_PREFIX}_truvari_tr/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_tr.txt
                mv ./${OUTPUT_PREFIX}_truvari_not_tr/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_not_tr.txt
            else
                # Assumed to be the only job running, so using all cores.
                rm -f ./${OUTPUT_PREFIX}_vcfdist_*
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_input.vcf.gz truth.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{single_sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_all/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_tr.vcf.gz truth_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{single_sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_tr/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_not_tr.vcf.gz truth_not_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{single_sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_not_tr/
                mv ./${OUTPUT_PREFIX}_vcfdist_all/precision-recall-summary.tsv ./~{sample_id}_${OUTPUT_PREFIX}_all.txt
                mv ./${OUTPUT_PREFIX}_vcfdist_tr/precision-recall-summary.tsv ./~{sample_id}_${OUTPUT_PREFIX}_tr.txt
                mv ./${OUTPUT_PREFIX}_vcfdist_not_tr/precision-recall-summary.tsv ./~{sample_id}_${OUTPUT_PREFIX}_not_tr.txt
            fi
            
            # Removing temporary files
            rm -f ${OUTPUT_PREFIX}_input.vcf.gz*
        }


        # Main program
        
        # Preprocessing the dipcall VCF
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z ~{single_sample_dipcall_vcf_gz} > truth.vcf.gz
        ${TIME_COMMAND} tabix -f truth.vcf.gz
        rm -f ~{single_sample_dipcall_vcf_gz}
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z truth.vcf.gz > truth_tr.vcf.gz &
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z truth.vcf.gz > truth_not_tr.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f truth_tr.vcf.gz &
        ${TIME_COMMAND} tabix -f truth_not_tr.vcf.gz &
        wait
        
        # Preprocessing the V1 cohort VCF
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ~{sample_id} --output-type z ~{v1_07_cohort_vcf_gz} > tmp1_07.vcf.gz
        ${TIME_COMMAND} tabix -f tmp1_07.vcf.gz
        rm -f ~{v1_07_cohort_vcf_gz}
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1_07.vcf.gz > v1_07.vcf.gz
        ${TIME_COMMAND} tabix -f v1_07.vcf.gz
        rm -f tmp1_*.vcf.gz
        
        # Localizing the personalized VCFs
        MIN_N_SAMPLES=$(echo ~{sep="," min_n_samples} | tr ',' ' ')
        for M in ${MIN_N_SAMPLES}; do
            gsutil -m cp ~{remote_input_dir}/${M}_samples/~{sample_id}_kanpig.vcf.gz ./${M}_tmp1.vcf.gz &
        done
        wait
        for M in ${MIN_N_SAMPLES}; do
            tabix -f ${M}_tmp1.vcf.gz &
        done
        wait
        for M in ${MIN_N_SAMPLES}; do
            ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z ${M}_tmp1.vcf.gz > ${M}.vcf.gz &
        done
        wait
        for M in ${MIN_N_SAMPLES}; do
            tabix -f ${M}.vcf.gz &
        done
        wait
        rm -f *_tmp1.vcf.gz*
        
        # Benchmarking
        if [ ~{bench_method} -eq 0 ]; then
            # Parallel for truvari (which is single-core).
            bench_thread v1_07.vcf.gz v1_07 &
            for M in ${MIN_N_SAMPLES}; do
                bench_thread ${M}.vcf.gz ${M} &
            done
            wait
        else
            # Sequential for vcfdist, since it takes too much RAM.
            bench_thread v1_07.vcf.gz v1_07
            for M in ${MIN_N_SAMPLES}; do
                bench_thread ${M}.vcf.gz ${M}
            done
        fi
    >>>
    
    output {
        Array[File] out_jsons = glob("*.txt")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_personalized_vcfs + 1
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
