version 1.0


# Same as `BenchCohortSamples.wdl`, but for every window in a BED file
# separately.
#
workflow BenchCohortSamples_windowed {
    input {
        Array[String] sample_ids = ["HG002", "HG00438", "HG005", "HG00621", "HG00673", "HG00733", "HG00735", "HG00741", "HG01071", "HG01106", "HG01109", "HG01123", "HG01175", "HG01243", "HG01258", "HG01358", "HG01361", "HG01891", "HG01928", "HG01952", "HG01978", "HG02055", "HG02080", "HG02109", "HG02145", "HG02148", "HG02257", "HG02486", "HG02559", "HG02572", "HG02622", "HG02630", "HG02717", "HG02723", "HG02818", "HG02886", "HG03098", "HG03453", "HG03486", "HG03492", "HG03516", "HG03540", "HG03579", "NA18906", "NA19240", "NA20129", "NA21309"]
        Int min_sv_length = 50
        Int max_sv_length = 10000
        Int bench_method
        
        Array[File] single_sample_kanpig_vcf_gz
        Array[File] single_sample_kanpig_annotated_vcf_gz
        
        File cohort_merged_07_vcf_gz
        File cohort_merged_07_tbi
        File cohort_merged_09_vcf_gz
        File cohort_merged_09_tbi
        
        File cohort_regenotyped_07_vcf_gz
        File cohort_regenotyped_07_tbi
        File cohort_regenotyped_09_vcf_gz
        File cohort_regenotyped_09_tbi
        
        Array[File] single_sample_dipcall_vcf_gz
        Array[File] single_sample_dipcall_bed
        File windows_bed
        File reference_fa
        File reference_fai
    }
    parameter_meta {
        bench_method: "0=truvari bench, 1=vcfdist."
        single_sample_dipcall_vcf_gz: "In the same order as `sample_ids`."
        single_sample_dipcall_bed: "In the same order as `sample_ids`."
        single_sample_kanpig_vcf_gz: "The output of kanpig without any further processing."
        single_sample_kanpig_annotated_vcf_gz: "The output of kanpig, annotated with `INFO/CALIBRATION_SENSITIVITY`."
    }
    
    
    call SubsetToSamples as cohort_merged_07 {
        input:
            cohort_vcf_gz = cohort_merged_07_vcf_gz,
            cohort_tbi = cohort_merged_07_tbi,
            sample_ids = sample_ids
    }
    call SubsetToSamples as cohort_merged_09 {
        input:
            cohort_vcf_gz = cohort_merged_09_vcf_gz,
            cohort_tbi = cohort_merged_09_tbi,
            sample_ids = sample_ids
    }
    call SubsetToSamples as cohort_regenotyped_07 {
        input:
            cohort_vcf_gz = cohort_regenotyped_07_vcf_gz,
            cohort_tbi = cohort_regenotyped_07_tbi,
            sample_ids = sample_ids
    }
    call SubsetToSamples as cohort_regenotyped_09 {
        input:
            cohort_vcf_gz = cohort_regenotyped_09_vcf_gz,
            cohort_tbi = cohort_regenotyped_09_tbi,
            sample_ids = sample_ids
    }
    scatter (i in range(length(sample_ids))) {
        call BenchSample {
            input:
                sample_id = sample_ids[i],
                min_sv_length = min_sv_length,
                max_sv_length = max_sv_length,
                bench_method = bench_method,
                single_sample_kanpig_vcf_gz = single_sample_kanpig_vcf_gz[i],
                single_sample_kanpig_annotated_vcf_gz = single_sample_kanpig_annotated_vcf_gz[i],
                cohort_merged_07_vcf_gz = cohort_merged_07.out_vcf_gz,
                cohort_merged_07_tbi = cohort_merged_07.out_tbi,
                cohort_merged_09_vcf_gz = cohort_merged_09.out_vcf_gz,
                cohort_merged_09_tbi = cohort_merged_09.out_tbi,
                cohort_regenotyped_07_vcf_gz = cohort_regenotyped_07.out_vcf_gz,
                cohort_regenotyped_07_tbi = cohort_regenotyped_07.out_tbi,
                cohort_regenotyped_09_vcf_gz = cohort_regenotyped_09.out_vcf_gz,
                cohort_regenotyped_09_tbi = cohort_regenotyped_09.out_tbi,
                single_sample_dipcall_vcf_gz = single_sample_dipcall_vcf_gz[i],
                single_sample_dipcall_bed = single_sample_dipcall_bed[i],
                windows_bed = windows_bed,
                reference_fa = reference_fa,
                reference_fai = reference_fai
        }
    }
    
    output {
        Array[Array[File]] out_jsons = BenchSample.out_jsons
    }
}


# Restricts the cohort VCF to records that occur in a given set of samples, to
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
        Int min_sv_length
        Int max_sv_length
        Int bench_method
        
        File single_sample_kanpig_vcf_gz
        File single_sample_kanpig_annotated_vcf_gz
        
        File cohort_merged_07_vcf_gz
        File cohort_merged_07_tbi
        File cohort_merged_09_vcf_gz
        File cohort_merged_09_tbi
        
        File cohort_regenotyped_07_vcf_gz
        File cohort_regenotyped_07_tbi
        File cohort_regenotyped_09_vcf_gz
        File cohort_regenotyped_09_tbi
        
        File single_sample_dipcall_vcf_gz
        File single_sample_dipcall_bed
        
        File windows_bed
        File reference_fa
        File reference_fai
        
        Int n_cpu = 4
        Int ram_size_gb = 30
    }
    parameter_meta {
        bench_method: "0=truvari bench, 1=vcfdist."
    }
    
    Int disk_size_gb = 10*( ceil(size(single_sample_kanpig_vcf_gz,"GB")) + ceil(size(single_sample_kanpig_annotated_vcf_gz,"GB")) + ceil(size(cohort_merged_07_vcf_gz,"GB")) + ceil(size(cohort_merged_09_vcf_gz,"GB")) + ceil(size(cohort_regenotyped_07_vcf_gz,"GB")) + ceil(size(cohort_regenotyped_09_vcf_gz,"GB")) + ceil(size(single_sample_dipcall_vcf_gz,"GB")) )
    
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
                ${TIME_COMMAND} bcftools view --include "(SVLEN>=~{min_sv_length} && SVLEN<=~{max_sv_length}) || (SVLEN<=-~{min_sv_length} && SVLEN>=-~{max_sv_length})" --output-type z ${INPUT_VCF_GZ} > ${OUTPUT_PREFIX}_input.vcf.gz
                tabix -f ${OUTPUT_PREFIX}_input.vcf.gz
            else
                cp ${INPUT_VCF_GZ} ${OUTPUT_PREFIX}_input.vcf.gz
                cp ${INPUT_VCF_GZ}.tbi ${OUTPUT_PREFIX}_input.vcf.gz.tbi
            fi
            
            i="0"
            while read ROW; do
                CHR=$(echo ${ROW} | cut -d , -f 1)
                START=$(echo ${ROW} | cut -d , -f 2)
                END=$(echo ${ROW} | cut -d , -f 3)
                echo -e "${CHR}\t${START}\t${END}" > ${OUTPUT_PREFIX}_window.bed
                cat ${OUTPUT_PREFIX}_window.bed
                
                if [ ~{bench_method} -eq 0 ]; then
                    # Assumed to be one of many jobs running in parallel, so
                    # using only one core.
                    rm -rf ./${OUTPUT_PREFIX}_truvari/
                    ${TIME_COMMAND} truvari bench --includebed ${OUTPUT_PREFIX}_window.bed -b truth.vcf.gz -c ${OUTPUT_PREFIX}_input.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari/
                    mv ./${OUTPUT_PREFIX}_truvari/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_${i}.txt
                else
                    # Assumed to be the only job running, so using all cores.
                    rm -f ./${OUTPUT_PREFIX}_vcfdist/
                    ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_input.vcf.gz truth.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ${OUTPUT_PREFIX}_window.bed --prefix ./${OUTPUT_PREFIX}_vcfdist/
                    mv ./${OUTPUT_PREFIX}_vcfdist/precision-recall-summary.tsv ./~{sample_id}_${OUTPUT_PREFIX}_${i}.txt
                fi
                i=$(( ${i} + 1 ))
            done < windows.csv
            
            # Removing temporary files
            rm -f ${OUTPUT_PREFIX}_input.vcf.gz*
        }


        # Main program
        cat ~{windows_bed} | tr '\t' ',' > windows.csv
        
        # Preprocessing: dipcall VCF.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z ~{single_sample_dipcall_vcf_gz} > truth.vcf.gz
        ${TIME_COMMAND} tabix -f truth.vcf.gz
        rm -f ~{single_sample_dipcall_vcf_gz}
        
        # Preprocessing the VCF after:
        # 1. intra-sample truvari -> kanpig
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z ~{single_sample_kanpig_vcf_gz} > sample_kanpig.vcf.gz
        ${TIME_COMMAND} tabix -f sample_kanpig.vcf.gz
        
        # Preprocessing the VCF after:
        # 1. intra-sample truvari -> kanpig
        # 2. scoring -> filtering 0.7 and 0.9
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "INFO/CALIBRATION_SENSITIVITY<=0.7" --output-type z ~{single_sample_kanpig_annotated_vcf_gz} > sample_07.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "INFO/CALIBRATION_SENSITIVITY<=0.9" --output-type z ~{single_sample_kanpig_annotated_vcf_gz} > sample_09.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f sample_07.vcf.gz &
        ${TIME_COMMAND} tabix -f sample_09.vcf.gz &
        wait
        
        # Preprocessing the VCF after:
        # 1. intra-sample truvari -> kanpig
        # 2. scoring -> filtering 0.7 and 0.9
        # 3. inter-sample truvari
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ~{sample_id} --output-type z ~{cohort_merged_07_vcf_gz} > tmp1_07.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ~{sample_id} --output-type z ~{cohort_merged_09_vcf_gz} > tmp1_09.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f tmp1_07.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp1_09.vcf.gz &
        rm -f ~{cohort_merged_07_vcf_gz} ~{cohort_merged_09_vcf_gz}
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1_07.vcf.gz > sample_cohort_merged_07.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1_09.vcf.gz > sample_cohort_merged_09.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f sample_cohort_merged_07.vcf.gz &
        ${TIME_COMMAND} tabix -f sample_cohort_merged_09.vcf.gz &
        wait
        rm -f tmp1_*.vcf.gz
        
        # Preprocessing the VCF after:
        # 1. intra-sample truvari -> kanpig
        # 2. scoring -> filtering 0.7 and 0.9
        # 3. inter-sample truvari
        # 4. re-genotyping
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ~{sample_id} --output-type z ~{cohort_regenotyped_07_vcf_gz} > tmp1_07.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ~{sample_id} --output-type z ~{cohort_regenotyped_09_vcf_gz} > tmp1_09.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f tmp1_07.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp1_09.vcf.gz &
        wait
        rm -f ~{cohort_regenotyped_07_vcf_gz} ~{cohort_regenotyped_09_vcf_gz}
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1_07.vcf.gz > sample_cohort_regenotyped_07.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1_09.vcf.gz > sample_cohort_regenotyped_09.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f sample_cohort_regenotyped_07.vcf.gz &
        ${TIME_COMMAND} tabix -f sample_cohort_regenotyped_09.vcf.gz &
        wait
        rm -f tmp1_*.vcf.gz
        
        # Benchmarking
        if [ ~{bench_method} -eq 0 ]; then
            # Parallel for truvari (which is single-core).
            bench_thread sample_kanpig.vcf.gz kanpig &
            bench_thread sample_07.vcf.gz 07 &
            bench_thread sample_09.vcf.gz 09 &
            bench_thread sample_cohort_merged_07.vcf.gz cohort_merged_07 &
            bench_thread sample_cohort_merged_09.vcf.gz cohort_merged_09 &
            bench_thread sample_cohort_regenotyped_07.vcf.gz cohort_regenotyped_07 &
            bench_thread sample_cohort_regenotyped_09.vcf.gz cohort_regenotyped_09 &
            wait
        else
            # Sequential for vcfdist, since it takes too much RAM.
            bench_thread sample_kanpig.vcf.gz kanpig
            bench_thread sample_07.vcf.gz 07
            bench_thread sample_09.vcf.gz 09
            bench_thread sample_cohort_merged_07.vcf.gz cohort_merged_07
            bench_thread sample_cohort_merged_09.vcf.gz cohort_merged_09
            bench_thread sample_cohort_regenotyped_07.vcf.gz cohort_regenotyped_07
            bench_thread sample_cohort_regenotyped_09.vcf.gz cohort_regenotyped_09
        fi
    >>>
    
    output {
        Array[File] out_jsons = glob("*.txt")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
