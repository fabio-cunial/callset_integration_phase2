version 1.0


# Given a subset of samples, the workflow benchmarks the calls of a cohort VCF
# limited to those samples, against per-sample dipcall truth VCFs.
#
workflow BenchCohortSamples {
    input {
        Array[String] sample_ids = ["HG002", "HG00438", "HG005", "HG00621", "HG00673", "HG00733", "HG00735", "HG00741", "HG01071", "HG01106", "HG01109", "HG01123", "HG01175", "HG01243", "HG01258", "HG01358", "HG01361", "HG01891", "HG01928", "HG01952", "HG01978", "HG02055", "HG02080", "HG02109", "HG02145", "HG02148", "HG02257", "HG02486", "HG02559", "HG02572", "HG02622", "HG02630", "HG02717", "HG02723", "HG02818", "HG02886", "HG03098", "HG03453", "HG03486", "HG03492", "HG03516", "HG03540", "HG03579", "NA18906", "NA19240", "NA20129", "NA21309"]
        Int min_sv_length
        
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
        File tandem_bed
        File reference_fai
    }
    parameter_meta {
        single_sample_dipcall_vcf_gz: "In the same order as `sample_ids`."
        single_sample_dipcall_bed: "In the same order as `sample_ids`."
        single_sample_kanpig_vcf_gz: "The output of kanpig without any further processing."
        single_sample_kanpig_annotated_vcf_gz: "The output of kanpig, annotated with `INFO/CALIBRATION_SENSITIVITY`."
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
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
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed
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
task BenchSample {
    input {
        String sample_id
        Int min_sv_length
        
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
        File tandem_bed
        File not_tandem_bed
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(single_sample_kanpig_vcf_gz,"GB")) + ceil(size(single_sample_kanpig_annotated_vcf_gz,"GB")) + ceil(size(cohort_merged_07_vcf_gz,"GB")) + ceil(size(cohort_merged_09_vcf_gz,"GB")) + ceil(size(cohort_regenotyped_07_vcf_gz,"GB")) + ceil(size(cohort_regenotyped_09_vcf_gz,"GB")) + ceil(size(single_sample_dipcall_vcf_gz,"GB")) )
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        if [ ~{min_sv_length} -ne 0 ]; then
            FILTER_STRING="--sizemin ~{min_sv_length} --sizefilt ~{min_sv_length}"
        else
            FILTER_STRING="--sizemin 0 --sizefilt 0"
        fi
        
        
        function bench_thread() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_PREFIX=$2
            
            # Extracting calls with POS inside and outside TRs
            ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > ${OUTPUT_PREFIX}_tr.vcf.gz
            ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > ${OUTPUT_PREFIX}_not_tr.vcf.gz
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_tr.vcf.gz
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_not_tr.vcf.gz
        
            # Benchmarking
            rm -rf ./${OUTPUT_PREFIX}_truvari_*
            ${TIME_COMMAND} truvari bench --includebed ~{single_sample_dipcall_bed} -b truth.vcf.gz -c ${INPUT_VCF_GZ} ${FILTER_STRING} --sizemax 1000000 -o ./${OUTPUT_PREFIX}_truvari_all/
            ${TIME_COMMAND} truvari bench --includebed ~{single_sample_dipcall_bed} -b truth_tr.vcf.gz -c ${OUTPUT_PREFIX}_tr.vcf.gz ${FILTER_STRING} --sizemax 1000000 -o ./${OUTPUT_PREFIX}_truvari_tr/
            ${TIME_COMMAND} truvari bench --includebed ~{single_sample_dipcall_bed} -b truth_not_tr.vcf.gz -c ${OUTPUT_PREFIX}_not_tr.vcf.gz ${FILTER_STRING} --sizemax 1000000 -o ./${OUTPUT_PREFIX}_truvari_not_tr/
            
            mv ./${OUTPUT_PREFIX}_truvari_all/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_all.json
            mv ./${OUTPUT_PREFIX}_truvari_tr/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_tr.json
            mv ./${OUTPUT_PREFIX}_truvari_not_tr/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_not_tr.json
        }


        # Main program
        
        # Preprocessing: dipcall VCF.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z ~{single_sample_dipcall_vcf_gz} > truth.vcf.gz
        ${TIME_COMMAND} tabix -f truth.vcf.gz
        rm -f ~{single_sample_dipcall_vcf_gz}
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z truth.vcf.gz > truth_tr.vcf.gz &
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z truth.vcf.gz > truth_not_tr.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f truth_tr.vcf.gz &
        ${TIME_COMMAND} tabix -f truth_not_tr.vcf.gz &
        wait
        
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
        bench_thread sample_kanpig.vcf.gz kanpig &
        bench_thread sample_07.vcf.gz 07 &
        bench_thread sample_09.vcf.gz 09 &
        bench_thread sample_cohort_merged_07.vcf.gz cohort_merged_07 &
        bench_thread sample_cohort_merged_09.vcf.gz cohort_merged_09 &
        bench_thread sample_cohort_regenotyped_07.vcf.gz cohort_regenotyped_07 &
        bench_thread sample_cohort_regenotyped_09.vcf.gz cohort_regenotyped_09 &
        wait
    >>>
    
    output {
        Array[File] out_jsons = glob("*.json")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
