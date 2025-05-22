version 1.0


# Given a subset of samples, the workflow benchmarks the calls of a cohort VCF
# limited to those samples, against per-sample dipcall truth VCFs.
#
workflow BenchCohortSamples {
    input {
        Array[String] sample_ids = ["HG002", "HG00438", "HG005", "HG00621", "HG00673", "HG00733", "HG00735", "HG00741", "HG01071", "HG01106", "HG01109", "HG01123", "HG01175", "HG01243", "HG01258", "HG01358", "HG01361", "HG01891", "HG01928", "HG01952", "HG01978", "HG02055", "HG02080", "HG02109", "HG02145", "HG02148", "HG02257", "HG02486", "HG02559", "HG02572", "HG02622", "HG02630", "HG02717", "HG02723", "HG02818", "HG02886", "HG03098", "HG03453", "HG03486", "HG03492", "HG03516", "HG03540", "HG03579", "NA18906", "NA19240", "NA20129", "NA21309"]
        File cohort_vcf_gz
        File cohort_tbi
        Array[File] dipcall_vcf_gz
        Int min_sv_length
    }
    parameter_meta {
        dipcall_vcf_gz: "In the same order as `sample_ids`."
    }

    call SubsetToSamples {
        input:
            cohort_vcf_gz = cohort_vcf_gz,
            cohort_tbi = cohort_tbi,
            sample_ids = sample_ids
    }
    
    scatter (i in range(length(sample_ids))) {
        call BenchSample {
            input:
                sample_id = sample_ids[i],
                samples_vcf_gz = SubsetToSamples.out_vcf_gz,
                samples_tbi = SubsetToSamples.out_tbi,
                dipcall_vcf_gz = dipcall_vcf_gz[i],
                min_sv_length = min_sv_length
        }
    }
    
    output {
        Array[File] jsons = BenchSample.out_json
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
    
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 5*ceil(size(cohort_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        
        INPUT_FILES=~{sep=',' sample_ids}
        echo ${INPUT_FILES} | tr ',' '\n' > list.txt
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file list.txt --output-type z ~{cohort_vcf_gz} > tmp1.vcf.gz
        ${TIME_COMMAND} tabix -f tmp1.vcf.gz
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="0/1" || GT="0|1" || GT="1/0" || GT="1|0" || GT="1/1" || GT="1|1")>0' --output-type z tmp1.vcf.gz > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = work_dir + "/out.vcf.gz"
        File out_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Truvari bench with default parameters
#
task BenchSample {
    input {
        String sample_id
        File samples_vcf_gz
        File samples_tbi
        File dipcall_vcf_gz
        Int min_sv_length
        
        Int n_cpu = 4
        Int ram_size_gb = 4
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 5*( ceil(size(samples_vcf_gz,"GB")) + ceil(size(dipcall_vcf_gz,"GB")) )
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        # Keeping only calls that occur in the sample
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ~{sample_id} --output-type z ~{samples_vcf_gz} > tmp1.vcf.gz
        ${TIME_COMMAND} tabix -f tmp1.vcf.gz
        rm -f ~{samples_vcf_gz}
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="0/1" || GT="0|1" || GT="1/0" || GT="1|0" || GT="1/1" || GT="1|1")>0' --output-type z tmp1.vcf.gz > sample.vcf.gz
        ${TIME_COMMAND} tabix -f sample.vcf.gz
        rm -f tmp1.vcf.gz
        
        # Truvari
        ${TIME_COMMAND} tabix -f ~{dipcall_vcf_gz}
        if [ ~{min_sv_length} -ne 0 ]; then
            FILTER_STRING="--sizemin ~{min_sv_length} --sizefilt ~{min_sv_length}"
        else
            FILTER_STRING="--sizemin 0 --sizefilt 0"
        fi
        ${TIME_COMMAND} truvari bench -b ~{dipcall_vcf_gz} -c sample.vcf.gz ${FILTER_STRING} --sizemax 1000000 -o ./truvari/
        mv ./truvari/summary.json ./~{sample_id}_summary.json
    >>>
    
    output {
        File out_json = work_dir + "/" + sample_id + "_summary.json"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
