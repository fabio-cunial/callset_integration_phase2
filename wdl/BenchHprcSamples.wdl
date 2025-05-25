version 1.0


# Like `BenchCohortSamples.wdl`, but studies the per-sample stage of the
# pipeline, where we have the truvari collapse of all callers and its re-
# genotyped version with kanpig.
#
workflow BenchHprcSamples {
    input {
        Array[String] sample_ids = ["HG002", "HG00438", "HG005", "HG00621", "HG00673", "HG00733", "HG00735", "HG00741", "HG01071", "HG01106", "HG01109", "HG01123", "HG01175", "HG01243", "HG01258", "HG01358", "HG01361", "HG01891", "HG01928", "HG01952", "HG01978", "HG02055", "HG02080", "HG02109", "HG02145", "HG02148", "HG02257", "HG02486", "HG02559", "HG02572", "HG02622", "HG02630", "HG02717", "HG02723", "HG02818", "HG02886", "HG03098", "HG03453", "HG03486", "HG03492", "HG03516", "HG03540", "HG03579", "NA18906", "NA19240", "NA20129", "NA21309"]
        
        Array[File] truvari_vcf_gz
        Array[File] truvari_tbi
        
        Array[File] kanpig_vcf_gz
        Array[File] kanpig_tbi
        
        Array[File] dipcall_vcf_gz
        Array[File] dipcall_bed
        
        File tandem_bed
        File reference_fai
        
        Int min_sv_length
    }
    parameter_meta {
        truvari_vcf_gz: "In the same order as `sample_ids`."
        kanpig_vcf_gz: "In the same order as `sample_ids`."
        dipcall_vcf_gz: "In the same order as `sample_ids`."
        dipcall_bed: "In the same order as `sample_ids`."
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    scatter (i in range(length(sample_ids))) {
        call BenchSample as bench_truvari {
            input:
                sample_id = sample_ids[i],
                sample_vcf_gz = truvari_vcf_gz[i],
                sample_tbi = truvari_tbi[i],
                dipcall_vcf_gz = dipcall_vcf_gz[i],
                dipcall_bed = dipcall_bed[i],
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                min_sv_length = min_sv_length,
                only_present = 0
        }
        call BenchSample as bench_kanpig {
            input:
                sample_id = sample_ids[i],
                sample_vcf_gz = kanpig_vcf_gz[i],
                sample_tbi = kanpig_tbi[i],
                dipcall_vcf_gz = dipcall_vcf_gz[i],
                dipcall_bed = dipcall_bed[i],
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                min_sv_length = min_sv_length,
                only_present = 1
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
        
        Int n_cpu = 1
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 10*ceil(size(tandem_bed,"GB"))
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        
        ${TIME_COMMAND} bedtools sort -i ~{tandem_bed} -faidx ~{reference_fai} > sorted.bed
        ${TIME_COMMAND} bedtools complement -i sorted.bed -L -g ~{reference_fai} > complement.bed
    >>>
    
    output {
        File sorted_bed = work_dir + "/sorted.bed"
        File complement_bed = work_dir + "/complement.bed"
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
        
        File sample_vcf_gz
        File sample_tbi
        
        File dipcall_vcf_gz
        File dipcall_bed
        
        File tandem_bed
        File not_tandem_bed
        
        Int min_sv_length
        Int only_present
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 10*( ceil(size(sample_vcf_gz,"GB")) + ceil(size(dipcall_vcf_gz,"GB")) )
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        if [ ~{min_sv_length} -ne 0 ]; then
            FILTER_STRING="--sizemin ~{min_sv_length} --sizefilt ~{min_sv_length}"
        else
            FILTER_STRING="--sizemin 0 --sizefilt 0"
        fi
        
        # Making sure the dipcall file is normed
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z ~{dipcall_vcf_gz} > truth.vcf.gz
        ${TIME_COMMAND} tabix -f truth.vcf.gz
        rm -f ~{dipcall_vcf_gz}
        
        # Keeping only calls that occur in the given sample
        if [ ~{only_present} -eq 1 ]; then
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="0/1" || GT="0|1" || GT="1/0" || GT="1|0" || GT="1/1" || GT="1|1")>0' --output-type z ~{sample_vcf_gz} > sample.vcf.gz
            ${TIME_COMMAND} tabix -f sample.vcf.gz
            rm -f ~{sample_vcf_gz}*
        else
            mv ~{sample_vcf_gz} sample.vcf.gz
            mv ~{sample_tbi} sample.vcf.gz.tbi
        fi
        
        # Extracting calls with POS inside and outside TRs
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z sample.vcf.gz > sample_tr.vcf.gz &
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z sample.vcf.gz > sample_not_tr.vcf.gz &
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z truth.vcf.gz > truth_tr.vcf.gz &
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z truth.vcf.gz > truth_not_tr.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f sample_tr.vcf.gz &
        ${TIME_COMMAND} tabix -f sample_not_tr.vcf.gz &
        ${TIME_COMMAND} tabix -f truth_tr.vcf.gz &
        ${TIME_COMMAND} tabix -f truth_not_tr.vcf.gz &
        wait
        
        # Benchmarking
        ${TIME_COMMAND} truvari bench --includebed ~{dipcall_bed} -b truth.vcf.gz -c sample.vcf.gz ${FILTER_STRING} --sizemax 1000000 -o ./truvari_all/ &
        ${TIME_COMMAND} truvari bench --includebed ~{dipcall_bed} -b truth_tr.vcf.gz -c sample_tr.vcf.gz ${FILTER_STRING} --sizemax 1000000 -o ./truvari_tr/ &
        ${TIME_COMMAND} truvari bench --includebed ~{dipcall_bed} -b truth_not_tr.vcf.gz -c sample_not_tr.vcf.gz ${FILTER_STRING} --sizemax 1000000 -o ./truvari_not_tr/ &
        wait
        mv ./truvari_all/summary.json ./~{sample_id}_summary_all.json
        mv ./truvari_tr/summary.json ./~{sample_id}_summary_tr.json
        mv ./truvari_not_tr/summary.json ./~{sample_id}_summary_not_tr.json
    >>>
    
    output {
        File all_json = work_dir + "/" + sample_id + "_summary_all.json"
        File tr_json = work_dir + "/" + sample_id + "_summary_tr.json"
        File not_tr_json = work_dir + "/" + sample_id + "_summary_not_tr.json"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
