version 1.0


#
workflow BuildTrainingResourceIntrasample {
    input {
        File regenotyped_vcf_gz
        File regenotyped_tbi
        File truth_vcf_gz
        File truth_vcf_gz_tbi
        File truth_bed
    }
    parameter_meta {
    }
    
    call BuildTrainingResourceIntrasampleImpl {
        input:
            regenotyped_vcf_gz = regenotyped_vcf_gz,
            regenotyped_tbi = regenotyped_tbi,
            truth_vcf_gz = truth_vcf_gz,
            truth_vcf_gz_tbi = truth_vcf_gz_tbi,
            truth_bed = truth_bed
    }
    
    output {
        File tp_comp_vcf = BuildTrainingResourceIntrasampleImpl.tp_comp_vcf
        File tp_comp_tbi = BuildTrainingResourceIntrasampleImpl.tp_comp_tbi
    }
}


#
task BuildTrainingResourceIntrasampleImpl {
    input {
        File regenotyped_vcf_gz
        File regenotyped_tbi
        File truth_vcf_gz
        File truth_vcf_gz_tbi
        File truth_bed
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    Int svlen_max = 1000000
    Int ram_gb = 32
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        TRUVARI_BENCH_SETTINGS="--sizemin 0 --sizemax ~{svlen_max} --sizefilt 0"
        
        mv ~{truth_vcf_gz} truth.vcf.gz
        mv ~{truth_vcf_gz_tbi} truth.vcf.gz.tbi
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS} --includebed ~{truth_bed} -b truth.vcf.gz -c ~{regenotyped_vcf_gz} -o truvari/
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z truvari/tp-comp.vcf.gz > tp_comp.vcf.gz
        tabix -f tp_comp.vcf.gz
    >>>

    output {
        File tp_comp_vcf = work_dir + "/tp_comp.vcf.gz"
        File tp_comp_tbi = work_dir + "/tp_comp.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: 16
        memory: ram_gb + "GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
