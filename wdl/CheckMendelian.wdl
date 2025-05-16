version 1.0


#
workflow CheckMendelian {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        File ped_csv
        String? include
        
        Int n_cpu = 2
        Int ram_size_gb = 16
    }
    parameter_meta {
        ped_csv: "In the format used by `bcftools +mendelian2`: `<ignored>,proband,father,mother,sex(1:male,2:female)`."
    }
    
    call CheckMendelianImpl {
        input:
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi,
            ped_csv = ped_csv,
            include = include,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb
    }
    
    output {
        File out_txt = CheckMendelianImpl.out_txt
    }
}


# Performance on 10'070 samples, 15x, GRCh38, stringent (_S) and lenient (_L):
#
# TOOL                      CPU_S   RAM_S   TIME_S  CPU_L   RAM_L   TIME_L
# bcftools +mendelian2      
#
task CheckMendelianImpl {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        File ped_csv
        String? include
        
        Int n_cpu = 2
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 2*ceil(size(intersample_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        export BCFTOOLS_PLUGINS=~{docker_dir}/bcftools-1.21/plugins
        if ~{defined(include)}
        then
            INCLUDE_FLAG="--include ~{include}"
        else
            INCLUDE_FLAG=" "
        fi
        ${TIME_COMMAND} bcftools +mendelian2 ~{intersample_vcf_gz} -P ~{ped_csv} ${INCLUDE_FLAG} > out.txt
    >>>

    output {
        File out_txt = work_dir + "/out.txt"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
