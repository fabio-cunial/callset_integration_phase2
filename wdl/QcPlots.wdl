version 1.0


#
workflow QcPlots {
    input {
        File intersample_vcf_gz
        File intersample_tbi
    }
    parameter_meta {
    }
    
    call TypeVsLength {
        input:
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi
    }
    call SampleVsType as geq1 {
        input:
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi,
            min_sv_length = 1
    }
    call SampleVsType as geq50 {
        input:
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi,
            min_sv_length = 50
    }
    
    output {
        File type_length_tsv = TypeVsLength.out_tsv
        File sample_type_all = geq1.out_tsv
        File sample_type_50 = geq50.out_tsv
    }
}


# Performance on 10'070 samples, 15x, GRCh38, stringent (_S) and lenient (_L):
#
# TOOL                      CPU_S   RAM_S   TIME_S  CPU_L   RAM_L   TIME_L
# bcftools query
#
task TypeVsLength {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        
        Int n_cpu = 2
        Int ram_size_gb = 16
        Int disk_size_gb = 500
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        ${TIME_COMMAND} bcftools query -f '%SVTYPE\t%SVLEN\n' ~{intersample_vcf_gz} > out.tsv
    >>>

    output {
        File out_tsv = work_dir + "/out.tsv"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Performance on 10'070 samples, 15x, GRCh38, stringent (_S) and lenient (_L):
#
# TOOL                      CPU_S   RAM_S   TIME_S  CPU_L   RAM_L   TIME_L
# bcftools query
# sort
# uniq
#
task SampleVsType {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        Int min_sv_length
        
        Int n_cpu = 8
        Int ram_size_gb = 16
        Int disk_size_gb = 500
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        QUERY_STRING="GT==\"alt\" && (SVLEN>=~{min_sv_length} || SVLEN<=-~{min_sv_length})"
        ${TIME_COMMAND} bcftools query -i "${QUERY_STRING}" -f '[%SAMPLE\t%SVTYPE\n]' ~{intersample_vcf_gz} > tmp1.txt
        ${TIME_COMMAND} sort --parallel ${N_THREADS} tmp1.txt > tmp2.txt
        rm -f tmp1.txt
        ${TIME_COMMAND} uniq -c tmp2.txt > out.tsv
    >>>

    output {
        File out_tsv = work_dir + "/out.tsv"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}