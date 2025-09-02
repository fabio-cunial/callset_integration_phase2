version 1.0


# 
#
workflow GetGtAdMatrix {
    input {
        File intersample_vcf_gz
        File intersample_tbi
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi
    }
    
    output {
        File matrix_tsv = Impl.matrix_tsv
    }
}


#
task Impl {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        
        Int n_cpu = 8
        Int ram_size_gb = 16
        Int disk_size_gb = 500
    }
    parameter_meta {
    }
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))    
        
        
        ${TIME_COMMAND} bcftools query --format '[%GT,%AD\t]\n' > matrix.tsv
    >>>
    
    output {
        File matrix_tsv = "matrix.tsv"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_squish"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
