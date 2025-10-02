version 1.0


# 
#
workflow AnalyzeGtAdMatrix {
    input {
        File matrix_tsv
        File script_java
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            matrix_tsv = matrix_tsv,
            script_java = script_java
    }
    
    output {
        Array[File] histograms = Impl.histograms
    }
}


#
task Impl {
    input {
        File matrix_tsv
        File script_java
        
        Int n_cpu = 8
        Int ram_size_gb = 16
        Int disk_size_gb = 500
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))    
        
        
        if ~{defined(script_java)}; then
            mv ~{script_java} ./AnalyzeGtAdMatrix.java
            javac AnalyzeGtAdMatrix.java
        else
            mv ~{docker_dir}/AnalyzeGtAdMatrix.java .
            javac AnalyzeGtAdMatrix.java
        fi
        ${TIME_COMMAND} java -Xmx$(( ~{ram_size_gb} - 2 ))G AnalyzeGtAdMatrix ~{matrix_tsv}
        ls -laht
    >>>
    
    output {
        Array[File] histograms = glob("*.txt")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_squish"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
