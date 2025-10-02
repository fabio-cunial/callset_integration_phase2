version 1.0


# 
#
workflow GetGenotypingPriors {
    input {
        Array[String] ids
        Array[File] dipcall_vcf_gz
        Array[File] dipcall_tbi
        Array[File] kanpig_vcf_gz
        Array[File] kanpig_tbi
        Array[File] dipcall_bed
        
        File script_java
        String truvari_bench_args = ""
    }
    parameter_meta {
    }
    
    scatter(i in range(length(ids))) {
        call Impl {
            input:
                id = ids[i],
                dipcall_vcf_gz = dipcall_vcf_gz[i],
                dipcall_tbi = dipcall_tbi[i],
                kanpig_vcf_gz = kanpig_vcf_gz[i],
                kanpig_tbi = kanpig_tbi[i],
                dipcall_bed = dipcall_bed[i],
                script_java = script_java,
                truvari_bench_args = truvari_bench_args
        }
    }
    
    output {
        Array[Array[File]] distributions = Impl.distributions
    }
}


#
task Impl {
    input {
        String id
        File dipcall_vcf_gz
        File dipcall_tbi
        File kanpig_vcf_gz
        File kanpig_tbi
        File dipcall_bed
        
        File script_java
        String truvari_bench_args
        
        Int n_cpu = 4
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 10*( ceil(size(dipcall_vcf_gz,"GB")) + ceil(size(kanpig_vcf_gz,"GB")) ) + 100
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))    
        
        
        if ~{defined(script_java)}; then
            mv ~{script_java} ./GetGenotypingPriors.java
            javac GetGenotypingPriors.java
        else
            mv ~{docker_dir}/GetGenotypingPriors.java .
            javac GetGenotypingPriors.java
        fi
        mv ~{dipcall_vcf_gz} dipcall.vcf.gz
        mv ~{dipcall_tbi} dipcall.vcf.gz.tbi
        mv ~{kanpig_vcf_gz} kanpig.vcf.gz
        mv ~{kanpig_tbi} kanpig.vcf.gz.tbi
        ${TIME_COMMAND} truvari bench --sizemin 50 --sizefilt 50 ~{truvari_bench_args} --base dipcall.vcf.gz --comp kanpig.vcf.gz --includebed ~{dipcall_bed} --output ./truvari_out/
        ${TIME_COMMAND} java GetGenotypingPriors ./truvari_out/tp-comp.vcf.gz ./truvari_out/tp-base.vcf.gz ./truvari_out/fp.vcf.gz ./~{id}_distributions
        ls -laht
    >>>
    
    output {
        Array[File] distributions = glob(id+"_distributions_*.csv")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
