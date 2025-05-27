version 1.0


#
workflow CheckDeNovo {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        
        String samples
        
        File tandem_track_bed
        File reference_fai
    }
    parameter_meta {
        samples: "Comma-separated: child1,parent1_1,parent1_2,child2,parent2_1,parent2_2,..."
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_track_bed,
            reference_fai = reference_fai
    }
    call GetMatrix {
        input:
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi,
            samples = samples,
            tandem_track_bed = ComplementBed.sorted_bed,
            tandem_track_complement_bed = ComplementBed.complement_bed
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


# Performance on 10'070 samples, 15x, GRCh38, stringent (_S) and lenient (_L):
#
# TOOL                      CPU_S   RAM_S   TIME_S  CPU_L   RAM_L   TIME_L
#  
#
task GetMatrix {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        
        String samples
        
        File tandem_track_bed
        File tandem_track_complement_bed
        
        Int n_cpu = 8
        Int ram_size_gb = 64
    }
    parameter_meta {
        samples: "Comma-separated: child1,parent1_1,parent1_2,child2,parent2_1,parent2_2,..."
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 10*ceil(size(intersample_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        ${TIME_COMMAND} bcftools query --samples ~{samples} -f '[%GT ]\n' ~{intersample_vcf_gz} > matrix_all.txt &
        ${TIME_COMMAND} bcftools query --regions-file ~{tandem_track_bed} --regions-overlap pos --samples ~{samples} -f '[%GT ]\n' ~{intersample_vcf_gz} > matrix_tr.txt &
        ${TIME_COMMAND} bcftools query --regions-file ~{tandem_track_complement_bed} --regions-overlap pos --samples ~{samples} -f '[%GT ]\n' ~{intersample_vcf_gz} > matrix_not_tr.txt &
        ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' --samples ~{samples} -f '[%GT ]\n' ~{intersample_vcf_gz} > matrix_all_50.txt &
        ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' --regions-file ~{tandem_track_bed} --regions-overlap pos --samples ~{samples} -f '[%GT ]\n' ~{intersample_vcf_gz} > matrix_tr_50.txt &
        ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' --regions-file ~{tandem_track_complement_bed} --regions-overlap pos --samples ~{samples} -f '[%GT ]\n' ~{intersample_vcf_gz} > matrix_not_tr_50.txt &
        wait
    >>>

    output {
        File matrix_all = "matrix_all.txt"
        File matrix_tr = "matrix_tr.txt"
        File matrix_not_tr = "matrix_not_tr.txt"
        File matrix_all_50 = "matrix_all_50.txt"
        File matrix_tr_50 = "matrix_tr_50.txt"
        File matrix_not_tr_50 = "matrix_not_tr_50.txt"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
