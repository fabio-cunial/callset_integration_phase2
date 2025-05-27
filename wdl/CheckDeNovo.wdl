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
    call GetSamples {
        input:
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi,
            samples = samples
    }
    call GetMatrix as all {
        input:
            intersample_vcf_gz = GetSamples.out_vcf_gz,
            intersample_tbi = GetSamples.out_tbi,
            tandem_track_bed = ComplementBed.sorted_bed,
            tandem_track_complement_bed = ComplementBed.complement_bed,
            region_mode = 0,
            only_50 = 0
    }
    call GetMatrix as all_50 {
        input:
            intersample_vcf_gz = GetSamples.out_vcf_gz,
            intersample_tbi = GetSamples.out_tbi,
            tandem_track_bed = ComplementBed.sorted_bed,
            tandem_track_complement_bed = ComplementBed.complement_bed,
            region_mode = 0,
            only_50 = 1
    }
    call GetMatrix as tr {
        input:
            intersample_vcf_gz = GetSamples.out_vcf_gz,
            intersample_tbi = GetSamples.out_tbi,
            tandem_track_bed = ComplementBed.sorted_bed,
            tandem_track_complement_bed = ComplementBed.complement_bed,
            region_mode = 1,
            only_50 = 0
    }
    call GetMatrix as tr_50 {
        input:
            intersample_vcf_gz = GetSamples.out_vcf_gz,
            intersample_tbi = GetSamples.out_tbi,
            tandem_track_bed = ComplementBed.sorted_bed,
            tandem_track_complement_bed = ComplementBed.complement_bed,
            region_mode = 1,
            only_50 = 1
    }
    call GetMatrix as not_tr {
        input:
            intersample_vcf_gz = GetSamples.out_vcf_gz,
            intersample_tbi = GetSamples.out_tbi,
            tandem_track_bed = ComplementBed.sorted_bed,
            tandem_track_complement_bed = ComplementBed.complement_bed,
            region_mode = 2,
            only_50 = 0
    }
    call GetMatrix as not_tr_50 {
        input:
            intersample_vcf_gz = GetSamples.out_vcf_gz,
            intersample_tbi = GetSamples.out_tbi,
            tandem_track_bed = ComplementBed.sorted_bed,
            tandem_track_complement_bed = ComplementBed.complement_bed,
            region_mode = 2,
            only_50 = 1
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
        docker: "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_denovo"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# To circumvent the following bug of `bcftools query --samples` or
# `--samples-file`, which occurs both with v1.21 and v1.18 in the cloud but
# that cannot be reproduced with v1.18 on an interactive VM or locally:
#
# free(): invalid next size (fast)
# Command terminated by signal 6
#
task GetSamples {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        
        String samples
        
        Int n_cpu = 4
        Int ram_size_gb = 32
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
        
        echo ~{samples} | tr ',' '\n' > samples.txt
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file samples.txt --output-type z ~{intersample_vcf_gz} > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>

    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_denovo"
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
        
        File tandem_track_bed
        File tandem_track_complement_bed
        
        Int region_mode
        Int only_50
        
        Int n_cpu = 4
        Int ram_size_gb = 32
    }
    parameter_meta {
        region_mode: "0: all; 1=TR; 2=not TR."
        only_50: "1=use only calls >=50bp."
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 10*ceil(size(intersample_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        if [ ~{region_mode} -eq 0 ]; then
            REGION_STRING=" "
        elif [ ~{region_mode} -eq 1 ]; then
            REGION_STRING="--regions-file ~{tandem_track_bed} --regions-overlap pos"
        else
            REGION_STRING="--regions-file ~{tandem_track_complement_bed} --regions-overlap pos"
        fi
        if [ ~{only_50} -eq 1 ]; then
            ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' ${REGION_STRING} -f '[%GT,]\n' ~{intersample_vcf_gz} > matrix.txt
        else
            ${TIME_COMMAND} bcftools query ${REGION_STRING} -f '[%GT,]\n' ~{intersample_vcf_gz} > matrix.txt
        fi
        
    >>>

    output {
        File matrix_all = "matrix.txt"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_denovo"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
