version 1.0


# 
#
workflow GetNCalls {
    input {        
        Array[String] vcf_ids
        Array[File] vcf_gz
        Array[File] tbi
        
        File tandem_bed
        File reference_fai
    }
    parameter_meta {
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    scatter (i in range(length(vcf_gz))) {
        call GetNCalls {
            input:
                vcf_id = vcf_ids[i],
                vcf_gz = vcf_gz[i],
                tbi = tbi[i],
                
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed
        }
    }
    
    output {
        Array[Array[File]] out_txt = GetNCalls.out_txt
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
    
    Int disk_size_gb = 10*ceil(size(tandem_bed,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        
        ${TIME_COMMAND} bedtools sort -i ~{tandem_bed} -faidx ~{reference_fai} > sorted.bed
        ${TIME_COMMAND} bedtools complement -i sorted.bed -L -g ~{reference_fai} > complement.bed
    >>>
    
    output {
        File sorted_bed = "sorted.bed"
        File complement_bed = "complement.bed"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


#
task GetNCalls {
    input {
        String vcf_id
        File vcf_gz
        File tbi
        
        File tandem_bed
        File not_tandem_bed
        
        Int n_cpu = 6
        Int ram_size_gb = 12
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 20*( ceil(size(vcf_gz,"GB")) ) + 100
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
            
        # Total
        ${TIME_COMMAND} bcftools view --threads 1 --no-header                                     ${INPUT_VCF_GZ} | wc -l > ${vcf_id}_all.txt &
        ${TIME_COMMAND} bcftools view --threads 1 --no-header --include 'SVLEN>=50 || SVLEN<=-50' ${INPUT_VCF_GZ} | wc -l > ${vcf_id}_all_50.txt &
        # Inside TR
        ${TIME_COMMAND} bcftools view --threads 1 --no-header --regions-file ~{tandem_bed} --regions-overlap pos                                     ${INPUT_VCF_GZ} | wc -l > ${vcf_id}_tr_all.txt &
        ${TIME_COMMAND} bcftools view --threads 1 --no-header --regions-file ~{tandem_bed} --regions-overlap pos --include 'SVLEN>=50 || SVLEN<=-50' ${INPUT_VCF_GZ} | wc -l > ${vcf_id}_tr_50.txt &
        # Outside TR
        ${TIME_COMMAND} bcftools view --threads 1 --no-header --regions-file ~{not_tandem_bed} --regions-overlap pos                                     ${INPUT_VCF_GZ} | wc -l > ${vcf_id}_not_tr_all.txt &
        ${TIME_COMMAND} bcftools view --threads 1 --no-header --regions-file ~{not_tandem_bed} --regions-overlap pos --include 'SVLEN>=50 || SVLEN<=-50' ${INPUT_VCF_GZ} | wc -l > ${vcf_id}_not_tr_50.txt &
        wait
    >>>
    
    output {
        Array[File] out_txt = glob("*.txt")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
