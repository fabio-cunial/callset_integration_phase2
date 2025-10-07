version 1.0


# Just counts the number of calls in each BED from an inter-sample VCF.
#
workflow CountTrSubstratifications {
    input {
        File vcf_gz
        File tbi
        Int only_50_bp
        File reference_fai
        
        Array[File] tandem_beds
        Array[String] bed_ids
        Int n_beds
    }
    parameter_meta {
    }
    
    scatter (i in range(n_beds)) {
        call Impl {
            input:
                only_50_bp = only_50_bp,
                
                vcf_gz = vcf_gz,
                tbi = tbi,
                reference_fai = reference_fai,
                
                bed_id = bed_ids[i],
                tandem_bed = tandem_beds[i]
        }
    }
    
    output {
        Array[Array[File]] out_txt = Impl.out_txt
    }
}


#
task Impl {
    input {
        Int only_50_bp
        
        File vcf_gz
        File tbi
        File reference_fai
        
        String bed_id
        File tandem_bed
        
        Int n_cpu = 4
        Int ram_size_gb = 32
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(vcf_gz,"GB") + size(tandem_bed,"GB")) )
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        # Sorting the BED file
        ${TIME_COMMAND} bedtools sort -i ~{tandem_bed} -faidx ~{reference_fai} > ~{bed_id}_sorted.bed
        rm -f ~{tandem_bed}
        
        # Counting
        if [ ~{only_50_bp} -ne 0 ]; then
            ${TIME_COMMAND} bcftools view --include 'SVLEN>=50 || SVLEN<=-50' --regions-file ~{bed_id}_sorted.bed --regions-overlap pos ~{vcf_gz} | wc -l > ~{bed_id}_count.txt
        else
            ${TIME_COMMAND} bcftools view --regions-file ~{bed_id}_sorted.bed --regions-overlap pos ~{vcf_gz} | wc -l > ~{bed_id}_count.txt
        fi
    >>>
    
    output {
        Array[File] out_txt = glob("*_count.txt")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
