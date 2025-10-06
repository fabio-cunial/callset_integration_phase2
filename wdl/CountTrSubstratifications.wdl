version 1.0


# Just counts the number of calls in each BED from an inter-sample VCF.
#
workflow CountTrSubstratifications {
    input {
        File vcf_gz
        File tbi
        Int only_50_bp
        File reference_fai
        
        Array[File] tandem_bed
        String bed_ids
    }
    parameter_meta {
        bed_ids: "Comma-separated"
    }
    
    call Impl {
        input:
            vcf_gz = vcf_gz,
            tbi = tbi,
            only_50_bp = only_50_bp,
            reference_fai = reference_fai,
            tandem_bed = tandem_bed,
            bed_ids = bed_ids
    }
    
    output {
        Array[File] out_txt = Impl.out_txt
    }
}


#
task Impl {
    input {
        File vcf_gz
        File tbi
        Int only_50_bp
        File reference_fai
        
        Array[File] tandem_bed
        String bed_ids
        
        Int n_cpu = 16
        Int ram_size_gb = 64
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(vcf_gz,"GB")) )
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        # Sorting BED files in parallel
        echo ~{sep="," tandem_bed} | tr ',' '\n' > bed_files.txt
        echo ~{bed_ids} | tr ',' '\n' > bed_ids.txt
        paste -d , bed_files.txt bed_ids.txt > bed_list.txt
        while read ROW; do
            BED_FILE=$(echo ${ROW} | cut -d , -f 1)
            BED_ID=$(echo ${ROW} | cut -d , -f 2)
            ${TIME_COMMAND} bedtools sort -i ${BED_FILE} -faidx ~{reference_fai} > ${BED_ID}_sorted.bed &
        done < bed_list.txt
        wait
        while read ROW; do
            BED_FILE=$(echo ${ROW} | cut -d , -f 1)
            rm -f ${BED_FILE}
        done < bed_list.txt
        
        # Counting in parallel
        while read ROW; do
            BED_ID=$(echo ${ROW} | cut -d , -f 2)
            if [ ~{only_50_bp} -ne 0 ]; then
                ${TIME_COMMAND} bcftools view --include 'SVLEN>=50 || SVLEN<=-50' --regions-file ${BED_ID}_sorted.bed --regions-overlap pos ~{vcf_gz} | wc -l > ${BED_ID}_count.txt &
            else
                ${TIME_COMMAND} bcftools view --regions-file ${BED_ID}_sorted.bed --regions-overlap pos ~{vcf_gz} | wc -l > ${BED_ID}_count.txt &
            fi    
        done < bed_list.txt
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
