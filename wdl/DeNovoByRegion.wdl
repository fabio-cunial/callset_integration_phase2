version 1.0


# 
#
workflow DeNovoByRegion {
    input {
        File ped_tsv
        Int n_trios
        
        File intersample_vcf_gz
        File intersample_tbi
        
        File reference_fa
        File reference_fai
        File tandem_bed
        
        Int only_50_bp
        Int max_region_length = 10000
    }
    parameter_meta {
        intersample_vcf_gz: "Assumed to contain only a few samples, i.e. we do not subset to only the samples in the trios internally."
    }
    
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    scatter (i in range(n_trios)) {
        call Impl {
            input:
                ped_tsv = ped_tsv,
                ped_tsv_row = i+1,
                intersample_vcf_gz = intersample_vcf_gz,
                intersample_tbi = intersample_tbi,
            
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
            
                only_50_bp = only_50_bp,
                max_region_length = max_region_length
        }
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
task Impl {
    input {
        File ped_tsv
        Int ped_tsv_row
        File intersample_vcf_gz
        File intersample_tbi
        
        File tandem_bed
        File not_tandem_bed
        
        Int only_50_bp
        Int max_region_length
        
        Int n_cpu = 4
        Int ram_size_gb = 12
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(intersample_vcf_gz,"GB"))
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        # Preparing the VCF
        head -n ~{ped_tsv_row} ~{ped_tsv} | tail -n 1 > ped.tsv
        PROBAND_ID=$(cut -f 2 ped.tsv)
        FATHER_ID=$(cut -f 3 ped.tsv)
        MOTHER_ID=$(cut -f 4 ped.tsv)
        if [ ~{only_50_bp} -ne 0 ]; then
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'SVLEN>=50 || SVLEN<=-50' --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{intersample_vcf_gz} > tmp.vcf.gz
        else
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{intersample_vcf_gz} > tmp.vcf.gz
        fi
        tabix -f tmp.vcf.gz
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp.vcf.gz > trio.vcf.gz
        tabix -f trio.vcf.gz
        rm -f tmp.vcf.gz*
        
        # Counting
        ${TIME_COMMAND} bcftools +mendelian2 trio.vcf.gz -m a -P ped.tsv -O z -o annotated.vcf.gz
        tabix -f annotated.vcf.gz
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z annotated.vcf.gz > tr.vcf.gz &
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z annotated.vcf.gz > not_tr.vcf.gz &
        wait
        tabix -f tr.vcf.gz &
        tabix -f not_tr.vcf.gz &
        wait
        ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigRegions annotated.vcf.gz ~{max_region_length} > ${PROBAND_ID}_all.bed &
        ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigRegions tr.vcf.gz ~{max_region_length} > ${PROBAND_ID}_tr.bed &
        ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigRegions not_tr.vcf.gz ~{max_region_length} > ${PROBAND_ID}_not_tr.bed &
        wait
    >>>
    
    output {
        Array[File] out_bed = glob("*.bed")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
