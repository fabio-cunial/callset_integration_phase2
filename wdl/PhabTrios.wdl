version 1.0


# 
#
workflow PhabTrios {
    input {
        File ped_tsv
        Int n_trios
        
        File intersample_vcf_gz
        File intersample_tbi
        
        File reference_fa
        File reference_fai
        File tandem_bed
        
        String phab_align = "wfa"
    }
    parameter_meta {
        intersample_vcf_gz: "Assumed to contain only a few samples, i.e. we do not subset to only the samples in the trios internally."
        phab_align: "Possible values: mafft, wfa, poa."
    }
    
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    call GetRegions {
        input:
            intersample_vcf_gz = intersample_vcf_gz
    }
    scatter (i in range(n_trios)) {
        call PhabTrio {
            input:
                ped_tsv = ped_tsv,
                ped_tsv_row = i+1,
                intersample_vcf_gz = intersample_vcf_gz,
                intersample_tbi = intersample_tbi,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
            
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
            
                phab_align = phab_align,
                regions_bed = GetRegions.regions_bed
        }
    }
    
    output {
    }
}


#
task GetRegions {
    input {
        File intersample_vcf_gz
    }
    parameter_meta {
        intersample_vcf_gz: "A VCF that has been re-genotyped with kanpig and carries all the corresponding annotations."
    }

    String docker_dir = "/callset_integration"
    Int disk_size_gb = 5*( ceil(size(intersample_vcf_gz,"GB")) )
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        
        
        ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigRegions ~{intersample_vcf_gz} > regions.bed
    >>>

    output {
        File regions_bed = "regions.bed"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: 1
        memory: "8GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
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
task PhabTrio {
    input {
        File ped_tsv
        Int ped_tsv_row
        File intersample_vcf_gz
        File intersample_tbi
        File reference_fa
        File reference_fai
        
        File tandem_bed
        File not_tandem_bed
        
        String phab_align
        File regions_bed
        
        Int n_cpu
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
        intersample_vcf_gz: "Assumed to contain only a few samples, i.e. we do not subset to only the samples in the trios internally."
        phab_align: "Possible values: mafft, wfa, poa."
    }
    
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
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{intersample_vcf_gz} > tmp.vcf.gz
        tabix -f tmp.vcf.gz
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp.vcf.gz > trio.vcf.gz
        tabix -f trio.vcf.gz
        rm -f tmp.vcf.gz*
        
        # Harmonizing
        ${TIME_COMMAND} truvari phab --threads ${N_THREADS} --align ~{phab_align} --baseline trio.vcf.gz --region ~{regions_bed} --reference ~{reference_fa} --output harmonized.vcf.gz
        tabix -f harmonized.vcf.gz
        rm -f trio.vcf.gz*
        
        # Benchmarking
        ${TIME_COMMAND} bcftools +mendelian2 harmonized.vcf.gz -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_all.txt
        ${TIME_COMMAND} bcftools query -f '%INFO/SVTYPE,%INFO/SVLEN,[%GT]\n' harmonized.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_all_gtmatrix.txt
        # Inside TRs
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_tr.vcf.gz
        tabix -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz
        ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_tr.vcf.gz -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr.txt
        ${TIME_COMMAND} bcftools query -f '%INFO/SVTYPE,%INFO/SVLEN,[%GT]\n' tmp_${OUTPUT_PREFIX}_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr_gtmatrix.txt
        rm -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz*
        # Outside TRs
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
        tabix -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
        ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr.txt
        ${TIME_COMMAND} bcftools query -f '%INFO/SVTYPE,%INFO/SVLEN,[%GT]\n' tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr_gtmatrix.txt
        rm -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz*
    >>>
    
    output {
        Array[File] out_txt = glob("*.txt")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
