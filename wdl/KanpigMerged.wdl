version 1.0


# Re-genotypes an inter-sample VCF.
#
workflow KanpigMerged {
    input {
        String sample_id
        Boolean is_male = false
        String sex = "F"
        
        File intersample_vcf_gz
        File intersample_tbi
        File alignments_bam
        File alignments_bai
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        
        String kanpig_params_multisample = "--sizemin 20 --sizemax 10000 --neighdist 500 --gpenalty 0.04 --hapsim 0.97"
    }
    parameter_meta {
        intersample_vcf_gz: "Assumed to have a single (artificial) sample column and all GTs equal to 0/1."
    }

    call GetRegenotypedVcfImpl {
        input:
            sample_id = sample_id,
            is_male = is_male,
            sex = sex,
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            kanpig_params_multisample = kanpig_params_multisample
    }
    
    output {
        File output_vcf_gz = GetRegenotypedVcfImpl.output_vcf_gz
        File output_tbi = GetRegenotypedVcfImpl.output_tbi
        File gts_txt = GetRegenotypedVcfImpl.gts_txt
    }
}


#
task GetRegenotypedVcfImpl {
    input {
        String sample_id
        Boolean is_male
        String sex
        
        File intersample_vcf_gz
        File intersample_tbi
        File alignments_bam
        File alignments_bai
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        
        String kanpig_params_multisample
        
        Int n_cpu = 16
        Int mem_gb = 64
    }
    parameter_meta {
        intersample_vcf_gz: "Assumed to have a single (artificial) sample column and all GTs equal to 0/1."
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    String output_prefix = sample_id + "_intersample_regenotyped"
    Int disk_size_gb = 100 + ceil(size(reference_fa,"GB")) + 5*ceil(size(intersample_vcf_gz,"GB")) + 2*ceil(size(alignments_bam,"GB"))
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_MEM_GB=~{mem_gb}
        EFFECTIVE_MEM_GB=$(( ${EFFECTIVE_MEM_GB} - 4 ))
        df -h
        
        # - Enforcing the right sample ID
        echo ~{sample_id} > samples.txt
        ${TIME_COMMAND} bcftools reheader --threads ${N_THREADS} --samples samples.txt ~{intersample_vcf_gz} > tmp1.vcf.gz
        rm -f ~{intersample_vcf_gz}
        tabix -f tmp1.vcf.gz
        
        # - Re-genotyping
        touch ~{alignments_bai}
        if [ ~{is_male} == "true" -o ~{sex} == "M" ]; then
            PLOIDY_BED=$(echo ~{ploidy_bed_male})
        else
            PLOIDY_BED=$(echo ~{ploidy_bed_female})
        fi
        export RUST_BACKTRACE="full"
        ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_multisample} --reference ~{reference_fa} --input tmp1.vcf.gz --reads ~{alignments_bam} --out tmp2.vcf.gz
        rm -f tmp1.vcf.gz
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp2.vcf.gz > ~{output_prefix}.vcf.gz
        tabix -f ~{output_prefix}.vcf.gz
        rm -f tmp2.vcf.gz
        
        # - Building the GT-only output file
        echo ~{sample_id} > ~{sample_id}_gts.txt
        bcftools view --no-header ~{output_prefix}.vcf.gz | cut -f 10 >> ~{sample_id}_gts.txt
    >>>

    output {
        File output_vcf_gz = work_dir + "/" + output_prefix + ".vcf.gz"
        File output_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
        File gts_txt = work_dir + "/" + sample_id + "_gts.txt"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
