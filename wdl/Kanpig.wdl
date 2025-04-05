version 1.0


#
workflow Kanpig {
    input {
        Boolean is_singlesample
        String sample_id
        Boolean is_male = false
        String sex = "F"
        File vcf_gz
        File vcf_gz_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        Int n_cpu = 16
        Int ram_size_gb = 64
        String kanpig_params_singlesample = "--sizemin 20 --sizemax 10000 --neighdist 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"
        String kanpig_params_multisample  = "--sizemin 50 --sizemax 10000 --neighdist  500 --gpenalty 0.04 --hapsim 0.97"
    }
    parameter_meta {
        ram_size_gb: "Suggested 2*n_cpu"
    }
    
    call KanpigImpl {
        input:
            is_singlesample = is_singlesample,
            sample_id = sample_id,
            is_male = is_male,
            sex = sex,
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            kanpig_params_singlesample = kanpig_params_singlesample,
            kanpig_params_multisample = kanpig_params_multisample
    }
    
    output {
        File regenotyped_kanpig = KanpigImpl.regenotyped_kanpig_vcf
        File regenotyped_kanpig_tbi = KanpigImpl.regenotyped_kanpig_tbi
    }
}


task KanpigImpl {
    input {
        Boolean is_singlesample
        String sample_id
        Boolean is_male
        String sex
        File vcf_gz
        File vcf_gz_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        Int n_cpu
        Int ram_size_gb
        String kanpig_params_singlesample
        String kanpig_params_multisample
    }
    parameter_meta {
        ram_size_gb: "Suggested 2*n_cpu"
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    String output_prefix = sample_id + ".kanpig"
    Int disk_size_gb = 200 + ceil(size(reference_fa,"GB")) + 100*ceil(size(vcf_gz,"GB")) + 2*ceil(size(alignments_bam,"GB"))
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_MEM_GB=$(( ~{ram_size_gb} - 2 ))
        df -h
        
        
        # Unpacks truvari's SUPP field into 3 INFO tags.
        #
        function transferSupp() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF_GZ=$2
            
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SUPP]\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                printf("%s",$1); \
                for (i=2; i<=NF-1; i++) printf("\t%s",$i); \
                if ($6=="0") printf("\t0\t0\t0");
                else if ($6=="1") printf("\t0\t0\t1");
                else if ($6=="2") printf("\t0\t1\t0");
                else if ($6=="3") printf("\t0\t1\t1");
                else if ($6=="4") printf("\t1\t0\t0");
                else if ($6=="5") printf("\t1\t0\t1");
                else if ($6=="6") printf("\t1\t1\t0");
                else if ($6=="7") printf("\t1\t1\t1");
                printf("\n"); \
            }' | bgzip -c > annotations.tsv.gz
            tabix -f -s1 -b2 -e2 annotations.tsv.gz
            rm -f header.txt
            echo '##INFO=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by PAV">' > header.txt
            echo '##INFO=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by Sniffles2">' >> header.txt
            echo '##INFO=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> header.txt
            ${TIME_COMMAND} bcftools annotate --annotations annotations.tsv.gz --header-lines header.txt --columns CHROM,POS,ID,REF,ALT,INFO/SUPP_PAV,INFO/SUPP_SNIFFLES,INFO/SUPP_PBSV ${INPUT_VCF_GZ} --output-type z > ${OUTPUT_VCF_GZ}
            tabix -f ${OUTPUT_VCF_GZ}
        }


        # Main program
        touch ~{alignments_bai}
    
        # Formatting the merged VCF
        HAS_SUPP=$(bcftools view --header-only ~{vcf_gz} | grep '##FORMAT=<ID=SUPP,' && echo 1 || echo 0)
        if [ ${HAS_SUPP} -eq 0 ]; then
            mv ~{vcf_gz} formatted.vcf.gz
            mv ~{vcf_gz_tbi} formatted.vcf.gz.tbi
        else
            transferSupp ~{vcf_gz} formatted.vcf.gz
        fi

        # Making sure there is just one occurrence of '##fileformat=' in the
        # header (otherwise kanpig complains).
        rm -f cleaned.vcf*
        echo '##fileformat=VCFv4.2' > cleaned.vcf
        bcftools view --header-only formatted.vcf.gz | grep -v '##fileformat=' >> cleaned.vcf
        bcftools view --no-header formatted.vcf.gz >> cleaned.vcf
        bgzip -@ ${N_THREADS} cleaned.vcf
        tabix -f cleaned.vcf.gz
        rm -f formatted.vcf.gz*

        # Kanpig
        if [ ~{is_male} == "true" -o ~{sex} == "M" ]; then
            PLOIDY_BED=$(echo ~{ploidy_bed_male})
        else
            PLOIDY_BED=$(echo ~{ploidy_bed_female})
        fi
        if [ ~{is_singlesample} == "true" ]; then
            KANPIG_PARAMS=$(echo ~{kanpig_params_singlesample})
        else
            KANPIG_PARAMS=$(echo ~{kanpig_params_multisample})
        fi
        export RUST_BACKTRACE="full"
        ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ${KANPIG_PARAMS} --reference ~{reference_fa} --input cleaned.vcf.gz --reads ~{alignments_bam} --out tmp1.vcf.gz
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp1.vcf.gz > ~{output_prefix}.vcf.gz
        tabix -f ~{output_prefix}.vcf.gz
    >>>

    output {
        File regenotyped_kanpig_vcf = work_dir + "/" + output_prefix + ".vcf.gz"
        File regenotyped_kanpig_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
