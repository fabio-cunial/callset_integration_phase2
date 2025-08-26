version 1.0


# Measures Mendelian concordance after some parametrizations of kanpig
# inter-sample.
#
workflow BenchCohortTriosSquish {
    input {
        File ped_tsv
        Int n_trios
        Int only_50_bp
        
        File cohort_regenotyped_07_vcf_gz
        File cohort_regenotyped_07_tbi
        
        File squish_vcf_gz
        File squish_tbi
        File ab_vcf_gz
        File ab_tbi
        File maxhom_vcf_gz
        File maxhom_tbi
        File fpenalty1_vcf_gz
        File fpenalty1_tbi
        File fpenalty2_vcf_gz
        File fpenalty2_tbi
        File gpenalty_vcf_gz
        File gpenalty_tbi
        File fnmax1_vcf_gz
        File fnmax1_tbi
        File fnmax2_vcf_gz
        File fnmax2_tbi
        
        File tandem_bed
        File reference_fai
    }
    parameter_meta {
        ped_tsv: "In the format used by `bcftools +mendelian2`: `<ignored>,proband,father,mother,sex(1:male,2:female)`."
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    scatter (i in range(n_trios)) {
        call BenchTrio {
            input:
                ped_tsv = ped_tsv,
                ped_tsv_row = i+1,
                only_50_bp = only_50_bp,
                
                cohort_regenotyped_07_vcf_gz = cohort_regenotyped_07_vcf_gz,
                cohort_regenotyped_07_tbi = cohort_regenotyped_07_tbi,
        
                squish_vcf_gz = squish_vcf_gz,
                squish_tbi = squish_tbi,
                ab_vcf_gz = ab_vcf_gz,
                ab_tbi = ab_tbi,
                maxhom_vcf_gz = maxhom_vcf_gz,
                maxhom_tbi = maxhom_tbi,
                fpenalty1_vcf_gz = fpenalty1_vcf_gz,
                fpenalty1_tbi = fpenalty1_tbi,
                fpenalty2_vcf_gz = fpenalty2_vcf_gz,
                fpenalty2_tbi = fpenalty2_tbi,
                gpenalty_vcf_gz = gpenalty_vcf_gz,
                gpenalty_tbi = gpenalty_tbi,
                fnmax1_vcf_gz = fnmax1_vcf_gz,
                fnmax1_tbi = fnmax1_tbi,
                fnmax2_vcf_gz = fnmax2_vcf_gz,
                fnmax2_tbi = fnmax2_tbi,
                
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed
        }
    }
    
    output {
        Array[Array[File]] out_txt = BenchTrio.out_txt
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


# Performance with 8 cores and 16GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# bcftools filter           100%        2G      1m
# bcftools +mendelian2      100%        50M     1m
#
task BenchTrio {
    input {
        File ped_tsv
        Int ped_tsv_row
        Int only_50_bp
        
        File cohort_regenotyped_07_vcf_gz
        File cohort_regenotyped_07_tbi
        
        File squish_vcf_gz
        File squish_tbi
        File ab_vcf_gz
        File ab_tbi
        File maxhom_vcf_gz
        File maxhom_tbi
        File fpenalty1_vcf_gz
        File fpenalty1_tbi
        File fpenalty2_vcf_gz
        File fpenalty2_tbi
        File gpenalty_vcf_gz
        File gpenalty_tbi
        File fnmax1_vcf_gz
        File fnmax1_tbi
        File fnmax2_vcf_gz
        File fnmax2_tbi
        
        File tandem_bed
        File not_tandem_bed
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        ped_tsv_row: "The row (one-based) in `ped_tsv` that corresponds to this trio."
    }
    
    Int disk_size_gb = 20*( ceil(size(cohort_regenotyped_07_vcf_gz,"GB")) ) + 100
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        # Evaluates an input VCF that contains only 3 samples (child, parents).
        # 
        function bench_thread() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_PREFIX=$2
            
            if [ ~{only_50_bp} -ne 0 ]; then
                ${TIME_COMMAND} bcftools +mendelian2 ${INPUT_VCF_GZ} -P ped.tsv --include 'SVLEN>=50 || SVLEN<=-50' > ${PROBAND_ID}_${OUTPUT_PREFIX}_all.txt
                ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' -f '[%GT,]\n' ${INPUT_VCF_GZ} > ${PROBAND_ID}_${OUTPUT_PREFIX}_all_gtmatrix.txt
                # Inside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_tr.vcf.gz -P ped.tsv --include 'SVLEN>=50 || SVLEN<=-50' > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr.txt
                ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' -f '[%GT,]\n' tmp_${OUTPUT_PREFIX}_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz*
                # Outside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz -P ped.tsv --include 'SVLEN>=50 || SVLEN<=-50' > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr.txt
                ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' -f '[%GT,]\n' tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz*
            else
                ${TIME_COMMAND} bcftools +mendelian2 ${INPUT_VCF_GZ} -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_all.txt
                ${TIME_COMMAND} bcftools query -f '[%GT,]\n' ${INPUT_VCF_GZ} > ${PROBAND_ID}_${OUTPUT_PREFIX}_all_gtmatrix.txt
                # Inside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_tr.vcf.gz -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr.txt
                ${TIME_COMMAND} bcftools query -f '[%GT,]\n' tmp_${OUTPUT_PREFIX}_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz*
                # Outside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr.txt
                ${TIME_COMMAND} bcftools query -f '[%GT,]\n' tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz*
            fi
        }
        

        # Main program
        head -n ~{ped_tsv_row} ~{ped_tsv} | tail -n 1 > ped.tsv
        PROBAND_ID=$(cut -f 2 ped.tsv)
        FATHER_ID=$(cut -f 3 ped.tsv)
        MOTHER_ID=$(cut -f 4 ped.tsv)
        
        # Preprocessing the VCFs
        ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{cohort_regenotyped_07_vcf_gz} > tmp1.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{squish_vcf_gz} > tmp2.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{ab_vcf_gz} > tmp3.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{maxhom_vcf_gz} > tmp4.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{fpenalty1_vcf_gz} > tmp5.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{fpenalty2_vcf_gz} > tmp6.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{gpenalty_vcf_gz} > tmp7.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{fnmax1_vcf_gz} > tmp8.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{fnmax2_vcf_gz} > tmp9.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f tmp1.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp2.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp3.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp4.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp5.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp6.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp7.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp8.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp9.vcf.gz &
        wait
        rm -f ~{cohort_regenotyped_07_vcf_gz} ~{squish_vcf_gz} ~{ab_vcf_gz} ~{maxhom_vcf_gz} ~{fpenalty1_vcf_gz} ~{fpenalty2_vcf_gz} ~{gpenalty_vcf_gz} ~{fnmax1_vcf_gz} ~{fnmax2_vcf_gz}
        ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp1.vcf.gz > trio_cohort_regenotyped_07.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp2.vcf.gz > trio_squish.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp3.vcf.gz > trio_ab.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp4.vcf.gz > trio_maxhom.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp5.vcf.gz > trio_fpenalty1.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp6.vcf.gz > trio_fpenalty2.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp7.vcf.gz > trio_gpenalty.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp8.vcf.gz > trio_fnmax1.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp9.vcf.gz > trio_fnmax2.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f trio_cohort_regenotyped_07.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_squish.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_ab.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_maxhom.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_fpenalty1.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_fpenalty2.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_gpenalty.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_fnmax1.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_fnmax2.vcf.gz &
        wait
        rm -f tmp*.vcf.gz*
        
        # Benchmarking
        bench_thread trio_cohort_regenotyped_07.vcf.gz cohort_regenotyped_07 &
        bench_thread trio_squish.vcf.gz squish &
        bench_thread trio_ab.vcf.gz ab &
        bench_thread trio_maxhom.vcf.gz maxhom &
        bench_thread trio_fpenalty1.vcf.gz fpenalty1 &
        bench_thread trio_fpenalty2.vcf.gz fpenalty2 &
        bench_thread trio_gpenalty.vcf.gz gpenalty &
        bench_thread trio_fnmax1.vcf.gz fnmax1 &
        bench_thread trio_fnmax2.vcf.gz fnmax2 &
        wait
    >>>
    
    output {
        Array[File] out_txt = glob("*.txt")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
