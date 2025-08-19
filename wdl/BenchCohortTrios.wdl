version 1.0


# Measures Mendelian concordance at every step of the SV merging pipeline.
#
workflow BenchCohortTrios {
    input {
        File ped_tsv
        Int only_50_bp
        
        Array[File] single_sample_kanpig_proband_vcf_gz
        Array[File] single_sample_kanpig_father_vcf_gz
        Array[File] single_sample_kanpig_mother_vcf_gz
        
        Array[File] single_sample_kanpig_annotated_proband_vcf_gz
        Array[File] single_sample_kanpig_annotated_father_vcf_gz
        Array[File] single_sample_kanpig_annotated_mother_vcf_gz
        
        File cohort_merged_07_vcf_gz
        File cohort_merged_07_tbi
        File cohort_merged_09_vcf_gz
        File cohort_merged_09_tbi
        
        File cohort_regenotyped_07_vcf_gz
        File cohort_regenotyped_07_tbi
        File cohort_regenotyped_09_vcf_gz
        File cohort_regenotyped_09_tbi
        
        File tandem_bed
        File reference_fai
    }
    parameter_meta {
        ped_tsv: "In the format used by `bcftools +mendelian2`: `<ignored>,proband,father,mother,sex(1:male,2:female)`."
        single_sample_kanpig_proband_vcf_gz: "The output of kanpig without any further processing."
        single_sample_kanpig_annotated_proband_vcf_gz: "The output of kanpig, annotated with `INFO/CALIBRATION_SENSITIVITY`."
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    call SubsetToSamples as cohort_merged_07 {
        input:
            cohort_vcf_gz = cohort_merged_07_vcf_gz,
            cohort_tbi = cohort_merged_07_tbi,
            ped_tsv = ped_tsv
    }
    call SubsetToSamples as cohort_merged_09 {
        input:
            cohort_vcf_gz = cohort_merged_09_vcf_gz,
            cohort_tbi = cohort_merged_09_tbi,
            ped_tsv = ped_tsv
    }
    call SubsetToSamples as cohort_regenotyped_07 {
        input:
            cohort_vcf_gz = cohort_regenotyped_07_vcf_gz,
            cohort_tbi = cohort_regenotyped_07_tbi,
            ped_tsv = ped_tsv
    }
    call SubsetToSamples as cohort_regenotyped_09 {
        input:
            cohort_vcf_gz = cohort_regenotyped_09_vcf_gz,
            cohort_tbi = cohort_regenotyped_09_tbi,
            ped_tsv = ped_tsv
    }
    scatter (i in range(length(single_sample_kanpig_proband_vcf_gz))) {
        call BenchTrio {
            input:
                ped_tsv = ped_tsv,
                ped_tsv_row = i+1,
                only_50_bp = only_50_bp,
                single_sample_kanpig_proband_vcf_gz = single_sample_kanpig_proband_vcf_gz[i],
                single_sample_kanpig_father_vcf_gz = single_sample_kanpig_father_vcf_gz[i],
                single_sample_kanpig_mother_vcf_gz = single_sample_kanpig_mother_vcf_gz[i],
                single_sample_kanpig_annotated_proband_vcf_gz = single_sample_kanpig_annotated_proband_vcf_gz[i],
                single_sample_kanpig_annotated_father_vcf_gz = single_sample_kanpig_annotated_father_vcf_gz[i],
                single_sample_kanpig_annotated_mother_vcf_gz = single_sample_kanpig_annotated_mother_vcf_gz[i],
                cohort_merged_07_vcf_gz = cohort_merged_07.out_vcf_gz,
                cohort_merged_07_tbi = cohort_merged_07.out_tbi,
                cohort_merged_09_vcf_gz = cohort_merged_09.out_vcf_gz,
                cohort_merged_09_tbi = cohort_merged_09.out_tbi,
                cohort_regenotyped_07_vcf_gz = cohort_regenotyped_07.out_vcf_gz,
                cohort_regenotyped_07_tbi = cohort_regenotyped_07.out_tbi,
                cohort_regenotyped_09_vcf_gz = cohort_regenotyped_09.out_vcf_gz,
                cohort_regenotyped_09_tbi = cohort_regenotyped_09.out_tbi,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed
        }
    }
    
    output {
        Array[Array[File]] out_txt = BenchTrio.out_txt
    }
}


# Restricts the cohort VCF to records that occur in a given set of samples, to
# speed up the following steps.
#
task SubsetToSamples {
    input {
        File cohort_vcf_gz
        File cohort_tbi
        File ped_tsv
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(cohort_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        
        cut -f 2 ~{ped_tsv} > tmp.txt
        cut -f 3 ~{ped_tsv} >> tmp.txt
        cut -f 4 ~{ped_tsv} >> tmp.txt
        sort tmp.txt | uniq > list.txt
        rm -f tmp.txt
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file list.txt --output-type z ~{cohort_vcf_gz} > tmp1.vcf.gz
        ${TIME_COMMAND} tabix -f tmp1.vcf.gz
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1.vcf.gz > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
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
task BenchTrio {
    input {
        File ped_tsv
        Int ped_tsv_row
        Int only_50_bp
        
        File single_sample_kanpig_proband_vcf_gz
        File single_sample_kanpig_father_vcf_gz
        File single_sample_kanpig_mother_vcf_gz
        
        File single_sample_kanpig_annotated_proband_vcf_gz
        File single_sample_kanpig_annotated_father_vcf_gz
        File single_sample_kanpig_annotated_mother_vcf_gz
        
        File cohort_merged_07_vcf_gz
        File cohort_merged_07_tbi
        File cohort_merged_09_vcf_gz
        File cohort_merged_09_tbi
        
        File cohort_regenotyped_07_vcf_gz
        File cohort_regenotyped_07_tbi
        File cohort_regenotyped_09_vcf_gz
        File cohort_regenotyped_09_tbi
        
        File tandem_bed
        File not_tandem_bed
        
        Int n_cpu = 8
        Int ram_size_gb = 64
    }
    parameter_meta {
        ped_tsv_row: "The row (one-based) in `ped_tsv` that corresponds to this trio."
    }
    
    Int disk_size_gb = 10*( ceil(size(single_sample_kanpig_proband_vcf_gz,"GB")) + ceil(size(single_sample_kanpig_father_vcf_gz,"GB")) + ceil(size(single_sample_kanpig_mother_vcf_gz,"GB")) + ceil(size(single_sample_kanpig_annotated_proband_vcf_gz,"GB")) + ceil(size(single_sample_kanpig_annotated_father_vcf_gz,"GB")) + ceil(size(single_sample_kanpig_annotated_mother_vcf_gz,"GB")) + ceil(size(cohort_merged_07_vcf_gz,"GB")) + ceil(size(cohort_merged_09_vcf_gz,"GB")) + ceil(size(cohort_regenotyped_07_vcf_gz,"GB")) + ceil(size(cohort_regenotyped_09_vcf_gz,"GB")) )
    
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
        
        # Evaluating the VCFs after:
        # 1. intra-sample truvari -> kanpig
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z ~{single_sample_kanpig_proband_vcf_gz} > proband_kanpig.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z ~{single_sample_kanpig_father_vcf_gz} > father_kanpig.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z ~{single_sample_kanpig_mother_vcf_gz} > mother_kanpig.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f proband_kanpig.vcf.gz &
        ${TIME_COMMAND} tabix -f father_kanpig.vcf.gz &
        ${TIME_COMMAND} tabix -f mother_kanpig.vcf.gz &
        wait
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --output-type z proband_kanpig.vcf.gz father_kanpig.vcf.gz mother_kanpig.vcf.gz > trio_kanpig.vcf.gz
        ${TIME_COMMAND} tabix -f trio_kanpig.vcf.gz
        rm -f proband_kanpig.vcf.gz* father_kanpig.vcf.gz* mother_kanpig.vcf.gz*
        
        # Evaluating the VCF after:
        # 1. intra-sample truvari -> kanpig
        # 2. scoring -> filtering 0.7 and 0.9
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "INFO/CALIBRATION_SENSITIVITY<=0.7" --output-type z ~{single_sample_kanpig_annotated_proband_vcf_gz} > proband_07.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "INFO/CALIBRATION_SENSITIVITY<=0.7" --output-type z ~{single_sample_kanpig_annotated_father_vcf_gz} > father_07.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "INFO/CALIBRATION_SENSITIVITY<=0.7" --output-type z ~{single_sample_kanpig_annotated_mother_vcf_gz} > mother_07.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "INFO/CALIBRATION_SENSITIVITY<=0.9" --output-type z ~{single_sample_kanpig_annotated_proband_vcf_gz} > proband_09.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "INFO/CALIBRATION_SENSITIVITY<=0.9" --output-type z ~{single_sample_kanpig_annotated_father_vcf_gz} > father_09.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "INFO/CALIBRATION_SENSITIVITY<=0.9" --output-type z ~{single_sample_kanpig_annotated_mother_vcf_gz} > mother_09.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f proband_07.vcf.gz &
        ${TIME_COMMAND} tabix -f father_07.vcf.gz &
        ${TIME_COMMAND} tabix -f mother_07.vcf.gz &
        ${TIME_COMMAND} tabix -f proband_09.vcf.gz &
        ${TIME_COMMAND} tabix -f father_09.vcf.gz &
        ${TIME_COMMAND} tabix -f mother_09.vcf.gz &
        wait
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --output-type z proband_07.vcf.gz father_07.vcf.gz mother_07.vcf.gz > trio_07.vcf.gz &
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --output-type z proband_09.vcf.gz father_09.vcf.gz mother_09.vcf.gz > trio_09.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f trio_07.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_09.vcf.gz &
        wait
        rm -f proband_07.vcf.gz* father_07.vcf.gz* mother_07.vcf.gz* proband_09.vcf.gz* father_09.vcf.gz* mother_09.vcf.gz*
        
        # Preprocessing the VCF after:
        # 1. intra-sample truvari -> kanpig
        # 2. scoring -> filtering 0.7 and 0.9
        # 3. inter-sample truvari
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{cohort_merged_07_vcf_gz} > tmp1_07.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{cohort_merged_09_vcf_gz} > tmp1_09.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f tmp1_07.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp1_09.vcf.gz &
        wait
        rm -f ~{cohort_merged_07_vcf_gz} ~{cohort_merged_09_vcf_gz}
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1_07.vcf.gz > trio_cohort_merged_07.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1_09.vcf.gz > trio_cohort_merged_09.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f trio_cohort_merged_07.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_cohort_merged_09.vcf.gz &
        wait
        rm -f tmp1_*.vcf.gz*
        
        # Preprocessing the VCF after:
        # 1. intra-sample truvari -> kanpig
        # 2. scoring -> filtering 0.7 and 0.9
        # 3. inter-sample truvari
        # 4. re-genotyping
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{cohort_regenotyped_07_vcf_gz} > tmp1_07.vcf.gz &
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ~{cohort_regenotyped_09_vcf_gz} > tmp1_09.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f tmp1_07.vcf.gz &
        ${TIME_COMMAND} tabix -f tmp1_09.vcf.gz &
        wait
        rm -f ~{cohort_regenotyped_07_vcf_gz} ~{cohort_regenotyped_09_vcf_gz}
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1_07.vcf.gz > trio_cohort_regenotyped_07.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1_09.vcf.gz > trio_cohort_regenotyped_09.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f trio_cohort_regenotyped_07.vcf.gz &
        ${TIME_COMMAND} tabix -f trio_cohort_regenotyped_09.vcf.gz &
        wait
        rm -f tmp1_*.vcf.gz*
        
        # Benchmarking
        bench_thread trio_kanpig.vcf.gz kanpig &
        bench_thread trio_07.vcf.gz 07 &
        bench_thread trio_09.vcf.gz 09 &
        bench_thread trio_cohort_merged_07.vcf.gz cohort_merged_07 &
        bench_thread trio_cohort_merged_09.vcf.gz cohort_merged_09 &
        bench_thread trio_cohort_regenotyped_07.vcf.gz cohort_regenotyped_07 &
        bench_thread trio_cohort_regenotyped_09.vcf.gz cohort_regenotyped_09 &
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
