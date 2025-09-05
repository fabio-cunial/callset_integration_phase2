version 1.0


# Measures Mendelian concordance after some parametrizations of kanpig
# inter-sample.
#
workflow BenchCohortTriosSquish2 {
    input {
        File ped_tsv
        Int n_trios
        Int only_50_bp
        
        Array[File] squish_vcf_gz
        Array[File] squish_tbi
        String squish_ids
        
        File tandem_bed
        File reference_fai
    }
    parameter_meta {
        ped_tsv: "In the format used by `bcftools +mendelian2`: `<ignored>,proband,father,mother,sex(1:male,2:female)`."
        squish_ids: "Unique ID of each file in `squish_vcf_gz` (comma-separated)."
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
                
                squish_vcf_gz = squish_vcf_gz,
                squish_tbi = squish_tbi,
                squish_ids = squish_ids,
                
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
        
        Array[File] squish_vcf_gz
        Array[File] squish_tbi
        String squish_ids
        
        File tandem_bed
        File not_tandem_bed
        
        Int n_cpu = 22
        Int ram_size_gb = 128
    }
    parameter_meta {
        ped_tsv_row: "The row (one-based) in `ped_tsv` that corresponds to this trio."
        squish_ids: "Unique ID of each file in `squish_vcf_gz` (comma-separated)."
    }
    
    Int n_squish = length(squish_vcf_gz)
    Int disk_size_gb = 20*( ceil(size(squish_vcf_gz,"GB")) ) + 100
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        # Prepares an input VCF for benchmarking
        #
        function preprocess_thread() {
            local INPUT_VCF_GZ=$1
            local ID=$2
            
            ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ${INPUT_VCF_GZ} > tmp_${ID}.vcf.gz
            tabix -f tmp_${ID}.vcf.gz
            ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp_${ID}.vcf.gz > ${ID}.vcf.gz
            tabix -f ${ID}.vcf.gz
            rm -f ${INPUT_VCF_GZ} ${INPUT_VCF_GZ}.tbi tmp_${ID}.vcf.gz tmp_${ID}.vcf.gz.tbi
        }
        
        
        # Evaluates an input VCF that contains only 3 samples (child, parents).
        # 
        function bench_thread() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_PREFIX=$2
            
            if [ ~{only_50_bp} -ne 0 ]; then
                ${TIME_COMMAND} bcftools +mendelian2 ${INPUT_VCF_GZ} -P ped.tsv --include 'SVLEN>=50 || SVLEN<=-50' > ${PROBAND_ID}_${OUTPUT_PREFIX}_all.txt
                ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' -f '[%GT,%AD\t]\n' ${INPUT_VCF_GZ} > ${PROBAND_ID}_${OUTPUT_PREFIX}_all_gtmatrix.txt
                # Inside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_tr.vcf.gz -P ped.tsv --include 'SVLEN>=50 || SVLEN<=-50' > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr.txt
                ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' -f '[%GT,%AD\t]\n' tmp_${OUTPUT_PREFIX}_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz*
                # Outside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz -P ped.tsv --include 'SVLEN>=50 || SVLEN<=-50' > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr.txt
                ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' -f '[%GT,%AD\t]\n' tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz*
            else
                ${TIME_COMMAND} bcftools +mendelian2 ${INPUT_VCF_GZ} -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_all.txt
                ${TIME_COMMAND} bcftools query -f '[%GT,%AD\t]\n' ${INPUT_VCF_GZ} > ${PROBAND_ID}_${OUTPUT_PREFIX}_all_gtmatrix.txt
                # Inside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_tr.vcf.gz -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr.txt
                ${TIME_COMMAND} bcftools query -f '[%GT,%AD\t]\n' tmp_${OUTPUT_PREFIX}_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz*
                # Outside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr.txt
                ${TIME_COMMAND} bcftools query -f '[%GT,%AD\t]\n' tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz*
            fi
        }
        



        # Main program
        head -n ~{ped_tsv_row} ~{ped_tsv} | tail -n 1 > ped.tsv
        PROBAND_ID=$(cut -f 2 ped.tsv)
        FATHER_ID=$(cut -f 3 ped.tsv)
        MOTHER_ID=$(cut -f 4 ped.tsv)
        echo ~{sep="," squish_vcf_gz} | tr ',' '\n' > files.txt
        echo ~{squish_ids} | tr ',' '\n' > ids.txt
        paste -d , files.txt ids.txt > list.txt
        cat list.txt
        
        # Preprocessing the VCFs
        while read ROW; do
            VCF_FILE=$(echo ${ROW} | cut -d , -f 1)
            ID=$(echo ${ROW} | cut -d , -f 2)
            preprocess_thread ${VCF_FILE} ${ID} &
        done < list.txt
        wait
        
        # Benchmarking
        while read ROW; do
            ID=$(echo ${ROW} | cut -d , -f 2)
            bench_thread ${ID}.vcf.gz ${ID} &
        done < list.txt
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
