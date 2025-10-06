version 1.0


# Same as `BenchCohortTriosSquish3`, just with more tandem BED files.
#
workflow BenchCohortTriosSquish4 {
    input {
        File ped_tsv
        Int n_trios
        Int only_50_bp
        
        File squish_vcf_gz
        File squish_tbi
        String squish_id
        
        Array[File] tandem_bed
        String bed_ids
        File reference_fai
        File? hardfilter_bed
    }
    parameter_meta {
        ped_tsv: "In the format used by `bcftools +mendelian2`: `<ignored>,proband,father,mother,sex(1:male,2:female)`."
        squish_id: "Unique ID of `squish_vcf_gz`"
        bed_ids: "Comma-separated"
        hardfilter_bed: "If defined, every VCF is hard-filtered to these regions before benchmarking."
    }
    
    scatter (i in range(n_trios)) {
        call BenchTrio {
            input:
                ped_tsv = ped_tsv,
                ped_tsv_row = i+1,
                only_50_bp = only_50_bp,
                
                squish_vcf_gz = squish_vcf_gz,
                squish_tbi = squish_tbi,
                squish_id = squish_id,
                
                tandem_bed = tandem_bed,
                bed_ids = bed_ids,
                reference_fai = reference_fai,
                hardfilter_bed = hardfilter_bed
        }
    }
    
    output {
        Array[Array[File]] out_txt = BenchTrio.out_txt
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
        
        File squish_vcf_gz
        File squish_tbi
        String squish_id
        
        Array[File] tandem_bed
        String bed_ids
        File reference_fai
        File? hardfilter_bed
        
        Int n_cpu = 16
        Int ram_size_gb = 64
    }
    parameter_meta {
        ped_tsv_row: "The row (one-based) in `ped_tsv` that corresponds to this trio."
        squish_id: "Unique ID of `squish_vcf_gz`"
        bed_ids: "Comma-separated"
        hardfilter_bed: "If defined, every VCF is hard-filtered to these regions before benchmarking."
    }
    
    Int disk_size_gb = 10*( ceil(size(squish_vcf_gz,"GB") + size(tandem_bed,"GB")) )
    
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
            
            if ~{defined(hardfilter_bed)}; then
                ${TIME_COMMAND} bcftools view --threads 1 --regions-file ~{hardfilter_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp0_${ID}.vcf.gz
                tabix -f tmp0_${ID}.vcf.gz
                HARDFILTERED_VCF_GZ="tmp0_${ID}.vcf.gz"
            else
                HARDFILTERED_VCF_GZ=${INPUT_VCF_GZ}
            fi
            ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z ${HARDFILTERED_VCF_GZ} > tmp1_${ID}.vcf.gz
            tabix -f tmp1_${ID}.vcf.gz
            ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp1_${ID}.vcf.gz > tmp2_${ID}.vcf.gz
            tabix -f tmp2_${ID}.vcf.gz
            truvari anno numneigh --sizemin 1 --refdist 1000 tmp2_${ID}.vcf.gz | bgzip > ${ID}.vcf.gz
            tabix -f ${ID}.vcf.gz
            rm -f ${INPUT_VCF_GZ} ${INPUT_VCF_GZ}.tbi tmp*_${ID}.vcf.gz*
        }
        
        
        # Evaluates an input VCF that contains only 3 samples (child, parents).
        # 
        function bench_thread() {
            local INPUT_VCF_GZ=$1
            local BED_ID=$2
            
            ${TIME_COMMAND} bcftools view --regions-file ${BED_ID}_sorted.bed --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${BED_ID}.vcf.gz
            tabix -f tmp_${BED_ID}.vcf.gz
            if [ ~{only_50_bp} -ne 0 ]; then
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${BED_ID}.vcf.gz -P ped.tsv --include 'SVLEN>=50 || SVLEN<=-50' > ${PROBAND_ID}_${BED_ID}.txt
                ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' -f '%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors,[%GT,%AD\t]\n' tmp_${BED_ID}.vcf.gz > ${PROBAND_ID}_${BED_ID}_gtmatrix.txt
            else
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${BED_ID}.vcf.gz -P ped.tsv > ${PROBAND_ID}_${BED_ID}_all.txt
                ${TIME_COMMAND} bcftools query -f '%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors,[%GT,%AD\t]\n' tmp_${BED_ID}.vcf.gz > ${PROBAND_ID}_${BED_ID}_gtmatrix.txt
            fi
            rm -f tmp_${BED_ID}.vcf.gz*
        }
        



        # Main program
        
        # Preprocessing the VCF
        head -n ~{ped_tsv_row} ~{ped_tsv} | tail -n 1 > ped.tsv
        PROBAND_ID=$(cut -f 2 ped.tsv)
        FATHER_ID=$(cut -f 3 ped.tsv)
        MOTHER_ID=$(cut -f 4 ped.tsv)
        source activate truvari5
        preprocess_thread ~{squish_vcf_gz} ~{squish_id}
        
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
        
        # Benchmarking
        while read ROW; do
            BED_ID=$(echo ${ROW} | cut -d , -f 2)
            bench_thread ~{squish_id}.vcf.gz ${BED_ID} &
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
