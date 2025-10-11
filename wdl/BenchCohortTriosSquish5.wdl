version 1.0


# Same as `BenchCohortTriosSquish4`, but just filtering by SVLEN and SVTYPE
# inside TRs.
#
workflow BenchCohortTriosSquish5 {
    input {
        File ped_tsv
        Int n_trios
        
        File squish_vcf_gz
        File squish_tbi
        String squish_id
        
        File tandems_bed
        Int min_svlen
        Array[Int] svlen_threhsolds
        
        File reference_fai
    }
    parameter_meta {
        ped_tsv: "In the format used by `bcftools +mendelian2`: `<ignored>,proband,father,mother,sex(1:male,2:female)`."
        squish_id: "Unique ID of `squish_vcf_gz`"
    }
    
    scatter (i in range(n_trios)) {
        call BenchTrio {
            input:
                ped_tsv = ped_tsv,
                ped_tsv_row = i+1,
                
                squish_vcf_gz = squish_vcf_gz,
                squish_tbi = squish_tbi,
                squish_id = squish_id,
                
                tandems_bed = tandems_bed,
                min_svlen = min_svlen,
                svlen_threhsolds = svlen_threhsolds,
                
                reference_fai = reference_fai
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
        
        File squish_vcf_gz
        File squish_tbi
        String squish_id
        
        File tandems_bed
        Int min_svlen
        Array[Int] svlen_threhsolds
        
        File reference_fai
        
        Int n_cpu = 16
        Int ram_size_gb = 64
    }
    parameter_meta {
        ped_tsv_row: "The row (one-based) in `ped_tsv` that corresponds to this trio."
        squish_id: "Unique ID of `squish_vcf_gz`"
    }
    
    Int disk_size_gb = 10*( ceil(size(squish_vcf_gz,"GB")) )
    
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
            
            ${TIME_COMMAND} bcftools view --threads 1 --regions-file ~{tandems_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp0_${ID}.vcf.gz
            tabix -f tmp0_${ID}.vcf.gz
            ${TIME_COMMAND} bcftools view --threads 1 --samples ${PROBAND_ID},${FATHER_ID},${MOTHER_ID} --output-type z tmp0_${ID}.vc.gz > tmp1_${ID}.vcf.gz
            tabix -f tmp1_${ID}.vcf.gz
            ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z tmp1_${ID}.vcf.gz > tmp2_${ID}.vcf.gz
            tabix -f tmp2_${ID}.vcf.gz
            truvari anno svinfo -m 1 tmp2_${ID}.vcf.gz | bgzip > ${ID}.vcf.gz
            tabix -f ${ID}.vcf.gz
            rm -f ${INPUT_VCF_GZ} ${INPUT_VCF_GZ}.tbi tmp*_${ID}.vcf.gz*
        }
        
        
        # Evaluates an input VCF that contains only 3 samples (child, parents).
        # 
        function bench_thread() {
            local INPUT_VCF_GZ=$1
            local SVLEN_FROM=$2
            local SVLEN_TO=$3
            
            ${TIME_COMMAND} bcftools view --include 'SVTYPE=="INS" && (SVLEN>='${SVLEN_FROM}' && SVLEN<'${SVLEN_TO}')' --output-type z ${INPUT_VCF_GZ} > tmp_${SVLEN_TO}.vcf.gz
            tabix -f tmp_${SVLEN_TO}.vcf.gz
            ${TIME_COMMAND} bcftools +mendelian2 tmp_${SVLEN_TO}.vcf.gz -P ped.tsv > ${PROBAND_ID}_${SVLEN_TO}_ins.txt
            ${TIME_COMMAND} bcftools query -f '%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors,[%GT,%AD\t]\n' tmp_${SVLEN_TO}.vcf.gz > ${PROBAND_ID}_${SVLEN_TO}_ins_gtmatrix.txt
            rm -f tmp_${SVLEN_TO}.vcf.gz*
            
            ${TIME_COMMAND} bcftools view --include 'SVTYPE=="DEL" && (SVLEN>='${SVLEN_FROM}' && SVLEN<'${SVLEN_TO}')' --output-type z ${INPUT_VCF_GZ} > tmp_${SVLEN_TO}.vcf.gz
            tabix -f tmp_${SVLEN_TO}.vcf.gz
            ${TIME_COMMAND} bcftools +mendelian2 tmp_${SVLEN_TO}.vcf.gz -P ped.tsv > ${PROBAND_ID}_${SVLEN_TO}_del.txt
            ${TIME_COMMAND} bcftools query -f '%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors,[%GT,%AD\t]\n' tmp_${SVLEN_TO}.vcf.gz > ${PROBAND_ID}_${SVLEN_TO}_del_gtmatrix.txt
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
        
        # Benchmarking
        SVLEN_THRESHOLDS=$(echo ~{sep="," svlen_threhsolds} | tr ',' ' ')
        SVLEN_FROM=~{min_svlen}
        for SVLEN_TO in ${SVLEN_THRESHOLDS}; do
            bench_thread ~{squish_id}.vcf.gz ${SVLEN_FROM} ${SVLEN_TO} &
            SVLEN_FROM=SVLEN_TO
        done
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
