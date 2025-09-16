version 1.0


# Given single-sample VCFs re-genotyped with kanpig for several trios, the
# program recomputes genotypes with beta binomial, then merges the VCFs of each
# trio separately and measures Mendelian error rate.
#
workflow TestBetaBinomial3 {
    input {
        File ped_tsv
        Int n_trios
        Int only_50_bp
        
        String remote_input_dir
        
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
                
                remote_input_dir = remote_input_dir,
                
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
# 
#
task BenchTrio {
    input {
        File ped_tsv
        Int ped_tsv_row
        Int only_50_bp
        
        String remote_input_dir
        
        File tandem_bed
        File not_tandem_bed
        
        Int n_cpu = 4
        Int ram_size_gb = 32
        Int disk_size_gb = 200
    }
    parameter_meta {
        ped_tsv_row: "The row (one-based) in `ped_tsv` that corresponds to this trio."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        # Re-genotypes a single VCF
        #
        function beta_binomial_thread() {
            local INPUT_VCF_GZ=$1
            local ID=$2
            
            # Re-genotyping
            ${TIME_COMMAND} python3 ~{docker_dir}/genotype-beta-binomial-mixture.py --kanpig-vcf ${INPUT_VCF_GZ} --output-prefix ${ID}
        
            # Annotating the input VCF
            bcftools view --header-only ${INPUT_VCF_GZ} > ${ID}_annotations.vcf
            cat ${ID}.annot.tsv >> ${ID}_annotations.vcf
            rm -f ${ID}.annot.tsv
            ${TIME_COMMAND} bgzip ${ID}_annotations.vcf
            tabix -f ${ID}_annotations.vcf.gz
            ${TIME_COMMAND} bcftools annotate --threads 1 --columns CHROM,POS,REF,ALT,FORMAT/GT,FORMAT/GQ,FORMAT/SQ --annotations ${ID}_annotations.vcf.gz --output-type z ${INPUT_VCF_GZ} > ${ID}_betabinomial.vcf.gz
            tabix -f ${ID}_betabinomial.vcf.gz
            rm -f ${ID}_annotations.vcf.gz
        }
        
        
        # Runs bcftools merge and truvari collapse on 3 VCFs.
        #
        function merge_thread() {
            local PROBAND_VCF_GZ=$1
            local FATHER_VCF_GZ=$2
            local MOTHER_VCF_GZ=$3
            local ID=$4
            
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --output-type z ${PROBAND_VCF_GZ} ${FATHER_VCF_GZ} ${MOTHER_VCF_GZ} > ${ID}_tmp1.vcf.gz
            tabix -f ${ID}_tmp1.vcf.gz
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z ${ID}_tmp1.vcf.gz > ${ID}_tmp2.vcf.gz
            tabix -f ${ID}_tmp2.vcf.gz
            rm -f ${ID}_tmp1.vcf.gz*
            ${TIME_COMMAND} truvari collapse --sizemin 0 --sizemax 1000000 --gt off --keep maxqual --input ${ID}_tmp2.vcf.gz --output ${ID}_tmp3.vcf
            rm -f ${ID}_tmp2.vcf.gz*
            ${TIME_COMMAND} bcftools sort --max-mem $(( ~{ram_size_gb} / 2 ))G --output-type z ${ID}_tmp3.vcf > ${ID}_merged.vcf.gz
            tabix -f ${ID}_merged.vcf.gz
            rm -f ${ID}_tmp3.vcf.gz*
        }
        
        
        # Evaluates an input VCF that contains only 3 samples (child, parents).
        # 
        function bench_thread() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_PREFIX=$2
            local PROBAND_ID=$3
            
            if [ ~{only_50_bp} -ne 0 ]; then
                ${TIME_COMMAND} bcftools +mendelian2 ${INPUT_VCF_GZ} -P ped.tsv --include 'SVLEN>=50 || SVLEN<=-50' > ${PROBAND_ID}_${OUTPUT_PREFIX}_all.txt
                ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' -f '%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors,[%GT,%AD\t]\n' ${INPUT_VCF_GZ} > ${PROBAND_ID}_${OUTPUT_PREFIX}_all_gtmatrix.txt
                # Inside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_tr.vcf.gz -P ped.tsv --include 'SVLEN>=50 || SVLEN<=-50' > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr.txt
                ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' -f '%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors,[%GT,%AD\t]\n' tmp_${OUTPUT_PREFIX}_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz*
                # Outside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz -P ped.tsv --include 'SVLEN>=50 || SVLEN<=-50' > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr.txt
                ${TIME_COMMAND} bcftools query --include 'SVLEN>=50 || SVLEN<=-50' -f '%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors,[%GT,%AD\t]\n' tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz*
            else
                ${TIME_COMMAND} bcftools +mendelian2 ${INPUT_VCF_GZ} -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_all.txt
                ${TIME_COMMAND} bcftools query -f '%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors,[%GT,%AD\t]\n' ${INPUT_VCF_GZ} > ${PROBAND_ID}_${OUTPUT_PREFIX}_all_gtmatrix.txt
                # Inside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_tr.vcf.gz -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr.txt
                ${TIME_COMMAND} bcftools query -f '%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors,[%GT,%AD\t]\n' tmp_${OUTPUT_PREFIX}_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_tr.vcf.gz*
                # Outside TRs
                ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                tabix -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz
                ${TIME_COMMAND} bcftools +mendelian2 tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz -P ped.tsv > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr.txt
                ${TIME_COMMAND} bcftools query -f '%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors,[%GT,%AD\t]\n' tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz > ${PROBAND_ID}_${OUTPUT_PREFIX}_not_tr_gtmatrix.txt
                rm -f tmp_${OUTPUT_PREFIX}_not_tr.vcf.gz*
            fi
        }




        # --------------------------- Main program -----------------------------
        
        # Downloading the VCFs
        head -n ~{ped_tsv_row} ~{ped_tsv} | tail -n 1 > ped.tsv
        PROBAND_ID=$(cut -f 2 ped.tsv)
        FATHER_ID=$(cut -f 3 ped.tsv)
        MOTHER_ID=$(cut -f 4 ped.tsv)
        gsutil -m cp ~{remote_input_dir}/${PROBAND_ID}_kanpig.vcf.gz* .
        gsutil -m cp ~{remote_input_dir}/${FATHER_ID}_kanpig.vcf.gz* .
        gsutil -m cp ~{remote_input_dir}/${MOTHER_ID}_kanpig.vcf.gz* .
        ${TIME_COMMAND} truvari anno numneigh --sizemin 1 --refdist 1000 ${PROBAND_ID}_kanpig.vcf.gz | bgzip > ${PROBAND_ID}.vcf.gz &
        ${TIME_COMMAND} truvari anno numneigh --sizemin 1 --refdist 1000 ${FATHER_ID}_kanpig.vcf.gz | bgzip > ${FATHER_ID}.vcf.gz &
        ${TIME_COMMAND} truvari anno numneigh --sizemin 1 --refdist 1000 ${MOTHER_ID}_kanpig.vcf.gz | bgzip > ${MOTHER_ID}.vcf.gz &
        wait
        tabix -f ${PROBAND_ID}.vcf.gz &
        tabix -f ${FATHER_ID}.vcf.gz &
        tabix -f ${MOTHER_ID}.vcf.gz &
        wait
        
        # Re-genotyping
        source activate pyro-kanpig
        beta_binomial_thread ${PROBAND_ID}.vcf.gz ${PROBAND_ID} &
        beta_binomial_thread ${FATHER_ID}.vcf.gz ${FATHER_ID} &
        beta_binomial_thread ${MOTHER_ID}.vcf.gz ${MOTHER_ID} &
        wait
        conda deactivate
        
        # Merging
        merge_thread ${PROBAND_ID}_betabinomial.vcf.gz ${FATHER_ID}_betabinomial.vcf.gz ${MOTHER_ID}_betabinomial.vcf.gz beta_binomial &
        merge_thread ${PROBAND_ID}.vcf.gz ${FATHER_ID}.vcf.gz ${MOTHER_ID}.vcf.gz kanpig &
        wait
        
        # Benchmarking
        bench_thread beta_binomial_merged.vcf.gz beta_binomial ${PROBAND_ID} &
        bench_thread kanpig_merged.vcf.gz kanpig ${PROBAND_ID} &
        wait
    >>>
    
    output {
        Array[File] out_txt = glob("*.txt")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_squish"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
