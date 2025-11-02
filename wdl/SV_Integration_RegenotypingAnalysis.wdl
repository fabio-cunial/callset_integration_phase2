version 1.0


# 
#
workflow RegenotypingAnalysis {
    input {
        File cohort_truvari_vcf_gz
        File cohort_truvari_tbi
        Array[Int] min_n_samples = [2, 3, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
        Int min_sv_length = 20
        
        Array[File] precision_recall_samples_scored_vcf_gz
        Array[File] precision_recall_samples_dipcall_vcf_gz
        Array[File] precision_recall_samples_dipcall_bed
        Int precision_recall_bench_method
        
        File mendelian_error_ped_tsv
        Int mendelian_error_n_trios

        File tandem_bed
        File reference_fa
        File reference_fai
    }
    parameter_meta {
        precision_recall_bench_method: "0=truvari bench, 1=vcfdist."
        precision_recall_samples_dipcall_vcf_gz: "In the same order as `precision_recall_samples`."
        precision_recall_samples_dipcall_bed: "In the same order as `sample_ids`."
    }
    
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    call CleanCohortVcf {
        input:
            cohort_truvari_vcf_gz = cohort_truvari_vcf_gz,
            cohort_truvari_tbi = cohort_truvari_tbi
    }
    call FilterCohortVcf_ByLength {
        input:
            cohort_truvari_vcf_gz = CleanCohortVcf.out_vcf_gz,
            cohort_truvari_tbi = CleanCohortVcf.out_tbi
    }
    scatter (i in range(length(min_n_samples))) {
        call FilterCohortVcf_ByNSamples {
            input:
                cohort_truvari_vcf_gz = FilterCohortVcf_ByLength.out_vcf_gz,
                cohort_truvari_tbi = FilterCohortVcf_ByLength.out_tbi,
                min_n_samples = min_n_samples[i]
        }
        scatter (i in range(length(precision_recall_samples))) {
            call PrecisionRecallAnalysis {
                input:
                --------->
                
                    sample_id = precision_recall_samples[i],
                    dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[i],
                    dipcall_bed = precision_recall_samples_dipcall_bed[i],
                
                
                    min_n_samples = min_n_samples,
                
                    v1_07_cohort_truvari_vcf_gz = v1.out_vcf_gz,
                    v1_07_cohort_truvari_tbi = v1.out_tbi,
                
                    min_sv_length = min_sv_length,
                    max_sv_length = max_sv_length,
                    bench_method = bench_method,
                
                    tandem_bed = ComplementBed.sorted_bed,
                    not_tandem_bed = ComplementBed.complement_bed,
                    reference_fa = reference_fa,
                    reference_fai = reference_fai
            }
        }
        scatter (i in range(mendelian_error_n_trios)) {
            call MendelianErrorAnalysis {
                input:
                --------->
                
                    sample_id = precision_recall_samples[i],
                    dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[i],
                    dipcall_bed = precision_recall_samples_dipcall_bed[i],
                
                
                    min_n_samples = min_n_samples,
                
                    v1_07_cohort_truvari_vcf_gz = v1.out_vcf_gz,
                    v1_07_cohort_truvari_tbi = v1.out_tbi,
                
                    min_sv_length = min_sv_length,
                    max_sv_length = max_sv_length,
                    bench_method = bench_method,
                
                    tandem_bed = ComplementBed.sorted_bed,
                    not_tandem_bed = ComplementBed.complement_bed,
                    reference_fa = reference_fa,
                    reference_fai = reference_fai
            }
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
        Int ram_size_gb = 4
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
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, SSD:
#
# TOOL                CPU     RAM     TIME
#
#
task CleanCohortVcf {
    input {
        File cohort_truvari_vcf_gz
        File cohort_truvari_tbi
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        cohort_truvari_vcf_gz: "The raw output of cohort-level truvari collapse."
    }
    
    Int disk_size_gb = 4*ceil(size(cohort_truvari_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        mv ~{cohort_truvari_vcf_gz} in.vcf.gz
        mv ~{cohort_truvari_tbi} in.vcf.gz.tbi
        
        # Computing the number of samples each record occurs in
        ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%COUNT(GT="alt")\n' in.vcf.gz | bgzip -c > annotations.tsv.gz
        tabix -s1 -b2 -e2 annotations.tsv.gz
        
        # Dropping genotypes
        ${TIME_COMMAND} bcftools view --drop-genotypes in.vcf.gz --output-type z out.vcf.gz
        rm -f in.vcf.gz* ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        # Annotating the number of samples
        echo '##INFO=<ID=N_SAMPLES,Number=1,Type=Integer,Description="Number of samples where the record was discovered">' > header.txt
        ${TIME_COMMAND} bcftools annotate --header-lines header.txt --annotations annotations.tsv.gz --columns CHROM,POS,~ID,REF,ALT,N_SAMPLES --output-type z in.vcf.gz > out.vcf.gz
        rm -f in.vcf.gz* ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        # Enforcing a distinct ID for every record
        bcftools view --header-only in.vcf.gz > header.txt
        N_ROWS=$(wc -l < header.txt)
        date
        (  head -n $(( ${N_ROWS} - 1 )) header.txt ; \
           echo '##INFO=<ID=ORIGINAL_ID,Number=1,Type=String,Description="Original ID from truvari collapse">' ; \
           echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" ; \
           bcftools view --no-header in.vcf.gz | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { gsub(/;/,"_",$3); printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s;ORIGINAL_ID=%s\tGT\t0/1\n",$1,$2,++i,$4,$5,$6,$7,$8,$3); }' \
        ) | bgzip --compress-level 1 > out.vcf.gz
        date
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, SSD:
#
# TOOL                CPU     RAM     TIME
#
#
task FilterCohortVcf_ByLength {
    input {
        File cohort_truvari_vcf_gz
        File cohort_truvari_tbi
        
        Int min_sv_length
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        cohort_truvari_vcf_gz: "Every record is assumed be already annotated with the correct SVLEN."
    }
    
    Int disk_size_gb = 4*ceil(size(cohort_truvari_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'ABS(SVLEN)>='~{min_sv_length} --output-type z ~{cohort_truvari_vcf_gz} > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Partitions `cohort_truvari_vcf_gz` into the subset of all records that occur
# in `>= min_n_samples`, and its complement.
#
# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, SSD:
#
# TOOL                CPU     RAM     TIME
# 
#
task FilterCohortVcf_ByNSamples {
    input {
        File cohort_truvari_vcf_gz
        File cohort_truvari_tbi
        
        Int min_n_samples
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 4*ceil(size(cohort_truvari_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'N_SAMPLES>='~{min_n_samples} --output-type z ~{cohort_truvari_vcf_gz} > frequent_~{min_n_samples}.vcf.gz &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'N_SAMPLES<'~{min_n_samples} --output-type z ~{cohort_truvari_vcf_gz} > infrequent_~{min_n_samples}.vcf.gz &
        wait
        ${TIME_COMMAND} tabix frequent_~{min_n_samples}.vcf.gz &
        ${TIME_COMMAND} tabix infrequent_~{min_n_samples}.vcf.gz &
        wait
    >>>
    
    output {
        File frequent_vcf_gz = "frequent_"+min_n_samples+".vcf.gz"
        File frequent_tbi = "frequent_"+min_n_samples+".vcf.gz.tbi"
        File infrequent_vcf_gz = "infrequent_"+min_n_samples+".vcf.gz"
        File infrequent_tbi = "infrequent_"+min_n_samples+".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999:
#
# TOOL                CPU     RAM     TIME
# 
#
task BuildPersonalizedVcf {
    input {
        File sample_vcf_gz
        Int min_sv_length
        
        File frequent_cohort_vcf_gz
        File frequent_cohort_tbi
        File infrequent_cohort_vcf_gz
        File infrequent_cohort_tbi
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        sample_vcf_gz: "The filename is assumed to be `SAMPLEID_kanpig.vcf.gz`. Every record is assumed be already annotated with the correct SVLEN, and to have been marked as present by kanpig."
        frequent_cohort_vcf_gz: "Assumed to already contain records >= min_sv_length"
        infrequent_cohort_vcf_gz: "Assumed to already contain records >= min_sv_length"
    }
    
    Int disk_size_gb = 4*ceil( size(frequent_cohort_vcf_gz,"GB") + size(infrequent_cohort_vcf_gz,"GB") )
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        INFINITY="1000000000"

        SAMPLE_ID=$(basename ~{sample_vcf_gz} .vcf.gz)
        SAMPLE_ID=${SAMPLE_ID%_*}
        
        mv ~{sample_vcf_gz} in.vcf.gz
        
        # Ensuring the required min length
        bcftools filter --threads ${N_THREADS} --include 'ABS(SVLEN)>='~{min_sv_length} --output-type z in.vcf.gz > out.vcf.gz
        rm -f in.vcf.gz* ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        # Selecting infrequent cohort records that occur in this sample
        ${TIME_COMMAND} truvari bench --sizemin 0 --sizemax ${INFINITY} --sizefilt 0 --pick multi --comp in.vcf.gz --base ~{infrequent_cohort_vcf_gz} --output ./truvari/ 
        mv ./truvari/tp-base.vcf.gz infrequent.vcf.gz
        tabix -f infrequent.vcf.gz
        rm -rf ./truvari/
        
        # Merging infrequent cohort records and frequent cohort records
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --rm-dups exact --output-type z infrequent.vcf.gz ~{frequent_cohort_vcf_gz} > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}




















#-------------------------------------------------------------------------------


# Restricts the cohort VCF to records that occur in a PED file, to speed up the
# following steps.
#
task FilterCohortVcfForTrios {
    input {
        File cohort_truvari_vcf_gz
        File cohort_truvari_tbi
        
        File mendelian_error_ped_tsv
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(cohort_truvari_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        # Computing the union of all samples
        cut -f 1 ~{mendelian_error_ped_tsv} >> tmp.txt
        cut -f 2 ~{mendelian_error_ped_tsv} >> tmp.txt
        cut -f 3 ~{mendelian_error_ped_tsv} >> tmp.txt
        sort tmp.txt | uniq > list.txt
        rm -f tmp.txt
        
        # Filtering
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file list.txt --output-type z ~{cohort_truvari_vcf_gz} > tmp1.vcf.gz
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




# Restricts a cohort VCF to records that occur in a given set of samples, to
# speed up the following steps.
#
task FilterCohortVcfForPrecisionRecall {
    input {
        File cohort_truvari_vcf_gz
        File cohort_truvari_tbi
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 5*ceil(size(cohort_truvari_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        
        INPUT_FILES=~{sep=',' precision_recall_samples}
        echo ${INPUT_FILES} | tr ',' '\n' > list.txt
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file list.txt --output-type z ~{cohort_truvari_vcf_gz} > tmp1.vcf.gz
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















# Remark: this uses `truvari bench` with default parameters.
#
# Performance with 16 cores and 32GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# truvari bench             
# vcfdist
#
task PrecisionRecallAnalysis {
    input {
        String sample_id
        File single_sample_dipcall_vcf_gz
        File single_sample_dipcall_bed
        
        String remote_input_dir
        Array[Int] min_n_samples
        
        File? v1_07_cohort_truvari_vcf_gz
        File? v1_07_cohort_truvari_tbi
        
        Int min_sv_length
        Int max_sv_length
        Int bench_method
        
        File tandem_bed
        File not_tandem_bed
        File reference_fa
        File reference_fai
        
        Int ram_size_gb = 32
        Int disk_size_gb = 512
    }
    parameter_meta {
    }
    
    Int n_personalized_vcfs = length(min_n_samples)
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        if [ ~{min_sv_length} -ne 0 ]; then
            FILTER_STRING_TRUVARI="--sizemin ~{min_sv_length} --sizefilt ~{min_sv_length} --sizemax ~{max_sv_length}"
            FILTER_STRING_VCFDIST="--sv-threshold ~{min_sv_length} --largest-variant ~{max_sv_length}"
        else
            FILTER_STRING_TRUVARI="--sizemin 0 --sizefilt 0 --sizemax ~{max_sv_length}"
            FILTER_STRING_VCFDIST="--sv-threshold 0 --largest-variant ~{max_sv_length}"
        fi
        # See https://github.com/TimD1/vcfdist/wiki/02-Parameters-and-Usage
        # Remark: `--max-supercluster-size` has to be >= `--largest-variant + 2`
        #         and we set it to 10002 to mimic kanpig's `--sizemax`.
        # Remark: we choose `--cluster gap` since it is faster. We choose 500
        #         to mimic kanpig inter-sample's `--neighdist` (the intra-
        #         sample value would be 1000, which might be too big).
        SV_STRING_VCFDIST="--cluster gap 500 --max-supercluster-size $((~{max_sv_length}+2)) --realign-query --realign-truth"
        
        
        function bench_thread() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_PREFIX=$2
            
            # Fltering by length, if needed.
            if [ ~{min_sv_length} -ne 0 ]; then
                truvari anno svinfo -m 1 ${INPUT_VCF_GZ} | bcftools view --include "(SVLEN>=~{min_sv_length} && SVLEN<=~{max_sv_length}) || (SVLEN<=-~{min_sv_length} && SVLEN>=-~{max_sv_length})" --output-type z > ${OUTPUT_PREFIX}_input.vcf.gz
                tabix -f ${OUTPUT_PREFIX}_input.vcf.gz
            else
                cp ${INPUT_VCF_GZ} ${OUTPUT_PREFIX}_input.vcf.gz
                cp ${INPUT_VCF_GZ}.tbi ${OUTPUT_PREFIX}_input.vcf.gz.tbi
            fi
            
            # Extracting calls with POS inside and outside TRs
            ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz > ${OUTPUT_PREFIX}_tr.vcf.gz
            ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz > ${OUTPUT_PREFIX}_not_tr.vcf.gz
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_tr.vcf.gz
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_not_tr.vcf.gz
        
            # Benchmarking
            if [ ~{bench_method} -eq 0 ]; then
                # Assumed to be one of many jobs running in parallel, so using
                # only one core.
                rm -rf ./${OUTPUT_PREFIX}_truvari_*
                ${TIME_COMMAND} truvari bench --includebed ~{single_sample_dipcall_bed} -b truth.vcf.gz -c ${OUTPUT_PREFIX}_input.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_all/
                ${TIME_COMMAND} truvari bench --includebed ~{single_sample_dipcall_bed} -b truth_tr.vcf.gz -c ${OUTPUT_PREFIX}_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_tr/
                ${TIME_COMMAND} truvari bench --includebed ~{single_sample_dipcall_bed} -b truth_not_tr.vcf.gz -c ${OUTPUT_PREFIX}_not_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_not_tr/
                mv ./${OUTPUT_PREFIX}_truvari_all/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_all.txt
                mv ./${OUTPUT_PREFIX}_truvari_tr/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_tr.txt
                mv ./${OUTPUT_PREFIX}_truvari_not_tr/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_not_tr.txt
            else
                # Assumed to be the only job running, so using all cores.
                rm -f ./${OUTPUT_PREFIX}_vcfdist_*
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_input.vcf.gz truth.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{single_sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_all/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_tr.vcf.gz truth_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{single_sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_tr/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_not_tr.vcf.gz truth_not_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{single_sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_not_tr/
                mv ./${OUTPUT_PREFIX}_vcfdist_all/precision-recall-summary.tsv ./~{sample_id}_${OUTPUT_PREFIX}_all.txt
                mv ./${OUTPUT_PREFIX}_vcfdist_tr/precision-recall-summary.tsv ./~{sample_id}_${OUTPUT_PREFIX}_tr.txt
                mv ./${OUTPUT_PREFIX}_vcfdist_not_tr/precision-recall-summary.tsv ./~{sample_id}_${OUTPUT_PREFIX}_not_tr.txt
            fi
            
            # Removing temporary files
            rm -f ${OUTPUT_PREFIX}_input.vcf.gz*
        }


        # Main program
        ls -laht
        df -h
        
        # Preprocessing the dipcall VCF
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z ~{single_sample_dipcall_vcf_gz} > truth.vcf.gz
        ${TIME_COMMAND} tabix -f truth.vcf.gz
        rm -f ~{single_sample_dipcall_vcf_gz}
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z truth.vcf.gz > truth_tr.vcf.gz &
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z truth.vcf.gz > truth_not_tr.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f truth_tr.vcf.gz &
        ${TIME_COMMAND} tabix -f truth_not_tr.vcf.gz &
        wait
        
        # Preprocessing the V1 cohort VCF
        if ~{defined(v1_07_cohort_truvari_vcf_gz)}
        then
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ~{sample_id} --output-type z ~{v1_07_cohort_truvari_vcf_gz} > tmp1_07.vcf.gz
            ${TIME_COMMAND} tabix -f tmp1_07.vcf.gz
            rm -f ~{v1_07_cohort_truvari_vcf_gz}
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z tmp1_07.vcf.gz > v1_07.vcf.gz
            ${TIME_COMMAND} tabix -f v1_07.vcf.gz
            rm -f tmp1_*.vcf.gz
        fi
        
        # Localizing the personalized VCFs
        MIN_N_SAMPLES=$(echo ~{sep="," min_n_samples} | tr ',' ' ')
        for M in ${MIN_N_SAMPLES}; do
            if [[ ${M} -eq 1 ]]; then
                SUFFIX=""
            else
                SUFFIX="s"
            fi
            gsutil -m cp ~{remote_input_dir}/${M}_sample${SUFFIX}/~{sample_id}_kanpig.vcf.gz ./${M}_tmp1.vcf.gz &
        done
        wait
        for M in ${MIN_N_SAMPLES}; do
            tabix -f ${M}_tmp1.vcf.gz &
        done
        wait
        for M in ${MIN_N_SAMPLES}; do
            ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z ${M}_tmp1.vcf.gz > ${M}.vcf.gz &
        done
        wait
        for M in ${MIN_N_SAMPLES}; do
            tabix -f ${M}.vcf.gz &
        done
        wait
        rm -f *_tmp1.vcf.gz*
        
        # Benchmarking
        if [ ~{bench_method} -eq 0 ]; then
            # Parallel for truvari (which is single-core).
            if ~{defined(v1_07_cohort_truvari_vcf_gz)}
            then
                bench_thread v1_07.vcf.gz v1_07 &
            fi
            for M in ${MIN_N_SAMPLES}; do
                bench_thread ${M}.vcf.gz ${M} &
            done
            wait
        else
            # Sequential for vcfdist, since it takes too much RAM.
            if ~{defined(v1_07_cohort_truvari_vcf_gz)}
            then
                bench_thread v1_07.vcf.gz v1_07
            fi
            for M in ${MIN_N_SAMPLES}; do
                bench_thread ${M}.vcf.gz ${M}
            done
        fi
    >>>
    
    output {
        Array[File] out_jsons = glob("*.txt")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_personalized_vcfs + 1
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
