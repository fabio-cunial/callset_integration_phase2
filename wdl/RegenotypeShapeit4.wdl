version 1.0


# Structure of `remote_output_dir`:
#
# ├── shapeit4/                for each sample, its records in the shapeit4 VCF;
# │   ├── precision_recall/             for each sample, precision/recall stats;
# │   └── mendelian/                     for each sample, mendelian error stats;
#
workflow RegenotypeShapeit4 {
    input {
        File shapeit4_vcf_gz
        File shapeit4_tbi
        Int min_sv_length = 20
        Int max_sv_length = 10000
        String limit_to_chromosome = "chr6"
        
        String remote_dir
        
        Int precision_recall_bench_method
        Array[String] precision_recall_samples
        Array[String] precision_recall_sex
        Array[File] precision_recall_bam
        Array[File] precision_recall_bai
        Array[File] precision_recall_samples_dipcall_vcf_gz
        Array[File] precision_recall_samples_dipcall_bed
        
        Array[String] mendelian_error_samples
        Array[String] mendelian_error_sex
        Array[File] mendelian_error_bam
        Array[File] mendelian_error_bai
        File mendelian_error_ped
        Int mendelian_error_n_trios
        
        File reference_fa
        File reference_fai
        File reference_agp
        File standard_chromosomes_bed
        File autosomes_bed
        File tandem_bed
        File ploidy_bed_male
        File ploidy_bed_female
    }
    parameter_meta {
        min_sv_length: "The program subsets `cohort_truvari_vcf_gz` to this SVLEN range before every other step, and only afterwards re-genotypes. Regardless of this variable, benchmarking results are always reported for >=20bp and >=50bp records separately."
        precision_recall_bench_method: "0=truvari bench, 1=vcfdist."
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    call PrepareCohortBcf {
        input:
            cohort_shapeit4_vcf_gz = shapeit4_vcf_gz,
            cohort_shapeit4_tbi = shapeit4_tbi,
            limit_to_chromosome = limit_to_chromosome,
            min_sv_length = min_sv_length,
            max_sv_length = max_sv_length
    }
    call SplitBcf as split_shapeit4 {
        input:
            cohort_bcf = PrepareCohortBcf.out_bcf,
            cohort_csi = PrepareCohortBcf.out_csi,
            samples = flatten([precision_recall_samples, mendelian_error_samples]),
            remote_dir = remote_dir+"/shapeit4",
            suffix = "shapeit4",
            keep_only_present = 0,
            remove_sample = 0
    }
    scatter (i in range(mendelian_error_n_trios)) {
        call BenchTrio as me_bench_shapeit4_20 {
            input:
                ped_tsv = mendelian_error_ped,
                ped_tsv_row = i+1,
                remote_indir = remote_dir+"/shapeit4",
                remote_outdir = remote_dir+"/shapeit4/mendelian",
                min_sv_length = 20,
                max_sv_length = max_sv_length,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                autosomes_bed = autosomes_bed,
                in_flag = [split_shapeit4.out_flag]
        }
        call BenchTrio as me_bench_shapeit4_50 {
            input:
                ped_tsv = mendelian_error_ped,
                ped_tsv_row = i+1,
                remote_indir = remote_dir+"/shapeit4",
                remote_outdir = remote_dir+"/shapeit4/mendelian",
                min_sv_length = 50,
                max_sv_length = max_sv_length,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                autosomes_bed = autosomes_bed,
                in_flag = [split_shapeit4.out_flag]
        }
    }
    scatter (j in range(length(precision_recall_samples))) {
        call PrecisionRecallAnalysis as pr_analysis_20 {
            input:
                sample_id = precision_recall_samples[j],
                sample_dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[j],
                sample_dipcall_bed = precision_recall_samples_dipcall_bed[j],
                chromosome = limit_to_chromosome,
                
                remote_dir = remote_dir,
            
                bench_method = precision_recall_bench_method,
                min_sv_length = 20,
                max_sv_length = max_sv_length,    
            
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_agp = reference_agp,
                standard_chromosomes_bed = standard_chromosomes_bed,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                
                in_flag_shapeit4 = split_shapeit4.out_flag
        }
        call PrecisionRecallAnalysis as pr_analysis_50 {
            input:
                sample_id = precision_recall_samples[j],
                sample_dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[j],
                sample_dipcall_bed = precision_recall_samples_dipcall_bed[j],
                chromosome = limit_to_chromosome,
                
                remote_dir = remote_dir,
            
                bench_method = precision_recall_bench_method,
                min_sv_length = 50,
                max_sv_length = max_sv_length,    
            
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_agp = reference_agp,
                standard_chromosomes_bed = standard_chromosomes_bed,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
                
                in_flag_shapeit4 = split_shapeit4.out_flag
        }
    }
    
    output {
    }
}




#----------------------- Personalized VCF construction -------------------------

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
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"

        
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


# Keeps all sample columns in the output, and enforces a distinct ID in every
# record.
#
# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, SSD:
#
# TOOL                  CPU     RAM     TIME
# bcftools view b       200%    600M    15m
# bcftools query        100%    600M    3m
# bcftools annotate     100%    1G      15m
#
task PrepareCohortBcf {
    input {
        File cohort_shapeit4_vcf_gz
        File cohort_shapeit4_tbi
        
        String limit_to_chromosome
        Int min_sv_length
        Int max_sv_length
        
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
        limit_to_chromosome: "all=do not limit to any chromosome."
    }
    
    Int disk_size_gb = 4*ceil(size(cohort_shapeit4_vcf_gz,"GB"))
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        mv ~{cohort_shapeit4_vcf_gz} in.vcf.gz
        mv ~{cohort_shapeit4_tbi} in.vcf.gz.tbi
        
        # Removing SNVs
        ${TIME_COMMAND} bcftools filter --exclude 'STRLEN(REF)=1 && STRLEN(ALT)=1' --output-type z in.vcf.gz > out.vcf.gz
        rm -f in.vcf.gz* ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        # Making sure SVLEN and SVTYPE are consistently annotated
        ${TIME_COMMAND} truvari anno svinfo --minsize 1 in.vcf.gz | bgzip > out.vcf.gz
        rm -f in.vcf.gz* ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        # Converting to BCF, to speed up all steps downstream.
        if [ ~{limit_to_chromosome} = "all" ]; then
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type b in.vcf.gz > out.bcf
        else
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --regions ~{limit_to_chromosome} --output-type b in.vcf.gz > out.bcf
        fi
        rm -f in.vcf.gz* ; mv out.bcf in.bcf ; bcftools index in.bcf
        
        # Enforcing a distinct ID in every record, and annotating every record
        # with the number of samples it occurs in.
        ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%ID\t%COUNT(GT="alt")\n' in.bcf | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { $3=++i; gsub(/;/,"_",$6); print $0 }' | bgzip -c > annotations.tsv.gz
        tabix -s1 -b2 -e2 annotations.tsv.gz
        echo '##INFO=<ID=N_DISCOVERY_SAMPLES,Number=1,Type=Integer,Description="Number of samples where the record was discovered">' > header.txt
        echo '##INFO=<ID=ORIGINAL_ID,Number=1,Type=String,Description="Original ID from truvari collapse">' >> header.txt
        ${TIME_COMMAND} bcftools annotate --header-lines header.txt --annotations annotations.tsv.gz --columns CHROM,POS,ID,REF,ALT,ORIGINAL_ID,N_DISCOVERY_SAMPLES --output-type z in.bcf > out.bcf
        bcftools index out.bcf
    >>>
    
    output {
        File out_bcf = "out.bcf"
        File out_csi = "out.bcf.csi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Writes to a separate file every sample column of a cohort BCF.
#
# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, SSD:
#
# TOOL                              CPU     RAM     TIME
# bcftools +split                   100%    500M    3m
# bcftools filter                   100%    200M    10s
# bcftools view --drop-genotypes    100%    15M     10s
#
task SplitBcf {
    input {
        File cohort_bcf
        File cohort_csi
        
        Array[String] samples
        String remote_dir
        String suffix
        
        Int keep_only_present
        Int remove_sample
        
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
        cohort_bcf: "Assumed to have all the sample columns in the original truvari collapse VCF."
        remote_dir: "The result of the split is stored in this bucket location."
        remove_sample: "1=The output files do not have FORMAT and SAMPLE columns. 0=The output files have the original FORMAT and SAMPLE columns."
    }
    
    Int disk_size_gb = 10*ceil(size(cohort_bcf,"GB"))
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        # Splitting
        echo ~{sep="," samples} | tr ',' '\n' > samples.txt
        ${TIME_COMMAND} bcftools +split --samples-file samples.txt --output-type b --output . ~{cohort_bcf}
        rm -f ~{cohort_bcf}
        ls -laht
        
        # Keeping only present records, then removing FORMAT and SAMPLE.
        for FILE in $(ls *.bcf); do
            SAMPLE_ID=$(basename ${FILE} .bcf)
            if [ ~{keep_only_present} -eq 1 ]; then
                ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type b ${FILE} > tmp.bcf
            else
                mv ${FILE} tmp.bcf
            fi
            bcftools index tmp.bcf
            if [ ~{remove_sample} -eq 1 ]; then
                ${TIME_COMMAND} bcftools view --drop-genotypes --output-type b tmp.bcf > ${SAMPLE_ID}_~{suffix}.bcf
                bcftools index ${SAMPLE_ID}_~{suffix}.bcf
                rm -f tmp.bcf*
            else
                mv tmp.bcf ${SAMPLE_ID}_~{suffix}.bcf
                mv tmp.bcf.csi ${SAMPLE_ID}_~{suffix}.bcf.csi
            fi
        done
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp '*_'~{suffix}'.bcf*' ~{remote_dir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading VCFs. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        # Fake output
        echo "done" > out.txt
    >>>
    
    output {
        File out_flag = "out.txt"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}




#------------------------------- Benchmarking ----------------------------------

# Performance with 4 cores and 32GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# truvari bench             
# vcfdist
#
task PrecisionRecallAnalysis {
    input {
        String sample_id
        File sample_dipcall_vcf_gz
        File sample_dipcall_bed
        String chromosome
        
        String remote_dir
        
        Int bench_method
        Int min_sv_length
        Int max_sv_length
        
        File reference_fa
        File reference_fai
        File reference_agp
        File standard_chromosomes_bed
        File tandem_bed
        File not_tandem_bed
        
        File in_flag_shapeit4
        
        Int n_cpu = 2
        Int ram_size_gb = 4
        Int disk_size_gb = 20
    }
    parameter_meta {
        chromosome: "all=do not restrict to a chromosome."        
        min_sv_length: "The input VCFs (truvari, kanpig and dipcall) are first hard-filtered based on SVLEN, and fed to the chosen benchmarking tool."
        bench_method: "0=truvari bench with default parameters; 1=vcfdist."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        # Returns a BED file that excludes every gap from the AGP file of the
        # reference.
        #
        function GetReferenceGaps() {
            local INPUT_AGP=$1
            local OUTPUT_BED=$2
            
            awk 'BEGIN { FS="\t"; OFS="\t"; } { if ($1=="chr1" || $1=="chr2" || $1=="chr3" || $1=="chr4" || $1=="chr5" || $1=="chr6" || $1=="chr7" || $1=="chr8" || $1=="chr9" || $1=="chr10" || $1=="chr11" || $1=="chr12" || $1=="chr13" || $1=="chr14" || $1=="chr15" || $1=="chr16" || $1=="chr17" || $1=="chr18" || $1=="chr19" || $1=="chr20" || $1=="chr21" || $1=="chr22" || $1=="chrX" || $1=="chrY" || $1=="chrM") print $0 }' ${INPUT_AGP} > in.bed
            awk 'BEGIN { FS="\t"; OFS="\t"; } { if ($5=="N") print $0 }' in.bed > out.bed
            mv out.bed in.bed
            bedtools sort -i in.bed -faidx ~{reference_fai} > out.bed
            mv out.bed in.bed
            bedtools complement -i in.bed -g ~{reference_fai} > out.bed
            mv out.bed in.bed
            
            mv in.bed ${OUTPUT_BED}
        }
        
        
        # Puts in canonical form a raw VCF from dipcall. This is identical to
        # `SV_Integration_BuildTrainingResource.wdl`.
        #
        function CanonizeDipcallVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local INPUT_TBI=$3
            local MIN_SV_LENGTH=$4
            local MAX_SV_LENGTH=$5
            local STANDARD_CHROMOSOMES_BED=$6
            local NOT_GAPS_BED=$7
            
            
            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_in.vcf.gz
            mv ${INPUT_TBI} ${SAMPLE_ID}_in.vcf.gz.tbi
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Removing SNVs, records that are not marked as present, records
            # with a FILTER, and records with unresolved REF/ALT.
            ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || COUNT(GT="alt")=0 || (FILTER!="PASS" && FILTER!=".") || REF="*" || ALT="*"' --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the standard chromosomes
            ${TIME_COMMAND} bcftools filter --regions-file ${STANDARD_CHROMOSOMES_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Removing records in reference gaps
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the dipcall BED
            ${TIME_COMMAND} bcftools filter --regions-file ~{sample_dipcall_bed} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Making sure SVLEN and SVTYPE are consistently annotated
            truvari anno svinfo --minsize 1 ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the given length range
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_truth.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_truth.vcf.gz.tbi
        }

        
        #
        function Benchmark() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local OUTPUT_PREFIX=$3
                
            if [ ${INPUT_VCF#*.} -eq vcf.gz ]; then
                cp ${INPUT_VCF} ${OUTPUT_PREFIX}_input.vcf.gz
                cp ${INPUT_VCF}.tbi ${OUTPUT_PREFIX}_input.vcf.gz.tbi
            else
                bcftools view --output-type z ${INPUT_VCF} > ${OUTPUT_PREFIX}_input.vcf.gz
                tabix ${OUTPUT_PREFIX}_input.vcf.gz
            fi
            
            # Extracting calls with POS inside and outside TRs
            ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz > ${OUTPUT_PREFIX}_tr.vcf.gz &
            ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz > ${OUTPUT_PREFIX}_not_tr.vcf.gz &
            wait
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_tr.vcf.gz
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_not_tr.vcf.gz
        
            # Benchmarking
            if [ ~{bench_method} -eq 0 ]; then
                rm -rf ./${OUTPUT_PREFIX}_truvari_*
                ${TIME_COMMAND} truvari bench --includebed ~{sample_dipcall_bed} -b ${SAMPLE_ID}_truth.vcf.gz -c ${OUTPUT_PREFIX}_input.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_all/ &
                ${TIME_COMMAND} truvari bench --includebed ~{sample_dipcall_bed} -b ${SAMPLE_ID}_truth_tr.vcf.gz -c ${OUTPUT_PREFIX}_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_tr/ &
                ${TIME_COMMAND} truvari bench --includebed ~{sample_dipcall_bed} -b ${SAMPLE_ID}_truth_not_tr.vcf.gz -c ${OUTPUT_PREFIX}_not_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_not_tr/ &
                wait
                mv ./${OUTPUT_PREFIX}_truvari_all/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_all.txt
                mv ./${OUTPUT_PREFIX}_truvari_tr/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_tr.txt
                mv ./${OUTPUT_PREFIX}_truvari_not_tr/summary.json ./${SAMPLE_ID}_${OUTPUT_PREFIX}_not_tr.txt
            else
                # Assumed to be the only job running, so using all cores.
                rm -f ./${OUTPUT_PREFIX}_vcfdist_*
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_input.vcf.gz ${SAMPLE_ID}_truth.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_all/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_tr.vcf.gz ${SAMPLE_ID}_truth_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_tr/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_not_tr.vcf.gz ${SAMPLE_ID}_truth_not_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_not_tr/
                mv ./${OUTPUT_PREFIX}_vcfdist_all/precision-recall-summary.tsv ./${SAMPLE_ID}_${OUTPUT_PREFIX}_all.txt
                mv ./${OUTPUT_PREFIX}_vcfdist_tr/precision-recall-summary.tsv ./${SAMPLE_ID}_${OUTPUT_PREFIX}_tr.txt
                mv ./${OUTPUT_PREFIX}_vcfdist_not_tr/precision-recall-summary.tsv ./${SAMPLE_ID}_${OUTPUT_PREFIX}_not_tr.txt
            fi
            
            # Removing temporary files
            rm -f ${OUTPUT_PREFIX}_input.vcf.gz*
        }


        # --------------------------- Main program -----------------------------
        
        FILTER_STRING_TRUVARI="--sizemin ~{min_sv_length} --sizefilt ~{min_sv_length} --sizemax ~{max_sv_length}"
        FILTER_STRING_VCFDIST="--sv-threshold ~{min_sv_length} --largest-variant ~{max_sv_length}"
        # See https://github.com/TimD1/vcfdist/wiki/02-Parameters-and-Usage
        # Remark: `--max-supercluster-size` has to be >= `--largest-variant + 2`
        #         and we set it to 10002 to mimic kanpig's `--sizemax`.
        # Remark: we choose `--cluster gap` since it is faster. We choose 500
        #         to mimic kanpig inter-sample's `--neighdist` (the intra-
        #         sample value would be 1000, which might be too big).
        SV_STRING_VCFDIST="--cluster gap 500 --max-supercluster-size $((~{max_sv_length}+2)) --realign-query --realign-truth"
        
        
        GetReferenceGaps ~{reference_agp} not_gaps.bed
        
        # Canonizing the dipcall VCF
        tabix -f ~{sample_dipcall_vcf_gz}
        CanonizeDipcallVcf ~{sample_id} ~{sample_dipcall_vcf_gz} ~{sample_dipcall_vcf_gz}.tbi ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
        rm -f ~{sample_dipcall_vcf_gz}
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ~{sample_id}_truth.vcf.gz > ~{sample_id}_truth_tr.vcf.gz &
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ~{sample_id}_truth.vcf.gz > ~{sample_id}_truth_not_tr.vcf.gz &
        wait
        ${TIME_COMMAND} tabix -f ~{sample_id}_truth_tr.vcf.gz &
        ${TIME_COMMAND} tabix -f ~{sample_id}_truth_not_tr.vcf.gz &
        wait
        
        # Localizing the VCF to benchmark
        gsutil -m cp ~{remote_dir}/shapeit4/~{sample_id}_shapeit4.'bcf*' .
        
        # Ensuring the right header for VCFDIST
        bcftools reheader --fai ~{reference_fai} ~{sample_id}_shapeit4.bcf > out.bcf
        rm -f ~{sample_id}_shapeit4.bcf* ; mv out.bcf ~{sample_id}_shapeit4.bcf ; bcftools index ~{sample_id}_shapeit4.bcf
        
        # Subsetting every VCF to the given chromosome, if specified.
        if [ ~{chromosome} != "all" ]; then
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth.vcf.gz > ~{sample_id}_truth_prime.vcf.gz
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth_tr.vcf.gz > ~{sample_id}_truth_tr_prime.vcf.gz
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth_not_tr.vcf.gz > ~{sample_id}_truth_not_tr_prime.vcf.gz
            
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type b ~{sample_id}_shapeit4.bcf > ~{sample_id}_shapeit4_prime.bcf
            
            rm -f ~{sample_id}_truth.vcf.gz* ; mv ~{sample_id}_truth_prime.vcf.gz ~{sample_id}_truth.vcf.gz ; bcftools index ~{sample_id}_truth.vcf.gz
            rm -f ~{sample_id}_truth_tr.vcf.gz* ; mv ~{sample_id}_truth_tr_prime.vcf.gz ~{sample_id}_truth_tr.vcf.gz ; bcftools index ~{sample_id}_truth_tr.vcf.gz
            rm -f ~{sample_id}_truth_not_tr.vcf.gz* ; mv ~{sample_id}_truth_not_tr_prime.vcf.gz ~{sample_id}_truth_not_tr.vcf.gz ; bcftools index ~{sample_id}_truth_not_tr.vcf.gz
            
            rm -f ~{sample_id}_shapeit4.bcf* ; mv ~{sample_id}_shapeit4_prime.bcf ~{sample_id}_shapeit4.bcf ; bcftools index ~{sample_id}_shapeit4.bcf
        fi
        
        # Keeping only records in the given length range
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ~{sample_id}_shapeit4.bcf > out.vcf.gz
        rm -f ~{sample_id}_shapeit4.bcf* ; mv out.vcf.gz ~{sample_id}_shapeit4.vcf.gz ; tabix -f ~{sample_id}_shapeit4.vcf.gz
        
        # Keeping only records that are genotyped as present. This is important,
        # since truvari bench does not consider GTs when matching.
        ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type z ~{sample_id}_shapeit4.vcf.gz > out.vcf.gz
        rm -f ~{sample_id}_shapeit4.vcf.gz* ; mv out.vcf.gz ~{sample_id}_shapeit4.vcf.gz ; tabix -f ~{sample_id}_shapeit4.vcf.gz
        
        # Benchmarking
        Benchmark ~{sample_id} ~{sample_id}_shapeit4.vcf.gz shapeit4_~{min_sv_length}bp
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp '*_shapeit4_*.txt' ~{remote_dir}/shapeit4/precision_recall/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading truvari benchmarks. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Remark: `bcftools +trio-dnm2` marks the following triplets (child, father,
# mother) as de novos:
#
# 0|0, 0|0, 1|1
# 0|0, 0|1, 1|1
# 1|1, 1|1, 0|0
# 1|1, 0|0, 1|0
#
# I.e. the disappearance of a record is considered a de novo, and the
# appearance of a record on a haplotype is considered a de novo.
#
# Performance with 8 cores and 16GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# bcftools merge            700%        50M     3s
# bcftools +mendelian2      100%        20M     2s
# bcftools +trio-dnm2       100%        20M     10s
#
task BenchTrio {
    input {
        File ped_tsv
        Int ped_tsv_row
        
        String remote_indir
        String remote_outdir
        
        Int min_sv_length
        Int max_sv_length
        
        File tandem_bed
        File not_tandem_bed
        File autosomes_bed
        
        Array[File] in_flag
        
        Int n_cpu = 8
        Int ram_size_gb = 8
        Int disk_size_gb = 20
    }
    parameter_meta {
        ped_tsv_row: "The row (one-based) in `ped_tsv` that corresponds to this trio."
        min_sv_length: "The input VCFs are first hard-filtered based on SVLEN, and then fed to the chosen benchmarking tool."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        # @param 1: a VCF that contains only 3 samples (child, parents).
        # 
        function Benchmark() {
            local INPUT_VCF_GZ=$1
            local SAMPLE_ID=$2
            local SUFFIX=$3
            
            # 1. Mendelian error
            ${TIME_COMMAND} bcftools +mendelian2 ${INPUT_VCF_GZ} --ped ped.tsv > ${SAMPLE_ID}_mendelian_${SUFFIX}.txt
            
            # 2. De novo rate
            
            # 2.1 +trio-dnm2
            ${TIME_COMMAND} bcftools +trio-dnm2 --use-NAIVE --chrX GRCh38 --ped ped.tsv --output-type z ${INPUT_VCF_GZ} > ${SAMPLE_ID}_annotated.vcf.gz
            NUMERATOR=$( bcftools view --no-header --include 'FORMAT/DNM[0]=1' ${SAMPLE_ID}_annotated.vcf.gz | wc -l )
            DENOMINATOR=$( bcftools view --no-header --include 'GT[0]="alt" && COUNT(GT="mis")=0' ${SAMPLE_ID}_annotated.vcf.gz | wc -l )
            echo -e "${NUMERATOR},${DENOMINATOR}" > ${SAMPLE_ID}_dnm1_${SUFFIX}.txt
            rm -f ${SAMPLE_ID}_annotated.vcf.gz*
            
            # 2.2 Simple count (and saving the whole matrix for future analysis)
            ${TIME_COMMAND} truvari anno numneigh --sizemin 1 --refdist 1000 ${INPUT_VCF_GZ} | bgzip > ${SAMPLE_ID}_annotated.vcf.gz
            ${TIME_COMMAND} bcftools query --format '[%GT,]0,0,0,0,0,%INFO/SVTYPE,%INFO/SVLEN,%INFO/NumNeighbors\n' ${SAMPLE_ID}_annotated.vcf.gz > ${SAMPLE_ID}_matrix.txt
            ${TIME_COMMAND} java -cp ~{docker_dir} CountDeNovoSimple ${SAMPLE_ID}_matrix.txt > ${SAMPLE_ID}_dnm2_${SUFFIX}.txt
            rm -f ${SAMPLE_ID}_annotated.vcf.gz*
        }
        

        # --------------------------- Main program -----------------------------
        
        head -n ~{ped_tsv_row} ~{ped_tsv} | tail -n 1 > ped.tsv
        PROBAND_ID=$(cut -f 2 ped.tsv)
        FATHER_ID=$(cut -f 3 ped.tsv)
        MOTHER_ID=$(cut -f 4 ped.tsv)
        
        # Localizing. These files may contain 0/0 records.
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir}/${PROBAND_ID}_'*.vcf.gz*' ~{remote_indir}/${FATHER_ID}_'*.vcf.gz*' ~{remote_indir}/${MOTHER_ID}_'*.vcf.gz*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                TEST=$(gsutil -m cp ~{remote_indir}/${PROBAND_ID}_'*.bcf*' ~{remote_indir}/${FATHER_ID}_'*.bcf*' ~{remote_indir}/${MOTHER_ID}_'*.bcf*' . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading VCF files. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            else 
                break
            fi
        done
        
        # Restricting to autosomes. This is just for simplicity and is not
        # strictly necessary.
        TEST=$(ls *.vcf.gz && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${PROBAND_ID}_*.bcf) > ${PROBAND_ID}_autosomes.vcf.gz &
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${FATHER_ID}_*.bcf) > ${FATHER_ID}_autosomes.vcf.gz &
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${MOTHER_ID}_*.bcf) > ${MOTHER_ID}_autosomes.vcf.gz &
        else
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${PROBAND_ID}_*.vcf.gz) > ${PROBAND_ID}_autosomes.vcf.gz &
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${FATHER_ID}_*.vcf.gz) > ${FATHER_ID}_autosomes.vcf.gz &
            ${TIME_COMMAND} bcftools filter --regions-file ~{autosomes_bed} --regions-overlap pos --output-type z $(ls ${MOTHER_ID}_*.vcf.gz) > ${MOTHER_ID}_autosomes.vcf.gz &
        fi
        wait
        tabix -f ${PROBAND_ID}_autosomes.vcf.gz &
        tabix -f ${FATHER_ID}_autosomes.vcf.gz &
        tabix -f ${MOTHER_ID}_autosomes.vcf.gz &
        wait
        echo ${PROBAND_ID}_autosomes.vcf.gz > list.txt
        echo ${FATHER_ID}_autosomes.vcf.gz >> list.txt
        echo ${MOTHER_ID}_autosomes.vcf.gz >> list.txt
        ls -laht
        
        # Remark: we keep all records from the three files, since otherwise
        # bcftools merge might transform a 0/0 (that would have disappeared)
        # into a ./., affecting the number of records over which Mendelian
        # error is computed.
        
        # Keeping only records in the given length range
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ${PROBAND_ID}_autosomes.vcf.gz > out.vcf.gz
        rm -f ${PROBAND_ID}_autosomes.vcf.gz* ; mv out.vcf.gz ${PROBAND_ID}_autosomes.vcf.gz ; tabix -f ${PROBAND_ID}_autosomes.vcf.gz
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ${FATHER_ID}_autosomes.vcf.gz > out.vcf.gz
        rm -f ${FATHER_ID}_autosomes.vcf.gz* ; mv out.vcf.gz ${FATHER_ID}_autosomes.vcf.gz ; tabix -f ${FATHER_ID}_autosomes.vcf.gz
        ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type z ${MOTHER_ID}_autosomes.vcf.gz > out.vcf.gz
        rm -f ${MOTHER_ID}_autosomes.vcf.gz* ; mv out.vcf.gz ${MOTHER_ID}_autosomes.vcf.gz ; tabix -f ${MOTHER_ID}_autosomes.vcf.gz
        
        # Merging records by ID, since the records in every VCF originate from
        # the same cohort VCF, which had distinct IDs.
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge id --output-type z --file-list list.txt > trio.vcf.gz
        ${TIME_COMMAND} tabix -f trio.vcf.gz
        rm -f ${PROBAND_ID}_* ${FATHER_ID}_* ${MOTHER_ID}_*
        ls -laht
        
        # Benchmarking
        Benchmark trio.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_all
        ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z trio.vcf.gz > tr.vcf.gz
        ${TIME_COMMAND} tabix -f tr.vcf.gz
        Benchmark tr.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_tr
        rm -f tr.vcf.gz*
        ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z trio.vcf.gz > not_tr.vcf.gz
        ${TIME_COMMAND} tabix -f not_tr.vcf.gz
        Benchmark not_tr.vcf.gz ${PROBAND_ID} ~{min_sv_length}bp_not_tr
        rm -f not_tr.vcf.gz*
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp '*.txt' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading stats. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
