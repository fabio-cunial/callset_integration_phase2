version 1.0


# Structure of `remote_output_dir`:
#
# ├── truvari/         for each sample, its records in the truvari collapse VCF;
# │   └── precision_recall/             for each sample, precision/recall stats;
# ├── X_samples
# │   ├── infrequent/             for each sample, the infrequent records in it;
# │   ├── kanpig/             for each sample, its personalized and re-GT'd VCF;
# │   └── precision_recall              for each sample, precision/recall stats;
#
workflow SV_Integration_RegenotypingAnalysis {
    input {
        File cohort_truvari_vcf_gz
        File cohort_truvari_tbi
        Array[Int] min_n_samples = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
        Int min_sv_length = 20
        Int max_sv_length = 10000
        
        String remote_dir
        
        Int precision_recall_bench_method
        Array[String] precision_recall_samples
        Array[String] precision_recall_sex
        Array[File] precision_recall_bam
        Array[File] precision_recall_bai
        Array[File] precision_recall_samples_dipcall_vcf_gz
        Array[File] precision_recall_samples_dipcall_bed
        
        File reference_fa
        File reference_fai
        File reference_agp
        File standard_chromosomes_bed
        File tandem_bed
        File ploidy_bed_male
        File ploidy_bed_female
    }
    parameter_meta {
        precision_recall_bench_method: "0=truvari bench, 1=vcfdist."
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    call PrepareCohortBcf {
        input:
            cohort_truvari_vcf_gz = cohort_truvari_vcf_gz,
            cohort_truvari_tbi = cohort_truvari_tbi
    }
    call FilterCohortBcf_ByLength {
        input:
            cohort_truvari_bcf = PrepareCohortBcf.out_bcf,
            cohort_truvari_csi = PrepareCohortBcf.out_csi,
            min_sv_length = min_sv_length
    }
    call SplitBcf as split_truvari {
        input:
            cohort_bcf = FilterCohortBcf_ByLength.out_bcf,
            cohort_csi = FilterCohortBcf_ByLength.out_csi,
            samples = precision_recall_samples,
            remote_dir = remote_dir+"/truvari",
            suffix = "truvari"
    }
    scatter (i in range(length(min_n_samples))) {
        call PartitionCohortBcf {
            input:
                cohort_truvari_bcf = FilterCohortBcf_ByLength.out_bcf,
                cohort_truvari_csi = FilterCohortBcf_ByLength.out_csi,
                min_n_samples = min_n_samples[i]
        }
        call SplitBcf as split_infrequent {
            input:
                cohort_bcf = PartitionCohortBcf.infrequent_bcf,
                cohort_csi = PartitionCohortBcf.infrequent_csi,
                samples = precision_recall_samples,
                remote_dir = remote_dir+"/"+min_n_samples[i]+"_samples/infrequent",
                suffix = "infrequent"
        }
        scatter (i in range(length(precision_recall_samples))) {
            call BuildPersonalizedVcf {
                input:
                    sample_id = precision_recall_samples[i],
                    frequent_cohort_bcf = PartitionCohortBcf.frequent_bcf,
                    frequent_cohort_csi = PartitionCohortBcf.frequent_csi,
                    remote_indir = remote_dir+"/"+min_n_samples[i]+"_samples/infrequent",
                    in_flag = split_infrequent.out_flag
            }
            call Kanpig {
                input:
                    sample_id = precision_recall_samples[i],
                    sex = precision_recall_sex[i],
                    personalized_vcf_gz = BuildPersonalizedVcf.out_vcf_gz,
                    personalized_tbi = BuildPersonalizedVcf.out_tbi,
                    alignments_bam = precision_recall_bam[i],
                    alignments_bai = precision_recall_bai[i],
                    remote_outdir = remote_dir+"/"+min_n_samples[i]+"_samples/kanpig",
                    reference_fa = reference_fa,
                    reference_fai = reference_fai,
                    ploidy_bed_male = ploidy_bed_male,
                    ploidy_bed_female = ploidy_bed_female
            }
            call PrecisionRecallAnalysis {
                input:
                    sample_id = precision_recall_samples[i],
                    sample_dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[i],
                    sample_dipcall_bed = precision_recall_samples_dipcall_bed[i],
                    
                    remote_dir = remote_dir,
                    min_n_samples = min_n_samples[i],
                
                    bench_method = precision_recall_bench_method,
                    min_sv_length = min_sv_length,
                    max_sv_length = max_sv_length,    
                
                    reference_fa = reference_fa,
                    reference_fai = reference_fai,
                    reference_agp = reference_agp,
                    standard_chromosomes_bed = standard_chromosomes_bed,
                    tandem_bed = ComplementBed.sorted_bed,
                    not_tandem_bed = ComplementBed.complement_bed,
                    
                    in_flag_truvari = split_truvari.out_flag,
                    in_flag_kanpig = Kanpig.out_flag
            }
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


# Keeps all sample columns in the output, but adds an INFO field that counts the
# number of samples each record occurs in, and enforces a distinct ID in every
# record.
#
# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, SSD:
#
# TOOL                CPU     RAM     TIME
#
#
task PrepareCohortBcf {
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
        
        # Converting to BCF, to speed up all steps downstream.
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --output-type b in.vcf.gz > out.bcf
        rm -f in.vcf.gz* ; mv out.bcf in.bcf ; bcftools index out.bcf
        
        # Annotating every record with the number of samples it occurs in.
        ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%COUNT(GT="alt")\n' in.bcf | bgzip -c > annotations.tsv.gz
        tabix -s1 -b2 -e2 annotations.tsv.gz
        echo '##INFO=<ID=N_SAMPLES,Number=1,Type=Integer,Description="Number of samples where the record was discovered">' > header.txt
        ${TIME_COMMAND} bcftools annotate --header-lines header.txt --annotations annotations.tsv.gz --columns CHROM,POS,~ID,REF,ALT,N_SAMPLES --output-type z in.bcf > out.bcf
        rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index in.bcf
        rm -f annotations.tsv.gz
        
        # Enforcing a distinct ID in every record
        bcftools view --header-only in.bcf > header.txt
        N_ROWS=$(wc -l < header.txt)
        date
        (  head -n $(( ${N_ROWS} - 1 )) header.txt ; \
           echo '##INFO=<ID=ORIGINAL_ID,Number=1,Type=String,Description="Original ID from truvari collapse">' ; \
           echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" ; \
           bcftools view --no-header in.bcf | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { \
               $8=printf("%s;ORIGINAL_ID=%s,$8,gsub(/;/,"_",$3)); \
               $3=printf("%d",++i); \
               print $0 \
           }' \
        ) | bcftools view --output-type b > out.bcf
        date
        rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index in.bcf
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


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, SSD:
#
# TOOL                CPU     RAM     TIME
#
#
task FilterCohortBcf_ByLength {
    input {
        File cohort_truvari_bcf
        File cohort_truvari_csi
        
        Int min_sv_length
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        cohort_truvari_bcf: "Every record is assumed be already annotated with the correct SVLEN."
    }
    
    Int disk_size_gb = 4*ceil(size(cohort_truvari_bcf,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'ABS(SVLEN)>='~{min_sv_length} --output-type b ~{cohort_truvari_bcf} > out.bcf
        ${TIME_COMMAND} bcftools index out.bcf
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


# Partitions `cohort_truvari_bcf` into the following BCFs:
# - The subset of all records that occur in `>= min_n_samples`. This file has
#   just one artificial sample column set to 0/1.
# - The subset of all records that occur in `< min_n_samples`. This file has
#   all the original sample columns kept intact.
#
# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, SSD:
#
# TOOL                CPU     RAM     TIME
# 
#
task PartitionCohortBcf {
    input {
        File cohort_truvari_bcf
        File cohort_truvari_csi
        
        Int min_n_samples
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 4*ceil(size(cohort_truvari_bcf,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Splitting
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'N_SAMPLES>='~{min_n_samples} --output-type b ~{cohort_truvari_bcf} > tmp.bcf &
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'N_SAMPLES<'~{min_n_samples} --output-type b ~{cohort_truvari_bcf} > infrequent_~{min_n_samples}.bcf &
        wait
        ${TIME_COMMAND} bcftools index tmp.bcf &
        ${TIME_COMMAND} bcftools index infrequent_~{min_n_samples}.bcf &
        wait
        
        # Enforcing a single sample in the frequent BCF
        bcftools view --header-only tmp.vcf.gz > header.txt
        N_ROWS=$(wc -l < header.txt)
        date
        (  head -n $(( ${N_ROWS} - 1 )) header.txt ; \
           echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" ; \
           bcftools view --no-header tmp.bcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\tGT\t0/1\n",$1,$2,$3,$4,$5,$6,$7,$8); }' \
        ) | bcftools view --output-type b > frequent_~{min_n_samples}.bcf
        bcftools index frequent_~{min_n_samples}.bcf
    >>>
    
    output {
        File frequent_bcf = "frequent_"+min_n_samples+".bcf"
        File frequent_csi = "frequent_"+min_n_samples+".bcf.csi"
        File infrequent_bcf = "infrequent_"+min_n_samples+".bcf"
        File infrequent_csi = "infrequent_"+min_n_samples+".bcf.csi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Writes to separate files all and only the records that are present in each
# sample column of a cohort BCF.
#
task SplitBcf {
    input {
        File cohort_bcf
        File cohort_csi
        
        Array[String] samples
        String remote_dir
        String suffix
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        cohort_bcf: "Assumed to have all the sample columns in the original truvari collapse VCF."
        remote_dir: "The result of the split is stored in this bucket location."
    }
    
    Int disk_size_gb = 4*ceil(size(cohort_bcf,"GB"))
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        
        # Splitting
        echo ~{sep="," samples} | tr ',' '\n' > samples.txt
        ${TIME_COMMAND} bcftools +split --samples-file samples.txt --output-type b --output . ~{cohort_bcf}
        rm -f ~{cohort_bcf}
        
        # Keeping only present records
        for FILE in $(ls *.bcf); do
            SAMPLE_ID=$(basename ${FILE} .bcf)
            ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type b ${FILE} > ${SAMPLE_ID}_~{suffix}.bcf
            bcftools index ${SAMPLE_ID}_~{suffix}.bcf
        done
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp '*_'~{suffix}'.bcf*' ~{remote_dir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading infrequent VCFs. Trying again..."
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


# Merges all the frequent cohort records with the infrequent cohort records
# that occur in a given sample.
#
# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999:
#
# TOOL                CPU     RAM     TIME
# 
#
task BuildPersonalizedVcf {
    input {
        String sample_id
        File frequent_cohort_bcf
        File frequent_cohort_csi
        String remote_indir
        
        File in_flag
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        frequent_cohort_bcf: "Assumed to have a single sample column."
        in_flag: "Just to flag to the sheduler that this task can start safely."
    }
    
    Int disk_size_gb = 4*ceil(size(frequent_cohort_bcf,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        INFINITY="1000000000"
        
        # Localizing infrequent cohort records that occur in this sample
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir}/~{sample_id}_infrequent.'bcf*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading ~{sample_id}_infrequent.bcf. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        # Changing the header of the frequent BCF
        echo ~{sample_id} > samples.txt
        ${TIME_COMMAND} bcftools reheader --threads ${N_THREADS} --samples samples.txt ~{frequent_cohort_bcf} > frequent.bcf
        bcftools index frequent.bcf
        rm -f ~{frequent_cohort_bcf}
        
        # Merging infrequent and frequent BCFs
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --rm-dups exact --output-type z ~{sample_id}_infrequent.bcf frequent.bcf > ~{sample_id}_personalized.vcf.gz
        ${TIME_COMMAND} tabix -f ~{sample_id}_personalized.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = sample_id+"_personalized.vcf.gz"
        File out_tbi = sample_id+"_personalized.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}


#
task Kanpig {
    input {
        String sample_id
        String sex
        File personalized_vcf_gz
        File personalized_tbi
        File alignments_bam
        File alignments_bai
        
        String remote_outdir
        
        String kanpig_params_singlesample = "--neighdist 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"
        File reference_fa
        File reference_fai
        File ploidy_bed_male
        File ploidy_bed_female
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 4*ceil( size(alignments_bam,"GB") )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export RUST_BACKTRACE="full"
        INFINITY="1000000000"
        
        if [ ~{sex} == "M" ]; then
            PLOIDY_BED=$(echo ~{ploidy_bed_male})
        else
            PLOIDY_BED=$(echo ~{ploidy_bed_female})
        fi
        
        mv ~{personalized_vcf_gz} in.vcf.gz
        mv ~{personalized_tbi} in.vcf.gz.tbi
        
        # Remark: kanpig needs --sizemin >= --kmer
        ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_singlesample} --sizemin 10 --sizemax ${INFINITY} --reference ~{reference_fa} --input in.vcf.gz --reads ~{alignments_bam} --out out.vcf
        rm -f in.vcf.gz* ; mv out.vcf in.vcf
        
        # Sorting
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z in.vcf > out.vcf.gz
        rm -f in.vcf ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        # Discarding records that are not marked as present by kanpig
        ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type z in.vcf.gz > out.vcf.gz
        rm -f in.vcf.gz* ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        mv in.vcf.gz ~{sample_id}_kanpig.vcf.gz
        mv in.vcf.gz.tbi ~{sample_id}_kanpig.vcf.gz.tbi
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m cp ~{sample_id}_kanpig.'vcf*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading ~{sample_id}_kanpig.vcf.gz. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        
        echo "done" > out.txt
    >>>
    
    output {
        File out_flag = "out.txt"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}




#------------------------------- Benchmarking ----------------------------------

# Remark: `truvari bench` is used with default parameters.
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
        File sample_dipcall_vcf_gz
        File sample_dipcall_bed
        
        String remote_dir
        String min_n_samples
        
        Int bench_method
        Int min_sv_length = 20
        Int max_sv_length = 10000
        
        File reference_fa
        File reference_fai
        File reference_agp
        File standard_chromosomes_bed
        File tandem_bed
        File not_tandem_bed
        
        File in_flag_truvari
        File in_flag_kanpig
        
        Int n_cpu = 4
        Int ram_size_gb = 32
        Int disk_size_gb = 50
    }
    parameter_meta {
    }
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        
        
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
            ${TIME_COMMAND} bcftools filter --regions-file ${SAMPLE_ID}.bed --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
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
            local INPUT_VCF_GZ=$1
            local OUTPUT_PREFIX=$2
            
            cp ${INPUT_VCF_GZ} ${OUTPUT_PREFIX}_input.vcf.gz
            cp ${INPUT_VCF_GZ}.tbi ${OUTPUT_PREFIX}_input.vcf.gz.tbi
            
            # Extracting calls with POS inside and outside TRs
            ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz > ${OUTPUT_PREFIX}_tr.vcf.gz &
            ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${OUTPUT_PREFIX}_input.vcf.gz > ${OUTPUT_PREFIX}_not_tr.vcf.gz &
            wait
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_tr.vcf.gz
            ${TIME_COMMAND} tabix -f ${OUTPUT_PREFIX}_not_tr.vcf.gz
        
            # Benchmarking
            if [ ~{bench_method} -eq 0 ]; then
                rm -rf ./${OUTPUT_PREFIX}_truvari_*
                ${TIME_COMMAND} truvari bench --includebed ~{sample_dipcall_bed} -b truth.vcf.gz -c ${OUTPUT_PREFIX}_input.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_all/ &
                ${TIME_COMMAND} truvari bench --includebed ~{sample_dipcall_bed} -b truth_tr.vcf.gz -c ${OUTPUT_PREFIX}_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_tr/ &
                ${TIME_COMMAND} truvari bench --includebed ~{sample_dipcall_bed} -b truth_not_tr.vcf.gz -c ${OUTPUT_PREFIX}_not_tr.vcf.gz ${FILTER_STRING_TRUVARI} -o ./${OUTPUT_PREFIX}_truvari_not_tr/ &
                wait
                mv ./${OUTPUT_PREFIX}_truvari_all/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_all.txt
                mv ./${OUTPUT_PREFIX}_truvari_tr/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_tr.txt
                mv ./${OUTPUT_PREFIX}_truvari_not_tr/summary.json ./~{sample_id}_${OUTPUT_PREFIX}_not_tr.txt
            else
                # Assumed to be the only job running, so using all cores.
                rm -f ./${OUTPUT_PREFIX}_vcfdist_*
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_input.vcf.gz truth.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_all/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_tr.vcf.gz truth_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_tr/
                ${TIME_COMMAND} vcfdist ${OUTPUT_PREFIX}_not_tr.vcf.gz truth_not_tr.vcf.gz ~{reference_fa} --max-threads ${N_THREADS} --max-ram $(( ~{ram_size_gb} - 2 )) ${SV_STRING_VCFDIST} ${FILTER_STRING_VCFDIST} --bed ~{sample_dipcall_bed} --prefix ./${OUTPUT_PREFIX}_vcfdist_not_tr/
                mv ./${OUTPUT_PREFIX}_vcfdist_all/precision-recall-summary.tsv ./~{sample_id}_${OUTPUT_PREFIX}_all.txt
                mv ./${OUTPUT_PREFIX}_vcfdist_tr/precision-recall-summary.tsv ./~{sample_id}_${OUTPUT_PREFIX}_tr.txt
                mv ./${OUTPUT_PREFIX}_vcfdist_not_tr/precision-recall-summary.tsv ./~{sample_id}_${OUTPUT_PREFIX}_not_tr.txt
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
        
        # Localizing the VCFs to benchmark
        gsutil -m cp ~{remote_dir}/truvari/~{sample_id}_truvari.vcf.'gz*' .
        gsutil -m cp ~{remote_dir}/~{min_n_samples}_samples/kanpig/~{sample_id}_kanpig.vcf.'gz*' .
        
        # Benchmarking
        Benchmark ~{sample_id}_truvari.vcf.gz ~{sample_id}_truvari
        Benchmark ~{sample_id}_kanpig.vcf.gz ~{sample_id}_kanpig
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m cp '*_truvari_*.txt' ~{remote_dir}/truvari/precision_recall/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading truvari benchmarks. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil -m cp '*_kanpig_*.txt' ~{remote_dir}/~{min_n_samples}_samples/precision_recall/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading kanpig benchmarks. Trying again..."
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
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}







# Performance with 8 cores and 16GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# bcftools filter           100%        2G      1m
# bcftools merge            200%        50M     1m
# truvari collapse          
# bcftools +mendelian2      100%        50M     1m
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
        
        Int truvari_collapse
        
        File tandem_bed
        File not_tandem_bed
        
        Int n_cpu = 8
        Int ram_size_gb = 16
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
        
        # Remark: this function overwrites the input file with its collapse.
        #
        function truvari_collapse() {
            local INPUT_VCF_GZ=$1
            local ID=$2
            
            truvari collapse --input ${INPUT_VCF_GZ} --sizemin 0 --sizemax 1000000 --keep common --gt all | bcftools sort --max-mem $(( ~{ram_size_gb} - 4 ))G --output-type z > collapse_${ID}.vcf.gz
            tabix -f collapse_${ID}.vcf.gz
            mv collapse_${ID}.vcf.gz ${INPUT_VCF_GZ}
            mv collapse_${ID}.vcf.gz.tbi ${INPUT_VCF_GZ}.tbi
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
        if [ ~{truvari_collapse} -ne 0 ]; then
            truvari_collapse trio_kanpig.vcf.gz trio_kanpig
        fi
        
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
        if [ ~{truvari_collapse} -ne 0 ]; then
            truvari_collapse trio_07.vcf.gz trio_07 &
            truvari_collapse trio_09.vcf.gz trio_09 &
            wait
        fi
        
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
