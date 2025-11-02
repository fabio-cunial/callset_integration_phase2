version 1.0


# 
#
workflow RegenotypingAnalysis {
    input {
        File cohort_truvari_vcf_gz
        File cohort_truvari_tbi
        Array[Int] min_n_samples = [2, 3, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
        Int min_sv_length = 20
        
        
        String infrequent_remote_dir
        
        Array[String] precision_recall_samples
        Array[String] precision_recall_sex
        Array[File] precision_recall_bam
        Array[File] precision_recall_bai
        String precision_recall_remote_dir
        
        
        Array[File] precision_recall_samples_dipcall_vcf_gz
        Array[File] precision_recall_samples_dipcall_bed
        Int precision_recall_bench_method
        
        File mendelian_error_ped_tsv
        Int mendelian_error_n_trios

        
        File reference_fa
        File reference_fai
        File tandem_bed
        File ploidy_bed_male
        File ploidy_bed_female
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
    call PrepareCohortBcf {
        input:
            cohort_truvari_vcf_gz = cohort_truvari_vcf_gz,
            cohort_truvari_tbi = cohort_truvari_tbi
    }
    call FilterCohortBcf_ByLength {
        input:
            cohort_truvari_bcf = PrepareCohortBcf.out_bcf,
            cohort_truvari_csi = PrepareCohortBcf.out_csi
    }
    scatter (i in range(length(min_n_samples))) {
        call PartitionCohortBcf {
            input:
                cohort_truvari_bcf = FilterCohortBcf_ByLength.out_bcf,
                cohort_truvari_csi = FilterCohortBcf_ByLength.out_csi,
                min_n_samples = min_n_samples[i]
        }
        call SplitInfrequentBcf {
            input:
                infrequent_cohort_bcf = PartitionCohortBcf.infrequent_bcf,
                infrequent_cohort_csi = PartitionCohortBcf.infrequent_csi,
                samples = precision_recall_samples,
                remote_outdir = infrequent_remote_dir
        }
        scatter (i in range(length(precision_recall_samples))) {
            call BuildPersonalizedVcf {
                input:
                    sample_id = precision_recall_samples[i],
                    frequent_cohort_bcf = PartitionCohortBcf.frequent_bcf,
                    frequent_cohort_csi = PartitionCohortBcf.frequent_csi,
                    infrequent_remote_dir = infrequent_remote_dir,
                    in_flag = SplitInfrequentBcf.out_flag
            }
            call Kanpig {
                input:
                    sample_id = precision_recall_samples[i],
                    sex = precision_recall_sex[i],
                    personalized_vcf_gz = BuildPersonalizedVcf.out_vcf_gz,
                    personalized_tbi = BuildPersonalizedVcf.out_tbi,
                    alignments_bam = precision_recall_bam[i],
                    alignments_bai = precision_recall_bai[i],
                    remote_dir = precision_recall_remote_dir,
                    reference_fa = reference_fa,
                    reference_fai = reference_fai,
                    ploidy_bed_male = ploidy_bed_male,
                    ploidy_bed_female = ploidy_bed_female
            }
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
                    
                    in_flag = Kanpig.out_flag
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
    
    Int disk_size_gb = 4*ceil(size(cohort_truvari_bcf,"GB"))
    
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
        rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index out.bcf
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
        rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index out.bcf
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


# Writes to separate files the records that are present in each sample column of
# the infrequent cohort BCF.
#
task SplitInfrequentBcf {
    input {
        File infrequent_cohort_bcf
        File infrequent_cohort_csi
        
        Array[String] samples
        String remote_outdir
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        infrequent_cohort_bcf: "Assumed to have all the sample columns in the original truvari collapse VCF."
        remote_outdir: "The result of the split is stored in this bucket location."
    }
    
    Int disk_size_gb = 4*ceil(size(infrequent_cohort_bcf,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        
        # Splitting
        echo ~{sep="," samples} | tr ',' '\n' > samples.txt
        ${TIME_COMMAND} bcftools +split --samples-file samples.txt --output-type b --output . ~{infrequent_cohort_bcf}
        rm -f ~{infrequent_cohort_bcf}
        
        # Keeping only present records
        for FILE in $(ls *.bcf); do
            SAMPLE_ID=$(basename ${FILE} .bcf)
            ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type b ${FILE} > ${SAMPLE_ID}_infrequent.bcf
            bcftools index ${SAMPLE_ID}_infrequent.bcf
        done
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp '*_infrequent.bcf*' ~{remote_outdir}/ && echo 0 || echo 1)
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
        String infrequent_remote_dir
        
        File in_flag
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        frequent_cohort_bcf: "Assumed to have a single sample column."
        infrequent_remote_indir: "Remote directory containing the infrequent cohort records that occur in each sample (with names `SAMPLEID_infrequent.bcf`)."
        flag: "Just to flag to the sheduler that this task can start."
    }
    
    Int disk_size_gb = 4*ceil(size(frequent_cohort_bcf,"GB"))
    
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
        
        # Localizing infrequent cohort records that occur in this sample
        while : ; do
            TEST=$(gsutil -m cp ~{infrequent_remote_dir}/~{sample_id}_infrequent.'bcf*' . && echo 0 || echo 1)
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
        
        String remote_dir
        
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
            TEST=$(gsutil -m cp ~{sample_id}_kanpig.'vcf*' ~{remote_dir}/personalized/ && echo 0 || echo 1)
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

















#----------------------- Benchmarking -------------------------


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
