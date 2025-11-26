version 1.0


# 
#
workflow Trgt2Kanpig {
    input {
        File sv_integration_chunk_tsv
        Array[String] precision_recall_samples
        Array[File] precision_recall_samples_dipcall_vcf_gz
        Array[File] precision_recall_samples_dipcall_bed
        
        String region = "chr6"
        String remote_outdir
        
        Int min_sv_length = 20
        Int max_sv_length = 10000
        String kanpig_params_singlesample = "--neighdist 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"
        Int bench_method = 1
        
        File reference_fa
        File reference_fai
        File standard_chromosomes_bed
        File reference_agp
        File tandem_bed
        File ploidy_bed_female
        File ploidy_bed_male
    }
    parameter_meta {
        sv_integration_chunk_tsv: "A subset of the rows of table `sv_integration_hg38`, without the header."
        region: "Only consider VCF records in this genomic region. Set to 'all' to disable."
        max_sv_length: "Calls above this length are deemed 'ultralong', are not given to kanpig re-genotyping, and are processed separately."
        remote_outdir: "Without final slash. Where the output of intra-sample truvari and kanpig is stored for each sample."
    }
    
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai
    }
    scatter (i in range(length(precision_recall_samples))) {
        call IntraSampleKanpig as kanpig_no_trgt {
            input:
                sv_integration_chunk_tsv = sv_integration_chunk_tsv,
                selected_row = i+1,
                use_trgt = 0,
                
                region = region,
                min_sv_length = min_sv_length,
                max_sv_length = max_sv_length,
                kanpig_params_singlesample = kanpig_params_singlesample,
                            
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                standard_chromosomes_bed = standard_chromosomes_bed,
                reference_agp = reference_agp,
                ploidy_bed_female = ploidy_bed_female,
                ploidy_bed_male = ploidy_bed_male
        }
        call PrecisionRecallAnalysis as bench_no_trgt {
            input:
                sample_id = precision_recall_samples[i],
                kanpig_vcf_gz = kanpig_no_trgt.kanpig_vcf_gz,
                kanpig_tbi = kanpig_no_trgt.kanpig_tbi,
                experiment_id = "no_trgt",
                sample_dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[i],
                sample_dipcall_bed = precision_recall_samples_dipcall_bed[i],
                
                chromosome = region,
                remote_outdir = remote_outdir,
                
                bench_method = bench_method,
                min_sv_length = min_sv_length,
                max_sv_length = max_sv_length,
            
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_agp = reference_agp,
                standard_chromosomes_bed = standard_chromosomes_bed,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed
        }
        
        call IntraSampleKanpig as kanpig_trgt {
            input:
                sv_integration_chunk_tsv = sv_integration_chunk_tsv,
                selected_row = i+1,
                use_trgt = 1,
                
                region = region,
                min_sv_length = min_sv_length,
                max_sv_length = max_sv_length,
                kanpig_params_singlesample = kanpig_params_singlesample,
                            
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                standard_chromosomes_bed = standard_chromosomes_bed,
                reference_agp = reference_agp,
                ploidy_bed_female = ploidy_bed_female,
                ploidy_bed_male = ploidy_bed_male
        }
        call PrecisionRecallAnalysis as bench_trgt {
            input:
                sample_id = precision_recall_samples[i],
                kanpig_vcf_gz = kanpig_trgt.kanpig_vcf_gz,
                kanpig_tbi = kanpig_trgt.kanpig_tbi,
                experiment_id = "trgt",
                sample_dipcall_vcf_gz = precision_recall_samples_dipcall_vcf_gz[i],
                sample_dipcall_bed = precision_recall_samples_dipcall_bed[i],
                
                chromosome = region,
                remote_outdir = remote_outdir,
                
                bench_method = bench_method,
                min_sv_length = min_sv_length,
                max_sv_length = max_sv_length,
            
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_agp = reference_agp,
                standard_chromosomes_bed = standard_chromosomes_bed,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed
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


#
task IntraSampleKanpig {
    input {
        File sv_integration_chunk_tsv
        Int selected_row
        Int use_trgt
        
        String region
        Int min_sv_length
        Int max_sv_length
        String kanpig_params_singlesample
        
        File reference_fa
        File reference_fai
        File standard_chromosomes_bed
        File reference_agp
        File ploidy_bed_female
        File ploidy_bed_male
        
        Int n_cpu = 6
        Int ram_size_gb = 8
        Int disk_size_gb = 256
    }
    parameter_meta {
        selected_row: "(One-based) A row of `sv_integration_chunk_tsv`."
    }
    
    String docker_dir = "/callset_integration"
    Int compression_level = 1
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        # @param 
        # $2 1=Localizes everything except the BAM. 2=Localizes just the BAM.
        # $3 A row of `sv_integration_chunk_tsv`.
        #
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local MODE=$2
            local LINE=$3
            
            ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 3)
            ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 4)
            PAV_BED=$(echo ${LINE} | cut -d , -f 5)
            PAV_TBI=$(echo ${LINE} | cut -d , -f 6)
            PAV_VCF_GZ=$(echo ${LINE} | cut -d , -f 7)
            PBSV_TBI=$(echo ${LINE} | cut -d , -f 8)
            PBSV_VCF_GZ=$(echo ${LINE} | cut -d , -f 9)
            SNIFFLES_TBI=$(echo ${LINE} | cut -d , -f 10)
            SNIFFLES_VCF_GZ=$(echo ${LINE} | cut -d , -f 11)
            TRGT_TBI=$(echo ${LINE} | cut -d , -f 12)
            TRGT_VCF_GZ=$(echo ${LINE} | cut -d , -f 13)
            
            if [ ${MODE} -eq 2 ]; then
                while : ; do
                    TEST=$(gsutil -m cp ${ALIGNED_BAM} ./${SAMPLE_ID}_aligned.bam && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <${ALIGNED_BAM}>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                while : ; do
                    TEST=$(gsutil -m cp ${ALIGNED_BAI} ./${SAMPLE_ID}_aligned.bam.bai && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <${ALIGNED_BAI}>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
            else
                while : ; do
                    TEST=$(gsutil -m cp ${PAV_VCF_GZ} ./${SAMPLE_ID}_pav.vcf.gz && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <${PAV_VCF_GZ}>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                while : ; do
                    TEST=$(gsutil -m cp ${PAV_TBI} ./${SAMPLE_ID}_pav.vcf.gz.tbi && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <${PAV_TBI}>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                while : ; do
                    TEST=$(gsutil -m cp ${PBSV_VCF_GZ} ./${SAMPLE_ID}_pbsv.vcf.gz && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <${PBSV_VCF_GZ}>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                while : ; do
                    TEST=$(gsutil -m cp ${PBSV_TBI} ./${SAMPLE_ID}_pbsv.vcf.gz.tbi && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <${PBSV_TBI}>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                while : ; do
                    TEST=$(gsutil -m cp ${SNIFFLES_VCF_GZ} ./${SAMPLE_ID}_sniffles.vcf.gz && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <${SNIFFLES_VCF_GZ}>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                while : ; do
                    TEST=$(gsutil -m cp ${SNIFFLES_TBI} ./${SAMPLE_ID}_sniffles.vcf.gz.tbi && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <${SNIFFLES_TBI}>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                if [ ~{use_trgt} -eq 1 ]; then
                    while : ; do
                        TEST=$(gsutil -m cp ${TRGT_VCF_GZ} ./${SAMPLE_ID}_trgt.vcf.gz && echo 0 || echo 1)
                        if [ ${TEST} -eq 1 ]; then
                            echo "Error downloading file <${TRGT_VCF_GZ}>. Trying again..."
                            sleep ${GSUTIL_DELAY_S}
                        else
                            break
                        fi
                    done
                    while : ; do
                        TEST=$(gsutil -m cp ${TRGT_TBI} ./${SAMPLE_ID}_trgt.vcf.gz.tbi && echo 0 || echo 1)
                        if [ ${TEST} -eq 1 ]; then
                            echo "Error downloading file <${TRGT_TBI}>. Trying again..."
                            sleep ${GSUTIL_DELAY_S}
                        else
                            break
                        fi
                    done
                fi
            fi
        }
        
        
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
        
        
        # Makes sure that SVLEN and SVTYPE are consistently annotated.
        #
        # Remark: for non-symbolic replacement records (where REF and ALT are
        # both >1bp), `truvari anno svinfo` sets SVLEN to ||REF|-|ALT|| and
        # SVTYPE to INS or DEL based only on SVLEN: this is incorrect, it
        # affects in particular INVs, TRGT VCFs, and VCFs with entire
        # haplotypes, and it leads to wrong SVLEN filtering downstream. The
        # procedure sets SVLEN=max(|REF|,|ALT|)-1 and SVTYPE=UNK in such cases.
        #
        function addSvlenSvtype() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_PREFIX=$2
            
            truvari anno svinfo --minsize 1 ${INPUT_VCF_GZ} | bgzip > ${OUTPUT_PREFIX}_tmp.vcf.gz
            tabix -f ${OUTPUT_PREFIX}_tmp.vcf.gz
            bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${OUTPUT_PREFIX}_tmp.vcf.gz | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                if (substr($5,1,1)!="<" && length($4)>1 && length($5)>1) {
                    printf("%s",$1); \
                    for (i=2; i<=5; i++) printf("\t%s",$i); \
                    max=length($5)>length($4)?length($5):length($4); \
                    printf("\t%d\tUNK\n",max-1); \
                } \
            }' | bgzip -c > ${OUTPUT_PREFIX}_annotations.tsv.gz
            tabix -f -s1 -b2 -e2 ${OUTPUT_PREFIX}_annotations.tsv.gz
            ${TIME_COMMAND} bcftools annotate --annotations ${OUTPUT_PREFIX}_annotations.tsv.gz --columns CHROM,POS,~ID,REF,ALT,INFO/SVLEN,INFO/SVTYPE --output-type z ${OUTPUT_PREFIX}_tmp.vcf.gz > ${OUTPUT_PREFIX}.vcf.gz
            tabix -f ${OUTPUT_PREFIX}.vcf.gz
            rm -f ${OUTPUT_PREFIX}_tmp.vcf.gz* ${OUTPUT_PREFIX}_annotations.tsv.gz*
        }
        
        
        # Puts in canonical form a raw VCF from an SV caller. The procedure
        # creates sorted output files `SAMPLEID_CALLERID_X.vcf.gz`, where X is:
        #
        # sv: non-BND records with length in [MIN_SV_LENGTH..MAX_SV_LENGTH], in
        #     canonical form.
        #
        function CanonizeVcf() {
            local INPUT_VCF_GZ=$1
            local INPUT_TBI=$2
            local SAMPLE_ID=$3
            local CALLER_ID=$4
            local MIN_SV_LENGTH=$5
            local MAX_SV_LENGTH=$6
            local STANDARD_CHROMOSOMES_BED=$7
            local NOT_GAPS_BED=$8
            
            # QUAL is used by truvari collapse to select a representation. We
            # assign values based on which representations we observed to be
            # more accurate in a few examples.
            if [ ${CALLER_ID} = 'trgt' ]; then
                QUAL="5"
            elif [ ${CALLER_ID} = 'pav' ]; then
                QUAL="4"
            elif [ ${CALLER_ID} = 'pbsv' ]; then
                QUAL="3"
            elif [ ${CALLER_ID} = 'sniffles' ]; then
                QUAL="2"
            fi
            
            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            mv ${INPUT_TBI} ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi
            
            # Subsetting to standard chromosomes, or to a given region if any.
            if [ ~{region} != "all" ]; then
                ${TIME_COMMAND} bcftools view --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ~{region} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            else
                ${TIME_COMMAND} bcftools filter --regions-file ${STANDARD_CHROMOSOMES_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            fi
            
            # Ensuring that SVLEN has the correct type for bcftools norm
            bcftools view --header-only ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz | sed 's/ID=SVLEN,Number=.,/ID=SVLEN,Number=A,/g' > ${SAMPLE_ID}_${CALLER_ID}_header.txt
            ${TIME_COMMAND} bcftools reheader --header ${SAMPLE_ID}_${CALLER_ID}_header.txt --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Removing SNVs, if any.
            if [ ${CALLER_ID} = 'pav' ]; then
                ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="SNV"' --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            fi
            
            # Removing non-ALT records from TRGT
            if [ ${CALLER_ID} = 'trgt' ]; then
                ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            fi
            
            # Removing records in reference gaps
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Making sure SVLEN and SVTYPE are consistently annotated
            addSvlenSvtype ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_out
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Isolating BNDs
            ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="BND"' --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Isolating ultra-long records and discarding short records
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1. Main VCF ------------------------------------------------------
            
            # 1.1 Sorting
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.2 Fixing symbolic records
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx5G FixSymbolicRecords ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ~{reference_fa} | bgzip > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.3 Fixing REF
            ${TIME_COMMAND} bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.4 Cleaning REF, ALT, QUAL, FILTER. 
            # - REF and ALT must be uppercase for XGBoost scoring downstream to
            #   work.
            # - QUAL is used by truvari collapse to select a representation. 
            #   Symbolic records are NOT given low quality (it was 1 in Phase 1)
            #   since e.g. all DEL records made by Sniffles are symbolic.
            # - We force every record to PASS, to rule out any filter-dependent
            #   effect in downstream tools.
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanRefAltQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${QUAL} | bgzip > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.5 Removing END, since its values may be inconsistent and make
            # GATK crash downstream.
            ${TIME_COMMAND} bcftools annotate --remove INFO/END --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.6 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --remove-duplicates --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_sv.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_sv.vcf.gz.tbi
        }
        
        
        # Collapses with truvari all files `SAMPLEID_CALLERID_sv.vcf.gz`,
        # creating an output file `SAMPLEID_sv.vcf.gz`.
        #
        function IntrasampleMerge_sv() {
            local SAMPLE_ID=$1
            
            # Remark: the order of the callers in `bcftools merge` affects the
            # value of the SAMPLE column emitted by `truvari collapse --intra`.
            if [ ~{use_trgt} -eq 1 ]; then
                TRGT_STRING="${SAMPLE_ID}_trgt_sv.vcf.gz"
bcftools view ${SAMPLE_ID}_trgt_sv.vcf.gz chr6:575456
            else
                TRGT_STRING=" "
            fi
            
            
bcftools view ${SAMPLE_ID}_pav_sv.vcf.gz chr6:575456
bcftools view ${SAMPLE_ID}_pbsv_sv.vcf.gz chr6:575456
bcftools view ${SAMPLE_ID}_sniffles_sv.vcf.gz chr6:575456

            
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type z ${TRGT_STRING} ${SAMPLE_ID}_pav_sv.vcf.gz ${SAMPLE_ID}_pbsv_sv.vcf.gz ${SAMPLE_ID}_sniffles_sv.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_*_sv.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_in.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0.90 --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_in.vcf > ${SAMPLE_ID}_in.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Ensuring that every record has a unique ID, to enable joining by
            # CHROM,POS,ID in downstream calls to `bcftools annotate`. Using
            # CHROM,POS,REF,ALT can make `bcftools annotate` segfault, and the
            # speed of joining by CHROM,POS,ID is independent of SVLEN.
            (bcftools view --header-only ${SAMPLE_ID}_in.vcf.gz ; bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { printf("%s\t%s\t%d-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,++i,$3,$4,$5,$6,$7,$8,$9,$10); }') | bgzip --compress-level 1 > ${SAMPLE_ID}_out.vcf.gz
            mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz            
            bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 5 || echo "0"
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_sv.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_sv.vcf.gz.tbi
        }
        
        
        function Kanpig() {
            local SAMPLE_ID=$1
            local SEX=$2
            local INPUT_VCF_GZ=$3
            local ALIGNMENTS_BAM=$4

            if [ ${SEX} == "M" ]; then
                PLOIDY_BED=$(echo ~{ploidy_bed_male})
            else
                PLOIDY_BED=$(echo ~{ploidy_bed_female})
            fi
            
            # Remark: kanpig needs --sizemin >= --kmer
            ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_singlesample} --sizemin 10 --sizemax ${INFINITY} --reference ~{reference_fa} --input ${INPUT_VCF_GZ} --reads ${ALIGNMENTS_BAM} --out ${SAMPLE_ID}_out.vcf
            rm -f ${INPUT_VCF_GZ}* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # Sorting
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_in.vcf > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            
            # Discarding records that are not marked as present by kanpig
            N_RECORDS_BEFORE=$( bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz.tbi )
            ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            N_RECORDS_AFTER=$( bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz.tbi )
            PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER} / ${N_RECORDS_BEFORE}" | bc )
            echo "Number of records that are marked as present by kanpig: ${N_RECORDS_AFTER}/${N_RECORDS_BEFORE} (${PERCENT}%)"
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_kanpig.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_kanpig.vcf.gz.tbi
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        export RUST_BACKTRACE="full"
        INFINITY="1000000000"
        
        GetReferenceGaps ~{reference_agp} not_gaps.bed
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        LINE=$( head -n ~{selected_row} chunk.csv | tail -n 1 )
        SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
        SEX=$(echo ${LINE} | cut -d , -f 2)
        
        # Merging
        LocalizeSample ${SAMPLE_ID} 1 ${LINE}
        if [ ~{use_trgt} -eq 1 ]; then
            CanonizeVcf ${SAMPLE_ID}_trgt.vcf.gz ${SAMPLE_ID}_trgt.vcf.gz.tbi ${SAMPLE_ID} trgt ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
        fi
        CanonizeVcf ${SAMPLE_ID}_pav.vcf.gz ${SAMPLE_ID}_pav.vcf.gz.tbi ${SAMPLE_ID} pav ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
        CanonizeVcf ${SAMPLE_ID}_pbsv.vcf.gz ${SAMPLE_ID}_pbsv.vcf.gz.tbi ${SAMPLE_ID} pbsv ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
        CanonizeVcf ${SAMPLE_ID}_sniffles.vcf.gz ${SAMPLE_ID}_sniffles.vcf.gz.tbi ${SAMPLE_ID} sniffles ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
        IntrasampleMerge_sv ${SAMPLE_ID}
        
        # Genotyping
        LocalizeSample ${SAMPLE_ID} 2 ${LINE}
        Kanpig ${SAMPLE_ID} ${SEX} ${SAMPLE_ID}_sv.vcf.gz ${SAMPLE_ID}_aligned.bam
    >>>
    
    output {
        File kanpig_vcf_gz = glob("*_kanpig.vcf.gz")[0]
        File kanpig_tbi = glob("*_kanpig.vcf.gz.tbi")[0]
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}



# Performance with 4 cores and 32GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# truvari bench             
# vcfdist                   200%        11G     6m
#
task PrecisionRecallAnalysis {
    input {
        String sample_id
        File kanpig_vcf_gz
        File kanpig_tbi
        String experiment_id
        File sample_dipcall_vcf_gz
        File sample_dipcall_bed
        
        String chromosome
        
        String remote_outdir
        
        Int bench_method
        Int min_sv_length
        Int max_sv_length
        
        File reference_fa
        File reference_fai
        File reference_agp
        File standard_chromosomes_bed
        File tandem_bed
        File not_tandem_bed
        
        Int n_cpu = 4
        Int ram_size_gb = 16
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
        
        # Subsetting the dipcall VCF to the given chromosome, if specified.
        if [ ~{chromosome} != "all" ]; then
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth.vcf.gz > ~{sample_id}_truth_prime.vcf.gz
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth_tr.vcf.gz > ~{sample_id}_truth_tr_prime.vcf.gz
            ${TIME_COMMAND} bcftools filter --regions ~{chromosome} --regions-overlap pos --output-type z ~{sample_id}_truth_not_tr.vcf.gz > ~{sample_id}_truth_not_tr_prime.vcf.gz
            
            rm -f ~{sample_id}_truth.vcf.gz* ; mv ~{sample_id}_truth_prime.vcf.gz ~{sample_id}_truth.vcf.gz ; bcftools index ~{sample_id}_truth.vcf.gz
            rm -f ~{sample_id}_truth_tr.vcf.gz* ; mv ~{sample_id}_truth_tr_prime.vcf.gz ~{sample_id}_truth_tr.vcf.gz ; bcftools index ~{sample_id}_truth_tr.vcf.gz
            rm -f ~{sample_id}_truth_not_tr.vcf.gz* ; mv ~{sample_id}_truth_not_tr_prime.vcf.gz ~{sample_id}_truth_not_tr.vcf.gz ; bcftools index ~{sample_id}_truth_not_tr.vcf.gz
        fi
        
        # Keeping only records that are genotyped as present
        ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type z ~{kanpig_vcf_gz} > out.vcf.gz
        rm -f ~{kanpig_vcf_gz} ; mv out.vcf.gz in.vcf.gz ; tabix -f in.vcf.gz
        
        # Benchmarking
        Benchmark ~{sample_id} in.vcf.gz ~{experiment_id}_~{min_sv_length}bp
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp '*_'~{experiment_id}_~{min_sv_length}bp'_*.txt' ~{remote_outdir}/precision_recall/ && echo 0 || echo 1)
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
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
