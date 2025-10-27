version 1.0


# Runs `PAV2SVs.wdl`, `Resolve.wdl`, `TruvariIntrasample.wdl`, `Kanpig.wdl` in
# the same VM for multiple samples.
#
workflow Workpackage1 {
    input {
        File sv_integration_chunk_tsv
        Int min_sv_length
        String remote_outdir
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        
        Int n_cpu = 6
        Int ram_size_gb = 8
        Int disk_size_gb = 100
        String kanpig_params_singlesample = "--sizemin 20 --sizemax 10000 --neighdist 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"
    }
    parameter_meta {
        sv_integration_chunk_tsv: "A subset of the rows of table `sv_integration_hg38`, without the header."
        remote_outdir: "Where the output of intra-sample truvari and kanpig is stored for each sample."
    }
    
    call Workpackage1Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            min_sv_length = min_sv_length,
            remote_outdir = remote_outdir,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb,
            kanpig_params_singlesample = kanpig_params_singlesample
    }
    
    output {
    }
}

 
#
task Workpackage1Impl {
    input {
        File sv_integration_chunk_tsv
        
        Int min_sv_length
        String kanpig_params_singlesample
        String remote_outdir
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        
        Int n_cpu
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
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
        # $3 A line of `sv_integration_chunk_tsv`.
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
            fi
        }
        
        
        # Deletes all and only the files downloaded by `LocalizeSample()`.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}_aligned.bam* ${SAMPLE_ID}_pav.vcf.gz* ${SAMPLE_ID}_pbsv.vcf.gz* ${SAMPLE_ID}_sniffles.vcf.gz*
        }
        
        
        # Puts in canonical form a raw VCF from an SV caller. The procedure
        # creates sorted output files `${SAMPLE_ID}_${CALLER_ID}_X.vcf.gz`,
        # where X is:
        #
        # sv: non-BND records with length in [MIN_SV_LENGTH..MAX_SV_LENGTH], in
        #     canonical form;
        # sv_ultralong: non-BND records with length `>MAX_SV_LENGTH`, devoid of
        #               sequence where possible to save space;
        # bnd: BND records, in their original form.
        #
        function CanonizeVcf() {
            local INPUT_VCF_GZ=$1
            local INPUT_TBI=$2
            local SAMPLE_ID=$3
            local CALLER_ID=$4
            local MIN_SV_LENGTH=$5
            local MAX_SV_LENGTH=$6
            
            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            mv ${INPUT_TBI} ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Removing SNVs, if any.
            if [ ${CALLER_ID} = 'pav' ]; then
                ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="SNV"' --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            fi
            
            # Making sure SVLEN and SVTYPE are consistently annotated
            truvari anno svinfo --minsize 1 ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Isolating BNDs
            ${TIME_COMMAND} bcftools filter --include 'SVTYPE="BND"' --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz
            tabix -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz
            ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="BND"' --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Isolating ultra-long calls and discarding short calls
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>'${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz
            tabix -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1. Main VCF ------------------------------------------------------
            
            # 1.1 Sorting
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.2 Fixing symbolic records
            java -cp ~{docker_dir} -Xmx5G FixSymbolicRecords ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ~{reference_fa} ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.3 Fixing REF
            ${TIME_COMMAND} bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.4 Cleaning REF, ALT, QUAL. 
            # - REF and ALT must be uppercase for XGBoost scoring downstream to
            #   work.
            # - QUAL is used by truvari collapse to select a representation. We
            #   assign values based on which representations we observed to be
            #   more accurate in a few examples. Symbolic records are NOT given
            #   low quality (it was 1 in Phase 1) since e.g. all DEL calls made
            #   by Sniffles are symbolic.
            if [ ${CALLER_ID} = 'pav' ]; then
                QUAL="4"
            elif [ ${CALLER_ID} = 'pbsv' ]; then
                QUAL="3"
            elif [ ${CALLER_ID} = 'sniffles' ]; then
                QUAL="2"
            fi
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanRefAltQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${QUAL} ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.5 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --remove-duplicates --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_sv.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_sv.vcf.gz.tbi
            
            # 2. BND VCF -------------------------------------------------------
            
            # 2.1 Sorting 
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz
            
            # 2.2 Fixing REF
            ${TIME_COMMAND} bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz
            
            # 2.3 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --remove-duplicates --output-type z ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz
            
            
            # 3. Ultralong VCF -------------------------------------------------
            
            # 3.1 Sorting
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz
            
            # 3.2 Removing sequence (lossless).
            ${TIME_COMMAND} java -cp ~{docker_dir} RemoveRefAlt ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz
            
            # 3.3 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --remove-duplicates --output-type z ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz ; tabix -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz
        }
        
        
        
        
        
        
        Symolic variants in ultralong VCF may be destroyed by bcftools merge!!!! Use our script.
        
        Remember that now SUPP order is PAV, PBSV, SNIFFLES!!!!!!!!!!
        

        
        # The procedure creates merged output files ${SAMPLE_ID}_X.vcf.gz`,
        # where X is:
        #
        # sv: non-BND records with length in [MIN_SV_LENGTH..MAX_SV_LENGTH],
        #     collapsed with truvari;
        # sv_ultralong: non-BND records with length `>MAX_SV_LENGTH`; collapsed
        #               with truvari but without sequence similarity.
        # bnd: BND records, not collapsed with truvari since `truvari collapse`
        #      does not yet work on BNDs.
        #
        function IntrasampleMerge() {
            local SAMPLE_ID=$1
            
            
            # 1. Main VCF ------------------------------------------------------
            
            # Remark: the order of the callers in `bcftools merge` affects the
            # value of the SAMPLE column emitted by `truvari collapse --intra`.
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples -output-type z ${SAMPLE_ID}_pav_sv.vcf.gz ${SAMPLE_ID}_pbsv_sv.vcf.gz ${SAMPLE_ID}_sniffles_sv.vcf.gz > ${SAMPLE_ID}_tmp1.vcf.gz
            tabix -f ${SAMPLE_ID}_tmp1.vcf.gz
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_tmp1.vcf.gz > ${SAMPLE_ID}_bcftools_merged.vcf.gz
            tabix -f ${SAMPLE_ID}_bcftools_merged.vcf.gz
            rm -f ${SAMPLE_ID}_tmp1.vcf.gz*
            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_bcftools_merged.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0.90 --pctsize 0.90 --sizemin 0 --sizemax -1 --output ${SAMPLE_ID}_tmp2.vcf
            rm -f ${SAMPLE_ID}_bcftools_merged.vcf.gz*
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_tmp2.vcf > ${SAMPLE_ID}_sv.vcf.gz
            tabix -f ${SAMPLE_ID}_sv.vcf.gz
            rm -f ${SAMPLE_ID}_tmp2.vcf
            
            # 2. BND VCF -------------------------------------------------------
            
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples -output-type z ${SAMPLE_ID}_pav_bnd.vcf.gz ${SAMPLE_ID}_pbsv_bnd.vcf.gz ${SAMPLE_ID}_sniffles_bnd.vcf.gz > ${SAMPLE_ID}_tmp1.vcf.gz
            tabix -f ${SAMPLE_ID}_tmp1.vcf.gz
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_tmp1.vcf.gz > ${SAMPLE_ID}_bnd.vcf.gz
            tabix -f ${SAMPLE_ID}_bnd.vcf.gz
            rm -f ${SAMPLE_ID}_tmp1.vcf.gz*
            
            # 3. Ultralong VCF -------------------------------------------------
            
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples -output-type z ${SAMPLE_ID}_pav_ultralong.vcf.gz ${SAMPLE_ID}_pbsv_ultralong.vcf.gz ${SAMPLE_ID}_sniffles_ultralong.vcf.gz > ${SAMPLE_ID}_tmp1.vcf.gz
            tabix -f ${SAMPLE_ID}_tmp1.vcf.gz
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_tmp1.vcf.gz > ${SAMPLE_ID}_bcftools_merged.vcf.gz
            tabix -f ${SAMPLE_ID}_bcftools_merged.vcf.gz
            rm -f ${SAMPLE_ID}_tmp1.vcf.gz*
            
            # Removing SVLEN from symbolic ALTs for `truvari collapse`.
            bcftools view --header-only ${SAMPLE_ID}_bcftools_merged.vcf.gz > ${SAMPLE_ID}_tmp2.vcf
            ${TIME_COMMAND} bcftools view --no-header ${SAMPLE_ID}_bcftools_merged.vcf.gz | awk '{ \
                $5=substr($5,1,4) ">"; \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' >> ${SAMPLE_ID}_tmp2.vcf
            rm -f ${SAMPLE_ID}_bcftools_merged.vcf.gz*
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} ${SAMPLE_ID}_tmp2.vcf
            tabix -f ${SAMPLE_ID}_tmp2.vcf.gz
            
            # Collapsing
            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_tmp2.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0 --pctsize 0.90 --sizemin 0 --sizemax -1 --output ${SAMPLE_ID}_tmp3.vcf
            rm -f ${SAMPLE_ID}_tmp2.vcf.gz*
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_tmp3.vcf > ${SAMPLE_ID}_ultralong.vcf.gz
            tabix -f ${SAMPLE_ID}_ultralong.vcf.gz
            rm -f ${SAMPLE_ID}_tmp3.vcf
        }
        
        
        # Unpacks truvari's SUPP field into 3 INFO tags.
        #
        function TransferSupp() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local OUTPUT_VCF_GZ=$3
            
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SUPP]\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                printf("%s",$1); \
                for (i=2; i<=NF-1; i++) printf("\t%s",$i); \
                if ($6=="0") printf("\t0\t0\t0");
                else if ($6=="1") printf("\t0\t0\t1");
                else if ($6=="2") printf("\t0\t1\t0");
                else if ($6=="3") printf("\t0\t1\t1");
                else if ($6=="4") printf("\t1\t0\t0");
                else if ($6=="5") printf("\t1\t0\t1");
                else if ($6=="6") printf("\t1\t1\t0");
                else if ($6=="7") printf("\t1\t1\t1");
                printf("\n"); \
            }' | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
            tabix -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by PAV">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by Sniffles2">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> ${SAMPLE_ID}_header.txt
            ${TIME_COMMAND} bcftools annotate --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns CHROM,POS,ID,REF,ALT,INFO/SUPP_PAV,INFO/SUPP_SNIFFLES,INFO/SUPP_PBSV ${INPUT_VCF_GZ} --output-type z > ${OUTPUT_VCF_GZ}
            tabix -f ${OUTPUT_VCF_GZ}
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt
        }
        
        
        function Kanpig() {
            local SAMPLE_ID=$1
            local SEX=$2
            local INPUT_VCF_GZ=$3
            local ALIGNMENTS_BAM=$4
    
            # Formatting the merged VCF
            HAS_SUPP=$(bcftools view --header-only ${INPUT_VCF_GZ} | grep '##FORMAT=<ID=SUPP,' && echo 1 || echo 0)
            if [ ${HAS_SUPP} -eq 0 ]; then
                mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_formatted.vcf.gz
                mv ${INPUT_VCF_GZ}.tbi ${SAMPLE_ID}_formatted.vcf.gz.tbi
            else
                TransferSupp ${SAMPLE_ID} ${INPUT_VCF_GZ} ${SAMPLE_ID}_formatted.vcf.gz
            fi

            # Making sure there is just one occurrence of '##fileformat=' in the
            # header (otherwise kanpig complains).
            echo '##fileformat=VCFv4.2' > ${SAMPLE_ID}_cleaned.vcf
            bcftools view --header-only ${SAMPLE_ID}_formatted.vcf.gz | grep -v '##fileformat=' >> ${SAMPLE_ID}_cleaned.vcf
            bcftools view --no-header ${SAMPLE_ID}_formatted.vcf.gz >> ${SAMPLE_ID}_cleaned.vcf
            bgzip -@ ${N_THREADS} ${SAMPLE_ID}_cleaned.vcf
            tabix -f ${SAMPLE_ID}_cleaned.vcf.gz
            rm -f ${SAMPLE_ID}_formatted.vcf.gz*

            # Kanpig
            if [ ${SEX} == "M" ]; then
                PLOIDY_BED=$(echo ~{ploidy_bed_male})
            else
                PLOIDY_BED=$(echo ~{ploidy_bed_female})
            fi
            export RUST_BACKTRACE="full"
            ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_singlesample} --reference ~{reference_fa} --input ${SAMPLE_ID}_cleaned.vcf.gz --reads ${ALIGNMENTS_BAM} --out ${SAMPLE_ID}_tmp1.vcf.gz
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_tmp1.vcf.gz > ${SAMPLE_ID}_kanpig.vcf.gz
            tabix -f ${SAMPLE_ID}_kanpig.vcf.gz
            rm -f ${SAMPLE_ID}_tmp1.vcf.gz*
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            SEX=$(echo ${LINE} | cut -d , -f 2)
            
            LocalizeSample ${SAMPLE_ID} 1 ${LINE}
            PAV2SVs ${SAMPLE_ID} ${SAMPLE_ID}_pav.vcf.gz
            rm -f ${SAMPLE_ID}_pav.vcf.gz*
            source activate truvari4
            Resolve ${SAMPLE_ID} pbsv ${SAMPLE_ID}_pbsv.vcf.gz
            rm -f ${SAMPLE_ID}_pbsv.vcf.gz*
            Resolve ${SAMPLE_ID} sniffles ${SAMPLE_ID}_sniffles.vcf.gz
            rm -f ${SAMPLE_ID}_sniffles.vcf.gz*
            Resolve ${SAMPLE_ID} pav ${SAMPLE_ID}_pav_sv.vcf.gz
            rm -f ${SAMPLE_ID}_pav_sv.vcf.gz*
            conda deactivate
            source activate truvari5
            TruvariIntrasample ${SAMPLE_ID} ${SAMPLE_ID}_pbsv_resolved.vcf.gz ${SAMPLE_ID}_sniffles_resolved.vcf.gz ${SAMPLE_ID}_pav_resolved.vcf.gz
            conda deactivate
            
            LocalizeSample ${SAMPLE_ID} 2 ${LINE}
            Kanpig ${SAMPLE_ID} ${SEX} ${SAMPLE_ID}_truvari_collapsed.vcf.gz ${SAMPLE_ID}_aligned.bam
            gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ${SAMPLE_ID}_kanpig.vcf.'gz*' ~{remote_outdir}
            
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done < chunk.csv
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
