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
        
        Int n_cpu = 8
        Int ram_size_gb = 16
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
        String remote_outdir
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        
        Int n_cpu
        Int ram_size_gb
        Int disk_size_gb
        String kanpig_params_singlesample
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int compression_level = 1
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 3)
            ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 4)
            PAV_BED=$(echo ${LINE} | cut -d , -f 5)
            PAV_TBI=$(echo ${LINE} | cut -d , -f 6)
            PAV_VCF_GZ=$(echo ${LINE} | cut -d , -f 7)
            PBSV_TBI=$(echo ${LINE} | cut -d , -f 8)
            PBSV_VCF_GZ=$(echo ${LINE} | cut -d , -f 9)
            SNIFFLES_TBI=$(echo ${LINE} | cut -d , -f 10)
            SNIFFLES_VCF_GZ=$(echo ${LINE} | cut -d , -f 11)
            
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
        }
        
        
        # Deletes all and only the files downloaded by `LocalizeSample()`.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}_aligned.bam* ${SAMPLE_ID}_pav.vcf.gz* ${SAMPLE_ID}_pbsv.vcf.gz* ${SAMPLE_ID}_sniffles.vcf.gz*
        }
        
        
        function PAV2SVs() {
            local SAMPLE_ID=$1
            local PAV_VCF_GZ=$2
            
            # Removing multiallelic records, if any.
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type v ${PAV_VCF_GZ} > ${SAMPLE_ID}_tmp1.vcf
        
            # Keeping only long calls
            ${TIME_COMMAND} java -cp ~{docker_dir} PAV2SVs ${SAMPLE_ID}_tmp1.vcf ~{min_sv_length} ${SAMPLE_ID}_pav_sv.vcf ${SAMPLE_ID}_pav_snp.vcf
            rm -f ${SAMPLE_ID}_tmp1.vcf ${SAMPLE_ID}_pav_snp.vcf
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} ${SAMPLE_ID}_pav_sv.vcf
            tabix -f ${SAMPLE_ID}_pav_sv.vcf.gz
        }
        
        
        function Resolve() {
            local SAMPLE_ID=$1
            local CALLER_ID=$2
            local INPUT_VCF_GZ=$3
            
            # Removing multiallelic records from the input
            bcftools norm --multiallelics - --output-type z ${INPUT_VCF_GZ} > ${SAMPLE_ID}_tmp1.vcf.gz
            tabix -f ${SAMPLE_ID}_tmp1.vcf.gz
        
            # Step 1 - clean up the VCF
            # - Assigns quality scores to each SV caller's result
            #  - pav 4
            #  - pbsv 3
            #  - sniffles 2
            # - Resolves any symbolic alts (e.g. `<DEL>` with the sequence from
            #   the reference)
            #   - symbolic variants are given quality score of 1
            # - Filters out `<CNV>` from pbsv and `<INS>` from sniffles
            # - Filters out variants greater than 100Mbp
            # - Fills in blank genotypes with `0`
            # - Filters out BND variants
            # - Turn some deletion/insertion pairs into inversions
            # The quality scores are set based on which variant representations
            # we believe are generally more accurate with higher being better.
            python ~{docker_dir}/resolve.py ${SAMPLE_ID}_tmp1.vcf.gz ~{reference_fa} | bcftools norm --check-ref s --fasta-ref ~{reference_fa} -N -m-any | bcftools view -i "SVTYPE != 'BND'" -O z -o ${SAMPLE_ID}_${CALLER_ID}_resolved.vcf.gz
            tabix -f ${SAMPLE_ID}_${CALLER_ID}_resolved.vcf.gz
            rm -f ${SAMPLE_ID}_tmp1.vcf.gz*
        
            # $inversion_guesser.py$ creates a DEL with wrong ALT when an INS
            # is followed by a DEL that cancels it out. I.e. this input:
            #
            # chr2    14066975        pbsv.INS.669    T       TATATATATGATATATATATCATATATATGATATATATGATATATATATCATATATATG
            # chr2    14066975        pbsv.DEL.670    TATATATATGATATATATATCATATATATGATATATATGATATATATATCATATATATG
            #
            # gives the following output:
            #
            # chr2    14066975        pbsv.DEL.670    TATATATATGATATATATATCATATATATGATATATATGATATATATATCATATATATG     TATATATATGATATATATATCATATATATGATATATATGATATATATATCATATATATG
            #python ~{docker_dir}/inversion_guesser.py -i $prename -o $outname
        }
        
        
        function TruvariIntrasample() {
            local SAMPLE_ID=$1
            local PBSV_VCF_GZ=$2
            local SNIFFLES_VCF_GZ=$3
            local PAV_VCF_GZ=$4
            
            # Step 1 - Merging
            # Pastes the samples together in the order of the preferred
            # genotypes. That is to say, this creates a three sample VCF with
            # sample columns from pbsv, sniffles, pav_sv
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples -O z -o ${SAMPLE_ID}_tmp1.vcf.gz ${PBSV_VCF_GZ} ${SNIFFLES_VCF_GZ} ${PAV_VCF_GZ}
            tabix -f ${SAMPLE_ID}_tmp1.vcf.gz
        
            # Step 2 - Removing multiallelic records. We observed that they are
            # created in Step 1 sometimes, e.g.:
            #
            # 2024-03-19 20:52:28,548 [ERROR] Cannot compare multi-allelic
            # records. Please split
            # 2024-03-19 20:52:28,548 [ERROR] line
            # chr4	137168756	pbsv.INS.2751;chr4-137168757-DEL-52	A   ACGTATGTGTATACGTATACATATACGCGTATATACATACGTATACATATACG,A	4	PASS	SVTYPE=INS;SVLEN=52;SVANN=TANDEM;ID=chr4-137168757-DEL-52;TIG_REGION=h2tg007223l:91206-91206;QUERY_STRAND=+;HOM_REF=0,23;HOM_TIG=0,23;INVScore=0.981132;AC=2,1	GT:AD:DP:SAC	1/1:1,8,.:9:1,0,3,5	./.:.:.:.	0|2:.:.:.
            #
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_tmp1.vcf.gz > ${SAMPLE_ID}_bcftools_merged.vcf.gz
            tabix -f ${SAMPLE_ID}_bcftools_merged.vcf.gz
            rm -f ${SAMPLE_ID}_tmp1.vcf.gz*

            # Step 3 - Collapsing
            ${TIME_COMMAND} truvari collapse -i ${SAMPLE_ID}_bcftools_merged.vcf.gz --sizemin 0 --sizemax 1000000 --keep maxqual --gt het --intra --pctseq 0.90 --pctsize 0.90 --refdist 500 --output ${SAMPLE_ID}_tmp2.vcf
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_tmp2.vcf > ${SAMPLE_ID}_truvari_collapsed.vcf.gz
            tabix -f ${SAMPLE_ID}_truvari_collapsed.vcf.gz
            rm -f ${SAMPLE_ID}_tmp2.vcf
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
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z ${SAMPLE_ID}_tmp1.vcf.gz > ${SAMPLE_ID}_kanpig.vcf.gz
            tabix -f ${SAMPLE_ID}_kanpig.vcf.gz
            rm -f ${SAMPLE_ID}_tmp1.vcf.gz*
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            SEX=$(echo ${LINE} | cut -d , -f 2)
            
            LocalizeSample ${SAMPLE_ID} ${LINE}
            
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
            Kanpig ${SAMPLE_ID} ${SEX} ${SAMPLE_ID}_truvari_collapsed.vcf.gz ${ALIGNED_BAM}
            
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
