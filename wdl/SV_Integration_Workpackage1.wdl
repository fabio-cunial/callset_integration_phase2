version 1.0


# Runs `PAV2SVs.wdl`, `Resolve.wdl`, `TruvariIntrasample.wdl`, `Kanpig.wdl` in
# the same VM for multiple samples.
#
workflow SV_Integration_Workpackage1 {
    input {
        File sv_integration_chunk_tsv
        String region = "all"
        String remote_outdir
        
        Int min_sv_length = 20
        Int max_sv_length = 10000
        String kanpig_params_singlesample = "--neighdist 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"
        
        File training_resource_vcf_gz
        File training_resource_tbi
        File training_resource_bed
        
        File reference_fa
        File reference_fai
        File standard_chromosomes_bed
        File autosomes_bed
        File reference_agp
        File ploidy_bed_female
        File ploidy_bed_male
        
        String docker_image = "fcunial/callset_integration_phase2_workpackages:v3"
    }
    parameter_meta {
        sv_integration_chunk_tsv: "A subset of the rows of table `sv_integration_hg38`, without the header."
        region: "Only consider VCF records in this genomic region. Set to 'all' to disable."
        remote_outdir: "Without final slash. Where the output of intra-sample truvari and kanpig is stored for each sample."
        max_sv_length: "Calls above this length are deemed 'ultralong', are not given to kanpig re-genotyping, and are processed separately."
        training_resource_vcf_gz: "We assume that the training resource VCF has already been subset to the correct length range upstream."
        training_resource_bed: "Training resource calls can belong only to these regions. Typically a high-confidence dipcall BED, or a BED derived from intersecting multiple dipcall BEDs."
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            region = region,
            remote_outdir = remote_outdir,
            
            min_sv_length = min_sv_length,
            max_sv_length = max_sv_length,
            kanpig_params_singlesample = kanpig_params_singlesample,
            
            training_resource_vcf_gz = training_resource_vcf_gz,
            training_resource_tbi = training_resource_tbi,
            training_resource_bed = training_resource_bed,
            
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            standard_chromosomes_bed = standard_chromosomes_bed,
            autosomes_bed = autosomes_bed,
            reference_agp = reference_agp,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            
            docker_image = docker_image
    }
    
    output {
    }
}


# Memory bottlenecks (measured on a 32GB VM):
#
# FixSymbolicRecords           5 GB
# CleanRefAltQual            250 MB
# RemoveRefAlt               200 MB
# truvari collapse           100 MB
# kanpig                     400 MB
#
# Multicore bottlenecks (measured on a 16-CPU VM):
# 
# bcftools merge             300 %
# bgzip                      500 %
# kanpig                     900 %
#
# Kanpig runtimes:
#
# genome, 16 CPUs, 32GB       2 m
# chr6, 6 CPUs, 8GB          20 s
#
task Impl {
    input {
        File sv_integration_chunk_tsv
        String region
        String remote_outdir
        
        Int min_sv_length
        Int max_sv_length
        String kanpig_params_singlesample
        
        File training_resource_vcf_gz
        File training_resource_tbi
        File training_resource_bed
        
        File reference_fa
        File reference_fai
        File standard_chromosomes_bed
        File autosomes_bed
        File reference_agp
        File ploidy_bed_female
        File ploidy_bed_male
        
        String docker_image
        Int n_cpu = 6
        Int ram_size_gb = 8
        Int disk_size_gb = 256
        Int preemptible_number = 4
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
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
        
        
        # Puts in canonical form a raw VCF from an SV caller. The procedure
        # creates sorted output files `SAMPLEID_CALLERID_X.vcf.gz`, where X is:
        #
        # sv: non-BND records with length in [MIN_SV_LENGTH..MAX_SV_LENGTH], in
        #     canonical form;
        # sv_ultralong: non-BND records with length >MAX_SV_LENGTH, devoid of
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
            local STANDARD_CHROMOSOMES_BED=$7
            local NOT_GAPS_BED=$8
            
            # QUAL is used by truvari collapse to select a representation. We
            # assign values based on which representations we observed to be
            # more accurate in a few test examples.
            if [ ${CALLER_ID} = 'pav' ]; then
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
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            else
                ${TIME_COMMAND} bcftools filter --regions-file ${STANDARD_CHROMOSOMES_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            fi
            
            # Removing records in reference gaps
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Ensuring that SVLEN has the correct type for bcftools norm
            bcftools view --header-only ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz | sed 's/ID=SVLEN,Number=.,/ID=SVLEN,Number=A,/g' > ${SAMPLE_ID}_${CALLER_ID}_header.txt
            ${TIME_COMMAND} bcftools reheader --header ${SAMPLE_ID}_${CALLER_ID}_header.txt --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Removing SNVs, if any.
            if [ ${CALLER_ID} = 'pav' ]; then
                ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="SNV"' --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            fi
            
            # Making sure SVLEN and SVTYPE are consistently annotated
            ${TIME_COMMAND} truvari anno svinfo --minsize 1 ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            ${TIME_COMMAND} java -cp ~{docker_dir} FixTruvariAnnoSvlen ${SAMPLE_ID}_${CALLER_ID}_in.vcf | bgzip > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Isolating BNDs
            ${TIME_COMMAND} bcftools filter --include 'SVTYPE="BND"' --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz
            tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz
            ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="BND"' --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Isolating ultra-long records and discarding short records
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>'${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz
            tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1. Main VCF ------------------------------------------------------
            
            # 1.1 Sorting
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.2 Fixing symbolic records
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx5G FixSymbolicRecords ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ~{reference_fa} | bgzip > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.3 Fixing REF
            ${TIME_COMMAND} bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.4 Cleaning REF, ALT, QUAL, FILTER. 
            # - REF and ALT must be uppercase for XGBoost scoring downstream to
            #   work.
            # - QUAL is used by truvari collapse to select a representation. 
            #   Symbolic records are NOT given low quality (it was 1 in Phase 1)
            #   since e.g. all DEL records made by Sniffles are symbolic.
            # - We force every record to PASS, to rule out any filter-dependent
            #   effect in downstream tools.
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanRefAltQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${QUAL} | bgzip > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.5 Removing END, since its values may be inconsistent and make
            # GATK crash downstream.
            ${TIME_COMMAND} bcftools annotate --remove INFO/END --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 1.6 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --remove-duplicates --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_sv.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_sv.vcf.gz.tbi
            
            # 2. BND VCF -------------------------------------------------------
            
            # 2.1 Sorting 
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # Remark: we do not run the following command, since it seems to
            # destroy BNDs ALTs (example: N]chr5:181473415] ->
            # GNcNNNNNNNNNNNNNN ):
            #
            # bcftools norm --check-ref s --fasta-ref ~{reference_fa}
            # --do-not-normalize
            
            # 2.2 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --rm-dup exact --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 2.3 Adding QUAL and forcing every record to PASS
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${QUAL} | bgzip > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz.tbi
            
            # 3. Ultralong VCF -------------------------------------------------
            
            # 3.1 Sorting
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 3.2 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --remove-duplicates --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            # 3.3 Removing sequence (lossless).
            # QUAL is used by truvari collapse to select a representation.
            ${TIME_COMMAND} java -cp ~{docker_dir} RemoveRefAlt ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${QUAL} ~{reference_fai} | bgzip > ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz.tbi
        }
        
        
        # Collapses with truvari all files `SAMPLEID_CALLERID_sv.vcf.gz`,
        # creating an output file `SAMPLEID_sv.vcf.gz`.
        #
        function IntrasampleMerge_sv() {
            local SAMPLE_ID=$1
            
            # Remark: the order of the callers in `bcftools merge` affects the
            # value of the SAMPLE column emitted by `truvari collapse --intra`.
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type z ${SAMPLE_ID}_pav_sv.vcf.gz ${SAMPLE_ID}_pbsv_sv.vcf.gz ${SAMPLE_ID}_sniffles_sv.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_*_sv.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_in.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0.90 --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_in.vcf > ${SAMPLE_ID}_in.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Ensuring that every record has a unique ID, to enable joining by
            # CHROM,POS,ID in downstream calls to `bcftools annotate`. Using
            # CHROM,POS,REF,ALT can make `bcftools annotate` segfault, and the
            # speed of joining by CHROM,POS,ID is independent of SVLEN.
            (bcftools view --header-only ${SAMPLE_ID}_in.vcf.gz ; bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { printf("%s\t%s\t%d-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,++i,$3,$4,$5,$6,$7,$8,$9,$10); }') | bgzip --compress-level 1 > ${SAMPLE_ID}_out.vcf.gz
            mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz            
            bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 5 || echo "0"
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_sv.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_sv.vcf.gz.tbi
        }
        
        
        # Collapses with truvari all files `SAMPLEID_CALLERID_ultralong.vcf.gz`,
        # creating an output file `SAMPLEID_ultralong.vcf.gz`.
        #
        # Remark: `truvari collapse` is run without taking sequence similarity
        # into account, so different INS/DUP/CNV sequences of similar length at
        # similar POS may be wrongly collapsed. We tolerate this for speed
        # reasons.
        #
        function IntrasampleMerge_ultralong() {
            local SAMPLE_ID=$1
            
            # Remark: the order of the callers in `bcftools merge` affects the
            # value of the SAMPLE column emitted by `truvari collapse --intra`.
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type z ${SAMPLE_ID}_pav_ultralong.vcf.gz ${SAMPLE_ID}_pbsv_ultralong.vcf.gz ${SAMPLE_ID}_sniffles_ultralong.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_*_ultralong.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Removing SVLEN from symbolic ALTs, in order not to interfere with
            # `truvari collapse`.
            bcftools view --header-only ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf
            ${TIME_COMMAND} bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                if (substr($0,1,1)!="#" && substr($5,1,1)=="<") $5 = substr($5,1,4) ">"; \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' >> ${SAMPLE_ID}_out.vcf
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_in.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0 --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_in.vcf > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Adding SVLEN back into symbolic ALTs, to avoid overcollapse in the
            # cohort-level bcftools merge downstream.
            ${TIME_COMMAND} java -cp ~{docker_dir} AddSvlenToSymbolicAlt ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_ultralong.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_ultralong.vcf.gz.tbi
        }
        
        
        # Collapses with truvari all files `SAMPLEID_CALLERID_bnd.vcf.gz`,
        # creating an output file `SAMPLEID_bnd.vcf.gz`.
        #
        function IntrasampleMerge_bnd() {
            local SAMPLE_ID=$1
            
            # Remark: the order of the callers in `bcftools merge` affects the
            # value of the SAMPLE column emitted by `truvari collapse --intra`.
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type z ${SAMPLE_ID}_pav_bnd.vcf.gz ${SAMPLE_ID}_pbsv_bnd.vcf.gz ${SAMPLE_ID}_sniffles_bnd.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_*_bnd.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_in.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0.90 --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            
            # Remark: BNDs are not collapsed with truvari, since `truvari
            # collapse` does not work on BNDs.
            #${TIME_COMMAND} java -cp ~{docker_dir} CollapseSamples ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
            #rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            #${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --rm-dup exact --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            #rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_in.vcf > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_bnd.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_bnd.vcf.gz.tbi
        }
        
        
        # Copies truvari's SUPP field from SAMPLE to three tags in INFO. This is
        # necessary, since kanpig overwrites the SAMPLE field.
        #
        function CopySuppToInfo() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local OUTPUT_VCF_GZ=$3
            
            bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SUPP]\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
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
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by pav">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by sniffles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> ${SAMPLE_ID}_header.txt
            # Remark: the order of the callers is now the reverse of the one in
            # which they were bcftools-merged.
            ${TIME_COMMAND} bcftools annotate --annotations ${SAMPLE_ID}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns CHROM,POS,~ID,REF,ALT,INFO/SUPP_SNIFFLES,INFO/SUPP_PBSV,INFO/SUPP_PAV --output-type z ${INPUT_VCF_GZ} > ${OUTPUT_VCF_GZ}
            tabix -@ ${N_THREADS} -f ${OUTPUT_VCF_GZ}
            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt ${INPUT_VCF_GZ}*
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
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Discarding records that are not marked as present by kanpig
            N_RECORDS_BEFORE_KANPIG=$( bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz.tbi )
            ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            N_RECORDS_AFTER_KANPIG=$( bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz.tbi )
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_kanpig.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_kanpig.vcf.gz.tbi
            
            # Printing debug information
            PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER_KANPIG} / ${N_RECORDS_BEFORE_KANPIG}" | bc )
            echo "${N_RECORDS_AFTER},${N_RECORDS_BEFORE},${PERCENT},Number of records that are marked as ALT by kanpig" > ${SAMPLE_ID}_kanpig.csv
            N_HETS_IN_AUTOSOMES=$( bcftools query --format '%ID' --include 'GT="het"' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l )
            N_RECORDS_IN_AUTOSOMES=$( bcftools query --format '%ID' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l )
            PERCENT=$( echo "scale=2; 100 * ${N_HETS_IN_AUTOSOMES} / ${N_RECORDS_IN_AUTOSOMES}" | bc )
            echo "${N_HETS_IN_AUTOSOMES},${N_RECORDS_IN_AUTOSOMES},${PERCENT},Number of records in autosomes that are marked as HET by kanpig" >> ${SAMPLE_ID}_kanpig.csv
            ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigWindows ${SAMPLE_ID}_kanpig.vcf.gz | bgzip > ${SAMPLE_ID}_kanpig.bed.gz
        }
        
        
        # Copies the following kanpig fields from SAMPLE to INFO:
        #
        # KS_1, KS_2, SQ, GQ, DP, AD_NON_ALT, AD_ALL
        #
        # This is necessary, since XGBoost downstream uses only INFO fields.
        #
        function CopyKanpigFieldsToInfo() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            # Creating new header lines
            touch ${SAMPLE_ID}_header.txt
            for FIELD in SQ GQ DP
            do
                bcftools view --header-only ${INPUT_VCF_GZ} | grep ID="${FIELD}," | sed -e 's/FORMAT/INFO/g' >> ${SAMPLE_ID}_header.txt
            done
            echo '##INFO=<ID=AD_NON_ALT,Number=1,Type=Integer,Description="Coverage for non-alternate alleles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=AD_ALL,Number=1,Type=Integer,Description="Coverage for all alleles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=KS_1,Number=1,Type=Integer,Description="Kanpig score 1">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=KS_2,Number=1,Type=Integer,Description="Kanpig score 2">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=GT_COUNT,Number=1,Type=Integer,Description="GT converted to an integer in {0,1,2}.">' >> ${SAMPLE_ID}_header.txt
            
            # Copying fields from FORMAT to INFO. Every record is assumed to
            # have a distinct ID, which is enforced by the steps of the
            # pipeline upstream.
            bcftools query -f '%CHROM\t%POS\t%ID\t[%KS]\t[%SQ]\t[%GQ]\t[%DP]\t[%AD]\t[%GT]\t%INFO/SUPP_PBSV\t%INFO/SUPP_SNIFFLES\t%INFO/SUPP_PAV\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                KS_1=-1; KS_2=-1; \
                p=0; \
                for (i=1; i<=length($4); i++) { \
                    if (substr($4,i,1)==",") { p=i; break; } \
                } \
                if (p==0) { KS_1=$4; KS_2=$4; } \
                else { KS_1=substr($4,1,p-1); KS_2=substr($4,p+1); } \
                if (KS_1==".") KS_1=-1; \
                if (KS_2==".") KS_2=-1; \
                \
                SQ=$5; \
                if (SQ==".") SQ=-1; \
                \
                GQ=$6; \
                if (GQ==".") GQ=-1; \
                \
                DP=$7; \
                if (DP==".") DP=-1; \
                \
                AD_NON_ALT=-1; AD_ALL=1; \
                p=0; \
                for (i=1; i<=length($8); i++) { \
                    if (substr($8,i,1)==",") { p=i; break; } \
                } \
                if (p==0) { AD_NON_ALT=$8; AD_ALL=$8; } \
                else { AD_NON_ALT=substr($8,1,p-1); AD_ALL=substr($8,p+1); } \
                if (AD_NON_ALT==".") AD_NON_ALT=-1; \
                if (AD_ALL==".") AD_ALL=-1; \
                \
                GT_COUNT=-1; \
                if ($9=="0/0" || $9=="0|0" || $9=="./."  || $9==".|." || $9=="./0" || $9==".|0" || $9=="0/." || $9=="0|." || $9=="0" || $9==".") GT_COUNT=0; \
                else if ($9=="0/1" || $9=="0|1" || $9=="1/0" || $9=="1|0" || $9=="./1" || $9==".|1" || $9=="1/." || $9=="1|." || $9=="1") GT_COUNT=1; \
                else if ($9=="1/1" || $9=="1|1") GT_COUNT=2; \
                \
                printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$1,$2,$3,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,$10,$11,$12); \
            }' | bgzip -c > ${SAMPLE_ID}_format.tsv.gz
            tabix -@ ${N_THREADS} -s1 -b2 -e2 ${SAMPLE_ID}_format.tsv.gz
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_format.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns CHROM,POS,~ID,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,SUPP_PBSV,SUPP_SNIFFLES,SUPP_PAV --output-type z ${INPUT_VCF_GZ} > ${SAMPLE_ID}_out.vcf.gz
            mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 5 || echo "0"
            
            rm -f ${INPUT_VCF_GZ}*
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_kanpig.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_kanpig.vcf.gz.tbi
            
            # Removing temporary files
            rm -f ${SAMPLE_ID}_header.txt ${SAMPLE_ID}_format.tsv.gz 
        }
        
        
        # Extracts every record that has a stringent `truvari bench` match with
        # some records in the resource.
        #
        # Remark: we use `--pick single` to force every resource record to be
        # matched with at most one sample record, which is hopefully the
        # most similar to it. This is because we assume that using a
        # contaminated training set in XGBoost downstream is worse than using a
        # slightly smaller training set.
        #
        function GetTrainingRecords() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            ${TIME_COMMAND} truvari bench -b ~{training_resource_vcf_gz} -c ${INPUT_VCF_GZ} --includebed ~{training_resource_bed} --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick single -o ${SAMPLE_ID}_truvari/
            mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_training.vcf.gz
            mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz.tbi ${SAMPLE_ID}_training.vcf.gz.tbi
            
            # Removing temporary files
            rm -rf ./${SAMPLE_ID}_truvari/
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        export RUST_BACKTRACE="full"
        INFINITY="1000000000"
        truvari --help 1>&2
        ~{docker_dir}/kanpig --version 1>&2
        
        GetReferenceGaps ~{reference_agp} not_gaps.bed
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            SEX=$(echo ${LINE} | cut -d , -f 2)
            
            # Skipping the sample if it has already been processed
            TEST1=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}_kanpig.vcf.gz || echo "0" )
            TEST2=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}_ultralong.vcf.gz || echo "0" )
            TEST3=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}_bnd.vcf.gz || echo "0" )
            if [ ${TEST1} != "0" -a ${TEST2} != "0" -a ${TEST3} != "0" ]; then
                continue
            fi
            
            # Merging
            LocalizeSample ${SAMPLE_ID} 1 ${LINE}
            CanonizeVcf ${SAMPLE_ID}_pav.vcf.gz ${SAMPLE_ID}_pav.vcf.gz.tbi ${SAMPLE_ID} pav ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            CanonizeVcf ${SAMPLE_ID}_pbsv.vcf.gz ${SAMPLE_ID}_pbsv.vcf.gz.tbi ${SAMPLE_ID} pbsv ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            CanonizeVcf ${SAMPLE_ID}_sniffles.vcf.gz ${SAMPLE_ID}_sniffles.vcf.gz.tbi ${SAMPLE_ID} sniffles ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            IntrasampleMerge_sv ${SAMPLE_ID}
            IntrasampleMerge_ultralong ${SAMPLE_ID}
            IntrasampleMerge_bnd ${SAMPLE_ID}
            
            # Genotyping
            LocalizeSample ${SAMPLE_ID} 2 ${LINE}
            CopySuppToInfo ${SAMPLE_ID} ${SAMPLE_ID}_sv.vcf.gz ${SAMPLE_ID}_sv_supp.vcf.gz
            Kanpig ${SAMPLE_ID} ${SEX} ${SAMPLE_ID}_sv_supp.vcf.gz ${SAMPLE_ID}_aligned.bam
            CopyKanpigFieldsToInfo ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz
            
            # Marking training records
            GetTrainingRecords ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz
            
            # Converting to BCF and uploading            
            bcftools view --output-type b ${SAMPLE_ID}_ultralong.vcf.gz > ${SAMPLE_ID}_ultralong.bcf
            bcftools index --threads ${N_THREADS} ${SAMPLE_ID}_ultralong.bcf
            rm -f ${SAMPLE_ID}_ultralong.vcf.'gz*'
            bcftools view --output-type b ${SAMPLE_ID}_bnd.vcf.gz > ${SAMPLE_ID}_bnd.bcf
            bcftools index --threads ${N_THREADS} ${SAMPLE_ID}_bnd.bcf
            rm -f ${SAMPLE_ID}_bnd.vcf.'gz*'
            gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ${SAMPLE_ID}_kanpig.vcf.'gz*' ${SAMPLE_ID}_kanpig.bed.gz ${SAMPLE_ID}_kanpig.csv ${SAMPLE_ID}_training.vcf.'gz*' ${SAMPLE_ID}_ultralong.'bcf*' ${SAMPLE_ID}_bnd.'bcf*' ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done < chunk.csv
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
    }
}
