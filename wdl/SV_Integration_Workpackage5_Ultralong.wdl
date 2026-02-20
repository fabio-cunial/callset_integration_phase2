version 1.0


# Performs a simple bcftools merge and truvari collapse of all the ultralong or
# BND intra-sample VCFs.
#
workflow SV_Integration_Workpackage5_Ultralong {
    input {
        File sample_ids
        String suffix = "ultralong"
        Array[String] bi_samples_to_prefer_over_ha
        String truvari_matching_parameters = "--refdist 500 --pctseq 0 --pctsize 0.95 --pctovl 0.0"
        
        String remote_indir_bi
        String remote_indir_ha
        String remote_indir_bcm
        String remote_indir_uw
        String remote_indir_controls_15x
        String remote_indir_controls_30x
        
        Int n_expected_samples_bi
        Int n_expected_samples_ha
        Int n_expected_samples_bcm
        Int n_expected_samples_uw
        Int n_expected_samples_controls_15x
        Int n_expected_samples_controls_30x
        
        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        sample_ids: "Speficies the order of the samples to use in bcftools merge."
        remote_indir_bi: "Without final slash"
        remote_outdir: "Without final slash"
        suffix: "Denoting the type of intra-sample VCFs we want to merge: 'ultralong' or 'bnd'."
        truvari_matching_parameters: "Truvari's definition of a match. Identical to `SV_Integration_Workpackage7.wdl`, except for sequence similarity which is turned off for speed reasons."
    }
    
    call Impl {
        input:
            sample_ids = sample_ids,
            suffix = suffix,
            bi_samples_to_prefer_over_ha = bi_samples_to_prefer_over_ha,
            truvari_matching_parameters = truvari_matching_parameters,
            
            remote_indir_bi = remote_indir_bi,
            remote_indir_ha = remote_indir_ha,
            remote_indir_bcm = remote_indir_bcm,
            remote_indir_uw = remote_indir_uw,
            remote_indir_controls_15x = remote_indir_controls_15x,
            remote_indir_controls_30x = remote_indir_controls_30x,
            
            n_expected_samples_bi = n_expected_samples_bi,
            n_expected_samples_ha = n_expected_samples_ha,
            n_expected_samples_bcm = n_expected_samples_bcm,
            n_expected_samples_uw = n_expected_samples_uw,
            n_expected_samples_controls_15x = n_expected_samples_controls_15x,
            n_expected_samples_controls_30x = n_expected_samples_controls_30x,
            
            remote_outdir = remote_outdir,
            
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, HDD, ultralong VCFs:
#
# TOOL                          CPU     RAM     TIME
# gcloud storage cp             
# bcftools merge level 1        
# bcftools norm level 1         
# bcftools merge level 2        
# bcftools norm level 2         
#
# Peak disk usage (all input files): 
#
task Impl {
    input {
        File sample_ids
        String suffix
        Array[String] bi_samples_to_prefer_over_ha
        String truvari_matching_parameters
        
        String remote_indir_bi
        String remote_indir_ha
        String remote_indir_bcm
        String remote_indir_uw
        String remote_indir_controls_15x
        String remote_indir_controls_30x
        
        Int n_expected_samples_bi
        Int n_expected_samples_ha
        Int n_expected_samples_bcm
        Int n_expected_samples_uw
        Int n_expected_samples_controls_15x
        Int n_expected_samples_controls_30x
        
        Int n_files_per_merge = 100
        String remote_outdir
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 16
        Int disk_size_gb = 50
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
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        #
        function LocalizeFiles() {
            touch all_remote_files.txt
            
            # Ensuring that every input dataset has the expected number of
            # samples in the chunk.
            date 1>&2
            if [ ~{n_expected_samples_bi} -gt 0 ]; then
                gcloud storage ls -l ~{remote_indir_bi}/'*_'~{suffix}'.bcf' | tr -s ' ' | sed 's/^[ ]*//' > bi_files.txt
                N_FILES=$(wc -l < bi_files.txt)
                N_FILES=$(( ${N_FILES} - 1 ))
                if [ ${N_FILES} -ne ~{n_expected_samples_bi} ]; then
                    echo "ERROR: BI has ${N_FILES} files != ~{n_expected_samples_bi}"
                    exit 1
                fi
                head -n ${N_FILES} bi_files.txt >> all_remote_files.txt
            fi
        
            if [ ~{n_expected_samples_ha} -gt 0 ]; then
                gcloud storage ls -l ~{remote_indir_ha}/'*_'~{suffix}'.bcf' | tr -s ' ' | sed 's/^[ ]*//' > ha_files.txt
                N_FILES=$(wc -l < ha_files.txt)
                N_FILES=$(( ${N_FILES} - 1 ))
                if [ ${N_FILES} -ne ~{n_expected_samples_ha} ]; then
                    echo "ERROR: HA has ${N_FILES} files != ~{n_expected_samples_ha}"
                    exit 1
                fi
                head -n ${N_FILES} ha_files.txt >> all_remote_files.txt
            fi
        
            if [ ~{n_expected_samples_bcm} -gt 0 ]; then
                gcloud storage ls -l ~{remote_indir_bcm}/'*_'~{suffix}'.bcf' | tr -s ' ' | sed 's/^[ ]*//' > bcm_files.txt
                N_FILES=$(wc -l < bcm_files.txt)
                N_FILES=$(( ${N_FILES} - 1 ))
                if [ ${N_FILES} -ne ~{n_expected_samples_bcm} ]; then
                    echo "ERROR: BCM has ${N_FILES} files != ~{n_expected_samples_bcm}"
                    exit 1
                fi
                head -n ${N_FILES} bcm_files.txt >> all_remote_files.txt
            fi
        
            if [ ~{n_expected_samples_uw} -gt 0 ]; then
                gcloud storage ls -l ~{remote_indir_uw}/'*_'~{suffix}'.bcf' | tr -s ' ' | sed 's/^[ ]*//' > uw_files.txt
                N_FILES=$(wc -l < uw_files.txt)
                N_FILES=$(( ${N_FILES} - 1 ))
                if [ ${N_FILES} -ne ~{n_expected_samples_uw} ]; then
                    echo "ERROR: UW has ${N_FILES} files != ~{n_expected_samples_uw}"
                    exit 1
                fi
                head -n ${N_FILES} uw_files.txt >> all_remote_files.txt
            fi
        
            if [ ~{n_expected_samples_controls_15x} -gt 0 ]; then
                gcloud storage ls -l ~{remote_indir_controls_15x}/'*_'~{suffix}'.bcf' | tr -s ' ' | sed 's/^[ ]*//' > control_15x_files.txt
                N_FILES=$(wc -l < control_15x_files.txt)
                N_FILES=$(( ${N_FILES} - 1 ))
                if [ ${N_FILES} -lt ~{n_expected_samples_controls_15x} ]; then
                    echo "ERROR: CONTROLS_15X has ${N_FILES} files < ~{n_expected_samples_controls_15x}"
                    exit 1
                elif [ ${N_FILES} -gt ~{n_expected_samples_controls_15x} ]; then
                    echo "WARNING: CONTROLS_15X has ${N_FILES} files > ~{n_expected_samples_controls_15x}"
                fi
                head -n ${N_FILES} control_15x_files.txt >> all_remote_files.txt
            fi
        
            if [ ~{n_expected_samples_controls_30x} -gt 0 ]; then
                gcloud storage ls -l ~{remote_indir_controls_30x}/'*_'~{suffix}'.bcf' | tr -s ' ' | sed 's/^[ ]*//' > control_30x_files.txt
                N_FILES=$(wc -l < control_30x_files.txt)
                N_FILES=$(( ${N_FILES} - 1 ))
                if [ ${N_FILES} -lt ~{n_expected_samples_controls_30x} ]; then
                    echo "ERROR: CONTROLS_30X has ${N_FILES} files < ~{n_expected_samples_controls_30x}"
                    exit 1
                elif [ ${N_FILES} -gt ~{n_expected_samples_controls_30x} ]; then
                    echo "WARNING: CONTROLS_30X has ${N_FILES} files > ~{n_expected_samples_controls_30x}"
                fi
                head -n ${N_FILES} control_30x_files.txt >> all_remote_files.txt
            fi
            date 1>&2
        
            # Failing immediately if the files are too large WRT the available
            # disk. Otherwise the VM may get stuck forever, and this gets worse
            # with preemption.
            AVAILABLE_GB=$(df -h | grep "cromwell_root" | tr -s ' ' | cut -d ' ' -f 4)
            AVAILABLE_GB=${AVAILABLE_GB%G}
            AVAILABLE_GB=${AVAILABLE_GB%.*}
            REMOTE_GB=$(java -cp ~{docker_dir} SumFileSizes all_remote_files.txt)
            SLACK_GB="5"
            REMOTE_GB=$(( ${REMOTE_GB} + ${SLACK_GB} ))
            if [ ${REMOTE_GB} -gt ${AVAILABLE_GB} ]; then
                echo "ERROR: the remote files are larger than the available disk space. Remote files + slack: ${REMOTE_GB}GB. Available disk: ${AVAILABLE_GB}GB."
                exit 1
            fi
            rm -f *_files.txt
        
            # - Localizing all the single-sample VCFs.
            # - Handling samples that occur in multiple input datasets.
            date 1>&2
            mkdir ./input_files/
            if [ ~{n_expected_samples_bi} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_bi}/'*_'~{suffix}'.bcf*' ./input_files/
            fi
            if [ ~{n_expected_samples_ha} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_ha}/'*_'~{suffix}'.bcf*' ./input_files/
            fi
            if [ ~{n_expected_samples_bi} -gt 0 -a ~{n_expected_samples_ha} -gt 0 ]; then
                echo ~{sep="," bi_samples_to_prefer_over_ha} | tr ',' '\n' > bi_samples_to_prefer_over_ha.txt
                rm -f list.txt
                while read SAMPLE_ID; do
                    echo "~{remote_indir_bi}/${SAMPLE_ID}_"~{suffix}".bcf" >> list.txt
                    echo "~{remote_indir_bi}/${SAMPLE_ID}_"~{suffix}".bcf.csi" >> list.txt
                done < bi_samples_to_prefer_over_ha.txt
                cat list.txt | gcloud storage cp -I ./input_files/
            fi
            if [ ~{n_expected_samples_uw} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_uw}/'*_'~{suffix}'.bcf*' ./input_files/
            fi
            if [ ~{n_expected_samples_bcm} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_bcm}/'*_'~{suffix}'.bcf*' ./input_files/
            fi
            if [ ~{n_expected_samples_controls_15x} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_controls_15x}/'*_'~{suffix}'.bcf*' ./input_files/
            fi
            if [ ~{n_expected_samples_controls_30x} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_controls_30x}/'*_'~{suffix}'.bcf*' ./input_files/
            fi
            date 1>&2
            N_DOWNLOADED_SAMPLES=$(ls ./input_files/*.bcf | wc -l)
            N_SAMPLES=$(cat ~{sample_ids} | wc -l)
            if [ ${N_DOWNLOADED_SAMPLES} -lt ${N_SAMPLES} ]; then
                echo "ERROR: The number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is smaller than the number of samples specified (${N_SAMPLES})."
                exit 1
            fi
            df -h
        }
        
        
        # Trivial "hierarchical" bcftools merge with just two steps.
        #
        function BcftoolsMerge() {
            # Step 1
            rm -f list.txt
            while read SAMPLE_ID; do
                echo ./input_files/${SAMPLE_ID}_~{suffix}.bcf >> list.txt
            done < ~{sample_ids}
            split -l ~{n_files_per_merge} -d -a 4 list.txt list_
            N_LIST_FILES=$(ls list_* | wc -l)
            for LIST_FILE in $(ls list_* | sort -V); do
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list ${LIST_FILE} --output-type b --output ${LIST_FILE}_merged.bcf
                ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_merged.bcf
                df -h 1>&2
                xargs --arg-file=${LIST_FILE} --max-lines=1 --max-procs=${N_THREADS} rm -f
                ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ${LIST_FILE}_merged.bcf --output ${LIST_FILE}_normed.bcf
                ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_normed.bcf
                df -h 1>&2
                rm -f ${LIST_FILE}_merged.bcf* ; mv ${LIST_FILE}_normed.bcf ${LIST_FILE}_merged.bcf ; mv ${LIST_FILE}_normed.bcf.csi ${LIST_FILE}_merged.bcf.csi
            done
            
            # Step 2
            ls list_*.bcf | sort -V > list.txt
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type b --output merged.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f merged.bcf
            df -h 1>&2
            xargs --arg-file=list.txt --max-lines=1 --max-procs=${N_THREADS} rm -f
            ls -laht 1>&2
            
            # Making sure no multiallelic record is passed downstream
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b merged.bcf --output normed.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f normed.bcf
            rm -f merged.bcf* ; mv normed.bcf merged.bcf ; mv normed.bcf.csi merged.bcf.csi
        }

        
        # Sets QUAL to the number of samples where a record was discovered, to
        # simulate `--keep common` (which is slow on 10k samples) with `--keep
        # maxqual` in `truvari collapse`. See e.g.:
        #
        # https://github.com/ACEnglish/truvari/issues/220#issuecomment-
        # 2830920205
        #
        # Remark: we use the number of samples rather than AC, since we don't
        # trust genotypes at this stage, and since a record being discovered
        # independently in more samples is more informative than it being
        # genotyped more times in fewer samples.
        #
        function CopyNSamplesToQual() {
            local INPUT_BCF=$1
            local INPUT_CSI=$2

            mv ${INPUT_BCF} in.bcf
            mv ${INPUT_CSI} in.bcf.csi

            # Remark: we cannot join annotations just by ID at this stage,
            # since the IDs in the output of bcftools merge are not necessarily
            # all distinct.
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%COUNT(GT="alt")\n' in.bcf | bgzip -c > annotations.tsv.gz
            tabix -@ ${N_THREADS} -s1 -b2 -e2 annotations.tsv.gz
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations annotations.tsv.gz --columns CHROM,POS,~ID,REF,ALT,QUAL --output-type b in.bcf --output out.bcf
            rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
            rm -f annotations.tsv.gz
    
            mv in.bcf annotated.bcf
            mv in.bcf.csi annotated.bcf.csi
        }
        
        
        # Remark: in theory we should set `--gt all` to avoid collapsing records
        # that are present in the same sample, since we assume that intra-
        # sample merging has already done that upstream. In practice `--gt all`
        # is too slow on 10k samples.
        #
        function TruvariCollapse() {
            local INPUT_BCF=$1
            local INPUT_CSI=$2
            
            # Removing SVLEN from symbolic ALTs, in order not to interfere with
            # `truvari collapse`.
            if [ ~{suffix} = "ultralong" ]; then
                date 1>&2
                ( bcftools view --header-only ${INPUT_BCF} ; bcftools view --no-header ${INPUT_BCF} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    if (substr($0,1,1)!="#" && substr($5,1,1)=="<") $5 = substr($5,1,4) ">"; \
                    printf("%s",$1); \
                    for (i=2; i<=NF; i++) printf("\t%s",$i); \
                    printf("\n"); \
                }' ) | bgzip --compress-level 1 > in.vcf.gz
                date 1>&2
            else
                bcftools view --threads ${N_THREADS} --output-type z ${INPUT_BCF} --output in.vcf.gz
            fi
            bcftools index --threads ${N_THREADS} -f -t in.vcf.gz
            rm -f ${INPUT_BCF} ${INPUT_CSI}

            # Remark: we do not store `removed.vcf` since it's not needed and
            # it can be much bigger than the collapsed output.
            ${TIME_COMMAND} truvari collapse --sizemin 0 --sizemax ${INFINITY} --keep maxqual --gt off ~{truvari_matching_parameters} --input in.vcf.gz --output out.vcf --removed-output /dev/null
            df -h 1>&2
            ls -laht 1>&2
            rm -f in.vcf.gz* ; mv out.vcf in.vcf
        
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type b in.vcf --output out.bcf
            df -h 1>&2
            ls -laht 1>&2
            rm -f in.vcf ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
            
            # Dropping the IDs written by truvari collapse, since they can be
            # very long on a large cohort and needlessly inflate output size.
            date 1>&2
            ( bcftools view --header-only in.bcf ; bcftools view --no-header in.bcf | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { $3=sprintf("%d",i++); print $0 }' ) | bcftools view --output-type b --output out.bcf
            date 1>&2
            rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
                
            mv in.bcf truvari_collapsed.bcf
            mv in.bcf.csi truvari_collapsed.bcf.csi
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        LocalizeFiles
        BcftoolsMerge
        CopyNSamplesToQual merged.bcf merged.bcf.csi
        TruvariCollapse annotated.bcf annotated.bcf.csi
        
        gcloud storage mv truvari_collapsed.'bcf*' ~{remote_outdir}/
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
