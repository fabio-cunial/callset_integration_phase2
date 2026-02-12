version 1.0


# Performs a bcftools merge of all the VCFs of a given chromosome chunk.
#
workflow SV_Integration_Workpackage5 {
    input {
        Int chunk_id
        File sample_ids
        Array[String] bi_samples_to_prefer_over_ha
        
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
        
        Int merge_mode
        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        sample_ids: "Speficies the order of the samples to use in bcftools merge."
        remote_indir_bi: "Without final slash"
        remote_outdir: "Without final slash"
        merge_mode: "1: standard bcftools merge (CHROM,POS,REF,ALT). 2: bcftools merge by ID only."
    }
    
    call Impl {
        input:
            chunk_id = chunk_id,
            sample_ids = sample_ids,
            bi_samples_to_prefer_over_ha = bi_samples_to_prefer_over_ha,
            
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
            
            merge_mode = merge_mode,
            remote_outdir = remote_outdir,
            
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, first 30 MB chunk of chr6:
#
# TOOL                          CPU     RAM     TIME
# gcloud storage ls                             30 s
# gcloud storage cp             280%    100 M    1 m
#
# Merge by CHROM,POS,REF,ALT:
# bcftools merge level 1        100%    300 M    3 s          // 100 files
# bcftools norm level 1         300%     50 M    1 s
# bcftools merge level 2        170%      1 G   20 m          // 127 files
# bcftools norm level 2         300%     11 G    6 m
#
# Peak disk usage (all input files of chunk 0): 2 GB
#
# Merge by ID:
# bcftools merge level 1        400%    1.5G     5 s          // 100 files
# bcftools merge level 2        300%    2.5G    30 s          // 126 files
#
# Peak disk usage (all input files of chunk 0): 74G
#
task Impl {
    input {
        Int chunk_id
        File sample_ids
        Array[String] bi_samples_to_prefer_over_ha
        
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
        
        Int merge_mode
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
        function LocalizeChunkFiles() {
            touch all_remote_files.txt
            
            # Ensuring that every input dataset has the expected number of
            # samples in the chunk.
            date 1>&2
            if [ ~{n_expected_samples_bi} -gt 0 ]; then
                gcloud storage ls -l ~{remote_indir_bi}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > bi_files.txt
                N_FILES=$(wc -l < bi_files.txt)
                N_FILES=$(( ${N_FILES} - 1 ))
                if [ ${N_FILES} -ne ~{n_expected_samples_bi} ]; then
                    echo "ERROR: BI has ${N_FILES} files != ~{n_expected_samples_bi}"
                    exit 1
                fi
                head -n ${N_FILES} bi_files.txt >> all_remote_files.txt
            fi
        
            if [ ~{n_expected_samples_ha} -gt 0 ]; then
                gcloud storage ls -l ~{remote_indir_ha}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > ha_files.txt
                N_FILES=$(wc -l < ha_files.txt)
                N_FILES=$(( ${N_FILES} - 1 ))
                if [ ${N_FILES} -ne ~{n_expected_samples_ha} ]; then
                    echo "ERROR: HA has ${N_FILES} files != ~{n_expected_samples_ha}"
                    exit 1
                fi
                head -n ${N_FILES} ha_files.txt >> all_remote_files.txt
            fi
        
            if [ ~{n_expected_samples_bcm} -gt 0 ]; then
                gcloud storage ls -l ~{remote_indir_bcm}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > bcm_files.txt
                N_FILES=$(wc -l < bcm_files.txt)
                N_FILES=$(( ${N_FILES} - 1 ))
                if [ ${N_FILES} -ne ~{n_expected_samples_bcm} ]; then
                    echo "ERROR: BCM has ${N_FILES} files != ~{n_expected_samples_bcm}"
                    exit 1
                fi
                head -n ${N_FILES} bcm_files.txt >> all_remote_files.txt
            fi
        
            if [ ~{n_expected_samples_uw} -gt 0 ]; then
                gcloud storage ls -l ~{remote_indir_uw}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > uw_files.txt
                N_FILES=$(wc -l < uw_files.txt)
                N_FILES=$(( ${N_FILES} - 1 ))
                if [ ${N_FILES} -ne ~{n_expected_samples_uw} ]; then
                    echo "ERROR: UW has ${N_FILES} files != ~{n_expected_samples_uw}"
                    exit 1
                fi
                head -n ${N_FILES} uw_files.txt >> all_remote_files.txt
            fi
        
            if [ ~{n_expected_samples_controls_15x} -gt 0 ]; then
                gcloud storage ls -l ~{remote_indir_controls_15x}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > control_15x_files.txt
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
                gcloud storage ls -l ~{remote_indir_controls_30x}/chunk_~{chunk_id}/'*.bcf' | tr -s ' ' | sed 's/^[ ]*//' > control_30x_files.txt
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
        
            # - Localizing all the samples for the given chunk.
            # - Handling samples that occur in multiple input datasets.
            date 1>&2
            mkdir ./input_files/
            if [ ~{n_expected_samples_bi} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_bi}/chunk_~{chunk_id}/'*' ./input_files/
            fi
            if [ ~{n_expected_samples_ha} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_ha}/chunk_~{chunk_id}/'*' ./input_files/
            fi
            if [ ~{n_expected_samples_bi} -gt 0 -a ~{n_expected_samples_ha} -gt 0 ]; then
                echo ~{sep="," bi_samples_to_prefer_over_ha} | tr ',' '\n' > bi_samples_to_prefer_over_ha.txt
                rm -f list.txt
                while read SAMPLE_ID; do
                    echo "~{remote_indir_bi}/chunk_~{chunk_id}/${SAMPLE_ID}.bcf" >> list.txt
                    echo "~{remote_indir_bi}/chunk_~{chunk_id}/${SAMPLE_ID}.bcf.csi" >> list.txt
                done < bi_samples_to_prefer_over_ha.txt                
                cat list.txt | gcloud storage cp -I ./input_files/
            fi
            if [ ~{n_expected_samples_uw} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_uw}/chunk_~{chunk_id}/'*' ./input_files/
            fi
            if [ ~{n_expected_samples_bcm} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_bcm}/chunk_~{chunk_id}/'*' ./input_files/
            fi
            if [ ~{n_expected_samples_controls_15x} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_controls_15x}/chunk_~{chunk_id}/'*' ./input_files/
            fi
            if [ ~{n_expected_samples_controls_30x} -gt 0 ]; then
                ${TIME_COMMAND} gcloud storage cp ~{remote_indir_controls_30x}/chunk_~{chunk_id}/'*' ./input_files/
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
        
        
        # Trivial "hierarchical" merge with just two steps.
        #
        function MergeChunkFiles() {
            if [ ~{merge_mode} -eq 1 ]; then
                MERGE_FLAG="none"
            elif [ ~{merge_mode} -eq 2 ]; then
                MERGE_FLAG="id"
            else
                echo "ERROR: Merge mode unknown."
                exit 1
            fi
            
            # Step 1
            rm -f list.txt
            while read SAMPLE_ID; do
                echo ./input_files/${SAMPLE_ID}.bcf >> list.txt
            done < ~{sample_ids}
            split -l ~{n_files_per_merge} -d -a 4 list.txt list_
            N_LIST_FILES=$(ls list_* | wc -l)
            for LIST_FILE in $(ls list_* | sort -V); do
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge ${MERGE_FLAG} --file-list ${LIST_FILE} --output-type b --output ${LIST_FILE}_merged.bcf
                ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_merged.bcf
                df -h 1>&2
                xargs --arg-file=${LIST_FILE} --max-lines=1 --max-procs=${N_THREADS} rm -f
                if [ ~{merge_mode} -eq 1 ]; then
                    ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ${LIST_FILE}_merged.bcf --output ${LIST_FILE}_normed.bcf
                    ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_normed.bcf
                    df -h 1>&2
                    rm -f ${LIST_FILE}_merged.bcf* ; mv ${LIST_FILE}_normed.bcf ${LIST_FILE}_merged.bcf ; mv ${LIST_FILE}_normed.bcf.csi ${LIST_FILE}_merged.bcf.csi
                fi
            done
            
            # Step 2
            ls list_*.bcf | sort -V > list.txt
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge ${MERGE_FLAG} --file-list list.txt --output-type b --output ~{chunk_id}_merged.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ~{chunk_id}_merged.bcf
            df -h 1>&2
            xargs --arg-file=list.txt --max-lines=1 --max-procs=${N_THREADS} rm -f
            ls -laht 1>&2
            
            # Making sure no multiallelic record is passed downstream. This is
            # not needed when merging by ID, by construction.
            if [ ~{merge_mode} -eq 1 ]; then
                ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ~{chunk_id}_merged.bcf --output ~{chunk_id}_normed.bcf
                ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ~{chunk_id}_normed.bcf
                rm -f ~{chunk_id}_merged.bcf* ; mv ~{chunk_id}_normed.bcf ~{chunk_id}_merged.bcf ; mv ~{chunk_id}_normed.bcf.csi ~{chunk_id}_merged.bcf.csi
            fi
            
            # Removing records that are REF in all samples. This is not needed
            # in the standard merge, since at that step of the pipeline every
            # input record is ALT in some sample by construction.
            if [ ~{merge_mode} -eq 2 ]; then
                ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type b ~{chunk_id}_merged.bcf --output ~{chunk_id}_cleaned.bcf
                ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ~{chunk_id}_cleaned.bcf
                N_RECORDS=$(bcftools index --nrecords ~{chunk_id}_merged.bcf)
                N_ALT_RECORDS=$(bcftools index --nrecords ~{chunk_id}_cleaned.bcf)
                PERCENT=$( echo "scale=2; 100 * ${N_ALT_RECORDS} / ${N_RECORDS}" | bc )
                echo "${N_ALT_RECORDS},${N_RECORDS},${PERCENT},Number of cohort-VCF records that are marked as ALT in >=1 sample by kanpig" > ~{chunk_id}_n_alt.csv
                rm -f ~{chunk_id}_merged.bcf* ; mv ~{chunk_id}_cleaned.bcf ~{chunk_id}_merged.bcf ; mv ~{chunk_id}_cleaned.bcf.csi ~{chunk_id}_merged.bcf.csi
            fi
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        LocalizeChunkFiles
        MergeChunkFiles
        gcloud storage mv ~{chunk_id}_merged.bcf ~{remote_outdir}/chunk_~{chunk_id}.bcf
        gcloud storage mv ~{chunk_id}_merged.bcf.csi ~{remote_outdir}/chunk_~{chunk_id}.bcf.csi
        gcloud storage mv ~{chunk_id}_n_alt.csv ~{remote_outdir}/
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
