version 1.0


# Merges by ID the re-genotyped personalized BCFs of all samples, limited to a
# given chromosome chunk. Handles samples that occur in multiple centers.
#
workflow SV_Integration_Workpackage10 {
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
        
        String remote_outdir
    }
    parameter_meta {
        sample_ids: "Speficies the order of the samples used by bcftools merge."
        bi_samples_to_prefer_over_ha: "If a sample occurs both in BI and HA, the HA one is usually preferred, except for these samples."
        remote_indir_bi: "Without final slash"
        remote_outdir: "Without final slash"
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
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, one 10 MB chunk of chr6:
#
# TOOL                          CPU     RAM     TIME
# bcftools merge level 1        400%    1.5G    5s          // 100 files
# bcftools merge level 2        300%    2.5G    30s         // 126 files
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
        
        Int n_files_per_merge = 100
        String remote_outdir
        
        Int n_cpu = 16
        Int ram_size_gb = 8
        Int disk_size_gb = 200
    }
    parameter_meta {
        n_files_per_merge: "Number of BCFs to be merged in the first level of a two-step bcftools merge."
        n_cpu: "The main part that takes advantage of multiple cores is file download (which takes ~1h with 16 logical cores)."
        disk_size_gb: "The largest chunk on chr6 takes 500GB"
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"

        
        # Ensuring that every input dataset has the expected number of samples
        # in the chunk.
        gsutil ls -l ~{remote_indir_bi}/'*_chunk_'~{chunk_id}.bcf | tr -s ' ' | sed 's/^[ ]*//' > bi_files.txt
        N_FILES=$(wc -l < bi_files.txt)
        N_FILES=$(( ${N_FILES} - 1 ))
        if [ ${N_FILES} -ne ~{n_expected_samples_bi} ]; then
            echo "ERROR: BI has ${N_FILES} files != ~{n_expected_samples_bi}"
            exit 1
        fi
        head -n ${N_FILES} bi_files.txt >> all_remote_files.txt
        
        gsutil ls -l ~{remote_indir_ha}/'*_chunk_'~{chunk_id}.bcf | tr -s ' ' | sed 's/^[ ]*//' > ha_files.txt
        N_FILES=$(wc -l < ha_files.txt)
        N_FILES=$(( ${N_FILES} - 1 ))
        if [ ${N_FILES} -ne ~{n_expected_samples_ha} ]; then
            echo "ERROR: HA has ${N_FILES} files != ~{n_expected_samples_ha}"
            exit 1
        fi
        head -n ${N_FILES} ha_files.txt >> all_remote_files.txt
        
        gsutil ls -l ~{remote_indir_bcm}/'*_chunk_'~{chunk_id}.bcf | tr -s ' ' | sed 's/^[ ]*//' > bcm_files.txt
        N_FILES=$(wc -l < bcm_files.txt)
        N_FILES=$(( ${N_FILES} - 1 ))
        if [ ${N_FILES} -ne ~{n_expected_samples_bcm} ]; then
            echo "ERROR: BCM has ${N_FILES} files != ~{n_expected_samples_bcm}"
            exit 1
        fi
        head -n ${N_FILES} bcm_files.txt >> all_remote_files.txt
        
        gsutil ls -l ~{remote_indir_uw}/'*_chunk_'~{chunk_id}.bcf | tr -s ' ' | sed 's/^[ ]*//' > uw_files.txt
        N_FILES=$(wc -l < uw_files.txt)
        N_FILES=$(( ${N_FILES} - 1 ))
        if [ ${N_FILES} -ne ~{n_expected_samples_uw} ]; then
            echo "ERROR: UW has ${N_FILES} files != ~{n_expected_samples_uw}"
            exit 1
        fi
        head -n ${N_FILES} uw_files.txt >> all_remote_files.txt
        
        gsutil ls -l ~{remote_indir_controls_15x}/'*_chunk_'~{chunk_id}.bcf | tr -s ' ' | sed 's/^[ ]*//' > control_15x_files.txt
        N_FILES=$(wc -l < control_15x_files.txt)
        N_FILES=$(( ${N_FILES} - 1 ))
        if [ ${N_FILES} -ne ~{n_expected_samples_controls_15x} ]; then
            echo "ERROR: CONTROLS_15X has ${N_FILES} files != ~{n_expected_samples_controls_15x}"
            exit 1
        fi
        head -n ${N_FILES} control_15x_files.txt >> all_remote_files.txt
        
        gsutil ls -l ~{remote_indir_controls_30x}/'*_chunk_'~{chunk_id}.bcf | tr -s ' ' | sed 's/^[ ]*//' > control_30x_files.txt
        N_FILES=$(wc -l < control_30x_files.txt)
        N_FILES=$(( ${N_FILES} - 1 ))
        if [ ${N_FILES} -ne ~{n_expected_samples_controls_30x} ]; then
            echo "ERROR: CONTROLS_30X has ${N_FILES} files != ~{n_expected_samples_controls_30x}"
            exit 1
        fi
        head -n ${N_FILES} control_30x_files.txt >> all_remote_files.txt
        
        # Failing immediately if the files are too large WRT the available
        # disk. Otherwise the VM may get stuck forever, and this is even worse
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
        mkdir ./input_bcfs/
        gsutil -m cp ~{remote_indir_bi}/'*'_chunk_~{chunk_id}.'bcf*' ./input_bcfs/
        gsutil -m cp ~{remote_indir_ha}/'*'_chunk_~{chunk_id}.'bcf*' ./input_bcfs/
        echo ~{sep="," bi_samples_to_prefer_over_ha} | tr ',' '\n' > bi_samples_to_prefer_over_ha.txt
        rm -f list.txt
        while read SAMPLE_ID; do
            echo "~{remote_indir_bi}/${SAMPLE_ID}_chunk_~{chunk_id}.bcf" >> list.txt
            echo "~{remote_indir_bi}/${SAMPLE_ID}_chunk_~{chunk_id}.bcf.csi" >> list.txt
        done < bi_samples_to_prefer_over_ha.txt
        cat list.txt | gsutil -m cp -I ./input_bcfs/
        gsutil -m cp ~{remote_indir_uw}/'*'_chunk_~{chunk_id}.'bcf*' ./input_bcfs/
        gsutil -m cp ~{remote_indir_bcm}/'*'_chunk_~{chunk_id}.'bcf*' ./input_bcfs/
        gsutil -m cp ~{remote_indir_controls_15x}/'*'_chunk_~{chunk_id}.'bcf*' ./input_bcfs/
        gsutil -m cp ~{remote_indir_controls_30x}/'*'_chunk_~{chunk_id}.'bcf*' ./input_bcfs/
        N_DOWNLOADED_SAMPLES=$(ls ./input_bcfs/*_chunk_~{chunk_id}.bcf | wc -l)
        N_SAMPLES=$(cat ~{sample_ids} | wc -l)
        if [ ${N_DOWNLOADED_SAMPLES} -lt ${N_SAMPLES} ]; then
            echo "Error: the number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is smaller than the number of samples specified (${N_SAMPLES})."
            exit 1
        fi
        df -h
        
        # Merging records by ID, since the records in every BCF originate from
        # the same cohort BCF, which had distinct IDs.
        rm -f list.txt
        while read SAMPLE_ID; do
            echo ./input_bcfs/${SAMPLE_ID}_chunk_~{chunk_id}.bcf >> list.txt
        done < ~{sample_ids}
        split -l ~{n_files_per_merge} -d -a 4 list.txt list_
        for LIST_FILE in $(ls list_* | sort -V); do
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge id --output-type b --file-list ${LIST_FILE} > ${LIST_FILE}.bcf
            bcftools index ${LIST_FILE}.bcf
            df -h
            while read INPUT_BCF; do
                rm -f ${INPUT_BCF}*
            done < ${LIST_FILE}
        done
        ls list_*.bcf > list.txt
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge id --output-type b --file-list list.txt > chunk_~{chunk_id}.bcf
        bcftools index chunk_~{chunk_id}.bcf
        df -h
        ls -laht
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp chunk_~{chunk_id}.'bcf*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading merged BCF. Trying again..."
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
        preemptible: 3
    }
}
