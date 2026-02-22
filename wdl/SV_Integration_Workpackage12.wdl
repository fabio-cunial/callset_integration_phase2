version 1.0


# Performs a simple bcftools merge of all the ultralong or BND intra-sample
# VCFs. Stores in output per-chromosome VCFs that should then be split for
# parallel truvari collapse.
#
workflow SV_Integration_Workpackage12 {
    input {
        File sample_ids
        String suffix = "ultralong"
        Array[String] bi_samples_to_prefer_over_ha
        String chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        
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
    }
    
    call Impl {
        input:
            sample_ids = sample_ids,
            suffix = suffix,
            bi_samples_to_prefer_over_ha = bi_samples_to_prefer_over_ha,
            chromosomes = chromosomes,
            
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
# TOOL                           CPU     RAM     TIME
# gcloud storage cp                                3m            // Whole genome
# bcftools merge level 1        300%    600M      10s            // Whole genome
# bcftools norm level 1         300%    300M      10s            // Whole genome
# bcftools merge level 2        200%    3.5G       4m            // Per chr
# bcftools norm level 2         250%      4G       2m            // Per chr
#
# Peak disk usage (all input files): 10G
#
task Impl {
    input {
        File sample_ids
        String suffix
        Array[String] bi_samples_to_prefer_over_ha
        String chromosomes
        
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
        Int ram_size_gb = 8
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
            df -h 1>&2
        }
        
        
        cat << 'END' > chunk_by_chr.sh
#!/bin/bash
INPUT_BCF=$1
CHROMOSOME=$2
mkdir -p ./${CHROMOSOME}/
bcftools view --output-type b ${INPUT_BCF} ${CHROMOSOME} --output ./${CHROMOSOME}/${INPUT_BCF}
bcftools index -f ./${CHROMOSOME}/${INPUT_BCF}
END
        chmod +x chunk_by_chr.sh
        cat << 'END' > chunk_by_region.sh
#!/bin/bash
INPUT_BCF=$1
CHROMOSOME=$2
REGION=$3
CHUNK_ID=$4
bcftools view --regions ${REGION} --regions-overlap pos --output-type b ${INPUT_BCF} --output ./${CHROMOSOME}/chunk_${CHUNK_ID}.bcf
bcftools index -f ./${CHROMOSOME}/chunk_${CHUNK_ID}.bcf
END
        chmod +x chunk_by_region.sh        
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        echo ~{chromosomes} | tr ',' '\n' > chromosomes.txt
        LocalizeFiles
        
        # Trivial "hierarchical" bcftools merge with just two steps.
        # Step 1: merging a few samples at a time over the whole genome.
        rm -f list.txt
        while read SAMPLE_ID; do
            echo ./input_files/${SAMPLE_ID}_~{suffix}.bcf >> list.txt
        done < ~{sample_ids}
        split -l ~{n_files_per_merge} -d -a 4 list.txt list_
        N_LIST_FILES=$(ls list_* | wc -l)
        for LIST_FILE in $(ls list_* | sort -V); do
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list ${LIST_FILE} --output-type b --output ${LIST_FILE}_merged.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_merged.bcf
            xargs --arg-file=${LIST_FILE} --max-lines=1 --max-procs=${N_THREADS} rm -f
            rm -f ${LIST_FILE}
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ${LIST_FILE}_merged.bcf --output ${LIST_FILE}_normed.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_normed.bcf
            rm -f ${LIST_FILE}_merged.bcf*
            ${TIME_COMMAND} xargs --arg-file=chromosomes.txt --max-lines=1 --max-procs=${N_THREADS} ./chunk_by_chr.sh ./${LIST_FILE}_normed.bcf
            rm -f ${LIST_FILE}_normed.bcf*
        done
        rm -rf ./input_files/
        
        # Step 2: merging all samples over each chromosome.
        rm -f files_list.txt
        while read CHROMOSOME; do
            ls ./${CHROMOSOME}/*.bcf | sort -V > list.txt
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type b --output ./${CHROMOSOME}/merged.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ./${CHROMOSOME}/merged.bcf
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ./${CHROMOSOME}/merged.bcf --output ./${CHROMOSOME}/normed.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ./${CHROMOSOME}/normed.bcf
            mv ./${CHROMOSOME}/normed.bcf ${CHROMOSOME}.bcf
            mv ./${CHROMOSOME}/normed.bcf.csi ${CHROMOSOME}.bcf.csi
            echo "${CHROMOSOME}.bcf" >> files_list.txt
            echo "${CHROMOSOME}.bcf.csi" >> files_list.txt
            rm -rf ./${CHROMOSOME}/
        done < chromosomes.txt
        df -h 1>&2
        ls -laht 1>&2
        
        # Uploading
        date 1>&2
        cat files_list.txt | gcloud storage mv -I ~{remote_outdir}/
        date 1>&2
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
