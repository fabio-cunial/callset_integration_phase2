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
        String ultralong_chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        String bnd_chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM,chr1_KI270706v1_random,chr1_KI270707v1_random,chr1_KI270708v1_random,chr1_KI270709v1_random,chr1_KI270710v1_random,chr1_KI270711v1_random,chr1_KI270712v1_random,chr1_KI270713v1_random,chr1_KI270714v1_random,chr2_KI270715v1_random,chr2_KI270716v1_random,chr3_GL000221v1_random,chr4_GL000008v2_random,chr5_GL000208v1_random,chr9_KI270717v1_random,chr9_KI270718v1_random,chr9_KI270719v1_random,chr9_KI270720v1_random,chr11_KI270721v1_random,chr14_GL000009v2_random,chr14_GL000225v1_random,chr14_KI270722v1_random,chr14_GL000194v1_random,chr14_KI270723v1_random,chr14_KI270724v1_random,chr14_KI270725v1_random,chr14_KI270726v1_random,chr15_KI270727v1_random,chr16_KI270728v1_random,chr17_GL000205v2_random,chr17_KI270729v1_random,chr17_KI270730v1_random,chr22_KI270731v1_random,chr22_KI270732v1_random,chr22_KI270733v1_random,chr22_KI270734v1_random,chr22_KI270735v1_random,chr22_KI270736v1_random,chr22_KI270737v1_random,chr22_KI270738v1_random,chr22_KI270739v1_random,chrY_KI270740v1_random,chrUn_KI270302v1,chrUn_KI270304v1,chrUn_KI270303v1,chrUn_KI270305v1,chrUn_KI270322v1,chrUn_KI270320v1,chrUn_KI270310v1,chrUn_KI270316v1,chrUn_KI270315v1,chrUn_KI270312v1,chrUn_KI270311v1,chrUn_KI270317v1,chrUn_KI270412v1,chrUn_KI270411v1,chrUn_KI270414v1,chrUn_KI270419v1,chrUn_KI270418v1,chrUn_KI270420v1,chrUn_KI270424v1,chrUn_KI270417v1,chrUn_KI270422v1,chrUn_KI270423v1,chrUn_KI270425v1,chrUn_KI270429v1,chrUn_KI270442v1,chrUn_KI270466v1,chrUn_KI270465v1,chrUn_KI270467v1,chrUn_KI270435v1,chrUn_KI270438v1,chrUn_KI270468v1,chrUn_KI270510v1,chrUn_KI270509v1,chrUn_KI270518v1,chrUn_KI270508v1,chrUn_KI270516v1,chrUn_KI270512v1,chrUn_KI270519v1,chrUn_KI270522v1,chrUn_KI270511v1,chrUn_KI270515v1,chrUn_KI270507v1,chrUn_KI270517v1,chrUn_KI270529v1,chrUn_KI270528v1,chrUn_KI270530v1,chrUn_KI270539v1,chrUn_KI270538v1,chrUn_KI270544v1,chrUn_KI270548v1,chrUn_KI270583v1,chrUn_KI270587v1,chrUn_KI270580v1,chrUn_KI270581v1,chrUn_KI270579v1,chrUn_KI270589v1,chrUn_KI270590v1,chrUn_KI270584v1,chrUn_KI270582v1,chrUn_KI270588v1,chrUn_KI270593v1,chrUn_KI270591v1,chrUn_KI270330v1,chrUn_KI270329v1,chrUn_KI270334v1,chrUn_KI270333v1,chrUn_KI270335v1,chrUn_KI270338v1,chrUn_KI270340v1,chrUn_KI270336v1,chrUn_KI270337v1,chrUn_KI270363v1,chrUn_KI270364v1,chrUn_KI270362v1,chrUn_KI270366v1,chrUn_KI270378v1,chrUn_KI270379v1,chrUn_KI270389v1,chrUn_KI270390v1,chrUn_KI270387v1,chrUn_KI270395v1,chrUn_KI270396v1,chrUn_KI270388v1,chrUn_KI270394v1,chrUn_KI270386v1,chrUn_KI270391v1,chrUn_KI270383v1,chrUn_KI270393v1,chrUn_KI270384v1,chrUn_KI270392v1,chrUn_KI270381v1,chrUn_KI270385v1,chrUn_KI270382v1,chrUn_KI270376v1,chrUn_KI270374v1,chrUn_KI270372v1,chrUn_KI270373v1,chrUn_KI270375v1,chrUn_KI270371v1,chrUn_KI270448v1,chrUn_KI270521v1,chrUn_GL000195v1,chrUn_GL000219v1,chrUn_GL000220v1,chrUn_GL000224v1,chrUn_KI270741v1,chrUn_GL000226v1,chrUn_GL000213v1,chrUn_KI270743v1,chrUn_KI270744v1,chrUn_KI270745v1,chrUn_KI270746v1,chrUn_KI270747v1,chrUn_KI270748v1,chrUn_KI270749v1,chrUn_KI270750v1,chrUn_KI270751v1,chrUn_KI270752v1,chrUn_KI270753v1,chrUn_KI270754v1,chrUn_KI270755v1,chrUn_KI270756v1,chrUn_KI270757v1,chrUn_GL000214v1,chrUn_KI270742v1,chrUn_GL000216v2,chrUn_GL000218v1,chrEBV"
        
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
        sample_ids: "Specifies the order of the samples to use in bcftools merge."
        remote_indir_bi: "Without final slash"
        remote_outdir: "Without final slash"
        suffix: "Denoting the type of intra-sample VCFs we want to merge: 'ultralong' or 'bnd'."
    }
    
    call Impl {
        input:
            sample_ids = sample_ids,
            suffix = suffix,
            bi_samples_to_prefer_over_ha = bi_samples_to_prefer_over_ha,
            ultralong_chromosomes = ultralong_chromosomes,
            bnd_chromosomes = bnd_chromosomes,
            
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
# BndCanonize 
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
        String ultralong_chromosomes
        String bnd_chromosomes
        
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
        RAM_PER_THREAD_MB=$(( (~{ram_size_gb} * 1024 - 500) / ${N_THREADS} ))
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeFiles() {
            touch all_remote_files.txt

            local N_FILES
            
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
            local AVAILABLE_GB
            AVAILABLE_GB=$(df -h | grep "cromwell_root" | tr -s ' ' | cut -d ' ' -f 4)
            AVAILABLE_GB=${AVAILABLE_GB%G}
            AVAILABLE_GB=${AVAILABLE_GB%.*}
            local REMOTE_GB
            REMOTE_GB=$(java -cp ~{docker_dir} SumFileSizes all_remote_files.txt)
            local SLACK_GB="5"
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
                if [ -s bi_samples_to_prefer_over_ha.txt ]; then
                    rm -f list.txt
                    local SAMPLE_ID
                    while read -u 5 SAMPLE_ID; do
                        echo "~{remote_indir_bi}/${SAMPLE_ID}_"~{suffix}".bcf" >> list.txt
                        echo "~{remote_indir_bi}/${SAMPLE_ID}_"~{suffix}".bcf.csi" >> list.txt
                    done 5< bi_samples_to_prefer_over_ha.txt
                    cat list.txt | gcloud storage cp -I ./input_files/
                fi
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
            local N_DOWNLOADED_SAMPLES=$(find ./input_files -name '*.bcf' | wc -l)
            local N_SAMPLES=$(cat ~{sample_ids} | wc -l)
            if [ ${N_DOWNLOADED_SAMPLES} -lt ${N_SAMPLES} ]; then
                echo "ERROR: The number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is smaller than the number of samples specified (${N_SAMPLES})."
                exit 1
            fi
            df -h 1>&2
        }
        
        
        cat << 'END' > chunk_by_chr.sh
#!/bin/bash
set -euxo pipefail

INPUT_BCF=$1
CHROMOSOME=$2

mkdir -p ./${CHROMOSOME}/
bcftools view --output-type b ${INPUT_BCF} ${CHROMOSOME} --output ./${CHROMOSOME}/${INPUT_BCF}
bcftools index -f ./${CHROMOSOME}/${INPUT_BCF}
END
        chmod +x chunk_by_chr.sh


        cat << 'END' > bnd_canonize.sh
#!/bin/bash
set -euxo pipefail

DOCKER_DIR=$1
RAM_MB=$2
TMP_DIR=$3
INPUT_BCF=$4

FILENAME=$(basename ${INPUT_BCF} .bcf)
bcftools view --output-type z ${INPUT_BCF} --output ${TMP_DIR}/${FILENAME}.vcf.gz
java -cp ${DOCKER_DIR} -Xmx${RAM_MB}M BndCanonize ${TMP_DIR}/${FILENAME}.vcf.gz 1> ${TMP_DIR}/${FILENAME}_canonized.vcf 2> ${TMP_DIR}/${FILENAME}_canonize.log
bcftools sort --max-mem ${RAM_MB}M --output-type b ${TMP_DIR}/${FILENAME}_canonized.vcf --output ${TMP_DIR}/${FILENAME}_canonized.bcf
rm -f ${TMP_DIR}/${FILENAME}.vcf.gz ${TMP_DIR}/${FILENAME}_canonized.vcf ${INPUT_BCF}; mv ${TMP_DIR}/${FILENAME}_canonized.bcf ${INPUT_BCF}
bcftools index -f ${INPUT_BCF}
END
        chmod +x bnd_canonize.sh
        
        
        
        
        # ---------------------------- Main program ----------------------------

        LocalizeFiles
        rm -f list.txt
        while read -u 3 SAMPLE_ID; do
            echo ./input_files/${SAMPLE_ID}_~{suffix}.bcf >> list.txt
        done 3< ~{sample_ids}

        # Canonizing BNDs before merging. This is necessary, since the 
        # annotation workflow upstream just adds features to the existing 
        # records, without altering the set of records.
        #
        # Remark: canonization should actually be done before annotation, to 
        # save annotation runtime.
        #
        # Remark: this pipeline is still inelegant, since the features of a BND 
        # are not invariant to symmetry, i.e. it is not guaranteed that a BND 
        # record and its symmetrical record get the same score (although they 
        # likely get symmetrical features). Canonization after feature 
        # annotation and filtering is equivalent to keeping a BND if either of 
        # its two records is kept; canonization before feature annotation and 
        # filtering is equivalent to keeping a BND iff its canonical 
        # representation is kept. We should make BND features invariant to
        # symmetry.
        if [[ ~{suffix} == *bnd* || ~{suffix} == *BND* ]]; then
            mkdir ./canonize_dir
            ${TIME_COMMAND} xargs --arg-file=list.txt --max-lines=1 --max-procs=${N_THREADS} ./bnd_canonize.sh ~{docker_dir} ${RAM_PER_THREAD_MB} ./canonize_dir
            tar -czf ./bnd_canonize.tar.gz ./canonize_dir ; rm -rf ./canonize_dir
            gcloud storage mv bnd_canonize.tar.gz ~{remote_outdir}/bnd_canonize.tar.gz
            echo ~{bnd_chromosomes} | tr ',' '\n' > chromosomes.txt
        else
            echo ~{ultralong_chromosomes} | tr ',' '\n' > chromosomes.txt
        fi
        
        # Trivial "hierarchical" bcftools merge with just two steps.
        #
        # Remark: we overwrite bcftools' default merging criterion for INFO/DP
        # since it is `sum`, so it creates very large values whose downstream 
        # utility is unclear. The default for all other values is to copy from
        # the first sample, which is also of unclear utility but at least does
        # not create large values. QUAL is always set to the max.
        #
        # Step 1: merging a few samples at a time over the whole genome.
        split -l ~{n_files_per_merge} -d -a 4 list.txt list_
        N_LIST_FILES=$(ls list_* | wc -l)
        for LIST_FILE in $(ls list_* | sort -V); do
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --info-rules - --file-list ${LIST_FILE} --output-type b --output ${LIST_FILE}_merged.bcf
            xargs --arg-file=${LIST_FILE} --max-lines=1 --max-procs=${N_THREADS} rm -f
            rm -f ${LIST_FILE}
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ${LIST_FILE}_merged.bcf --output ${LIST_FILE}_normed.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_normed.bcf
            rm -f ${LIST_FILE}_merged.bcf*
            ${TIME_COMMAND} xargs --arg-file=chromosomes.txt --max-lines=1 --max-procs=${N_THREADS} ./chunk_by_chr.sh ${LIST_FILE}_normed.bcf
            rm -f ${LIST_FILE}_normed.bcf*
        done
        rm -rf ./input_files/
        
        # Step 2: merging all samples over each chromosome.
        # Only chromosomes with some variant are uploaded.
        while read -u 4 CHROMOSOME; do
            ls ./${CHROMOSOME}/*.bcf | sort -V > list.txt
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --info-rules - --file-list list.txt --output-type b --output ./${CHROMOSOME}/merged.bcf
            N_RECORDS=$(bcftools query --format '%ID\n' ./${CHROMOSOME}/merged.bcf | wc -l)
            if [ ${N_RECORDS} -eq 0 ]; then
                rm -f ./${CHROMOSOME}/merged.bcf
                continue
            fi
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ./${CHROMOSOME}/merged.bcf --output ./${CHROMOSOME}/normed.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ./${CHROMOSOME}/normed.bcf
            gcloud storage mv ./${CHROMOSOME}/normed.bcf ~{remote_outdir}/${CHROMOSOME}.bcf
            gcloud storage mv ./${CHROMOSOME}/normed.bcf.csi ~{remote_outdir}/${CHROMOSOME}.bcf.csi
            rm -rf ./${CHROMOSOME}/
        done 4< chromosomes.txt
        df -h 1>&2
        ls -laht 1>&2
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
