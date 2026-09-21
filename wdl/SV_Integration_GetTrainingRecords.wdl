version 1.0


# Isolates into a separate workflow just the `GetTrainingRecords()` function in
# `SV_Integration_Workpackage1.wdl`. 
#
# Purpose: experiment with a smaller training resource when plotting ROC curves.
#
workflow SV_Integration_GetTrainingRecords {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir
        
        File training_resource_vcf_gz
        File training_resource_tbi
        File training_resource_bed
        
        File reference_fai
        File reference_agp
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        sv_integration_chunk_tsv: "A subset of the rows of table `sv_integration_hg38`, without the header."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
        training_resource_vcf_gz: "We assume that the training resource VCF has already been subset to the correct length range upstream."
        training_resource_bed: "Training resource calls can belong only to these regions. Typically a high-confidence dipcall BED, or a BED derived from intersecting multiple dipcall BEDs."
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            
            training_resource_vcf_gz = training_resource_vcf_gz,
            training_resource_tbi = training_resource_tbi,
            training_resource_bed = training_resource_bed,
            
            reference_fai = reference_fai,
            reference_agp = reference_agp,
            
            docker_image = docker_image
    }
    
    output {
    }
}


# Truvari bench (measured on a 6 CPUs, 8GB VM):
#
# Sequential                  10 m
# Parallel                     5 m
#
task Impl {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir
        
        File training_resource_vcf_gz
        File training_resource_tbi
        File training_resource_bed
        
        File reference_fa
        File reference_fai
        File standard_chromosomes_bed
        File autosomes_bed
        File reference_agp
        
        String docker_image
        Int n_cpu = 6
        Int ram_size_gb = 8
        Int disk_size_gb = 10
        Int preemptible_number = 0
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
        function LocalizeSample() {
            local SAMPLE_ID=$1

            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}_kanpig.vcf.'gz*' .
        }
        
        
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -rf ${SAMPLE_ID}_*
        }
        
        
        # Builds a BED file that excludes every gap from the AGP file of
        # the reference.
        #
        function GetReferenceGaps() {
            # Computing non-gap regions
            awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    if ( ( $1=="chr1" || $1=="chr2" || $1=="chr3" || $1=="chr4" || $1=="chr5" || $1=="chr6" || $1=="chr7" || $1=="chr8" || $1=="chr9" || $1=="chr10" || \
                           $1=="chr11" || $1=="chr12" || $1=="chr13" || $1=="chr14" || $1=="chr15" || $1=="chr16" || $1=="chr17" || $1=="chr18" || $1=="chr19" || $1=="chr20" || \
                           $1=="chr21" || $1=="chr22" || $1=="chrX" || $1=="chrY" || $1=="chrM" \
                         ) && $5=="N" \
                       ) print $0 \
                 }' ~{reference_agp} > gaps_unsorted.bed
            bedtools sort -i gaps_unsorted.bed -faidx ~{reference_fai} > gaps.bed
            bedtools complement -L -i gaps.bed -g ~{reference_fai} > not_gaps.bed
            
            # Intersecting non-gap regions with the training BED
            bedtools sort -i ~{training_resource_bed} -faidx ~{reference_fai} > training_resource_sorted.bed
            rm -f training_not_gaps_beds.wsv
            local ID="0"
            local ROW
            while read -u 4 ROW || [ -n "${ROW}" ]; do
                ID=$(( ${ID} + 1 ))
                echo "${ROW}" > ${ID}.bed
                bedtools intersect -a ${ID}.bed -b training_resource_sorted.bed -sorted -g ~{reference_fai} > training_not_gaps_${ID}.bed
                if [ -s training_not_gaps_${ID}.bed ]; then
                    echo "${ID} training_not_gaps_${ID}.bed" >> training_not_gaps_beds.wsv
                else
                    rm -f training_not_gaps_${ID}.bed
                fi
                rm -f ${ID}.bed
            done 4< not_gaps.bed
            ls -lht *.bed 1>&2
            
            # Removing temporary files
            rm -f gaps_unsorted.bed training_resource_sorted.bed
        }
        
        
        cat << 'END' > truvari_bench.sh
#!/bin/bash
set -euxo pipefail

SAMPLE_ID=$1
INPUT_VCF_GZ=$2
TRAINING_RESOURCE_VCF_GZ=$3
INFINITY=$4
CHUNK_ID=$5
INCLUDE_BED=$6

truvari bench -b ${TRAINING_RESOURCE_VCF_GZ} -c ${INPUT_VCF_GZ} --includebed ${INCLUDE_BED} --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick single -o ${SAMPLE_ID}_truvari_${CHUNK_ID}/
END
        chmod +x truvari_bench.sh
        
        
        # Extracts every record that has a stringent `truvari bench` match with
        # some records in the resource.
        #
        # Remark: we use `--pick single` to force every resource record to be
        # matched with at most one sample record, which is hopefully the
        # most similar to it. This is because we assume that using a
        # contaminated training set in XGBoost downstream is worse than using a
        # slightly smaller training set. With `--pick multi` e.g. two records in
        # the sample VCF might be matched to the same record in the resource 
        # VCF (probably not good) and vice versa (good).
        #
        # Remark: multiple instances of `truvari bench` are run in parallel
        # using `not_gaps.bed`.
        #
        # Remark: in few anecdotal tests, `--pick multi` seems a bit faster than
        # `--pick single` (4m vs 5m with 6 hyperthreading cores).
        #
        # Remark: both the inputs and the output of the function are indexed
        # `.vcf.gz`, since they are needed by `truvari bench`.
        #
        function GetTrainingRecords() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            # Running in parallel
            ${TIME_COMMAND} xargs --arg-file=training_not_gaps_beds.wsv --max-lines=1 --max-procs=${N_THREADS} ./truvari_bench.sh ${SAMPLE_ID} ${INPUT_VCF_GZ} ~{training_resource_vcf_gz} ${INFINITY}
            
            # Concatenating outputs
            local ID
            local ROW
            rm -f ${SAMPLE_ID}_outputs.txt
            while read -u 4 ROW || [ -n "${ROW}" ]; do
                ID=$(echo ${ROW} | cut -d ' ' -f 1)
                echo ${SAMPLE_ID}_truvari_${ID}/tp-comp.vcf.gz >> ${SAMPLE_ID}_outputs.txt
            done 4< training_not_gaps_beds.wsv
            ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list ${SAMPLE_ID}_outputs.txt --output-type z --output ${SAMPLE_ID}_training.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_training.vcf.gz
            
            # Removing temporary files
            rm -rf ${SAMPLE_ID}_outputs.txt ./${SAMPLE_ID}_truvari_*/
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        truvari --help 1>&2
        
        GetReferenceGaps 
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE || [ -n "${LINE}" ]; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            LocalizeSample ${SAMPLE_ID}
            GetTrainingRecords ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz
            
            # Uploading
            gcloud storage mv ${SAMPLE_ID}_training.vcf.'gz*' ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done 3< chunk.csv
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
