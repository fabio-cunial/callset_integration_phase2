version 1.0


# Simple evaluations of ultralong calls. BND calls are analyzed by a separate
# workflow.
# 
# Structure of `remote_outdir`:
#
# ├── dipcall/                                                  for each sample;
# │   └── type_X/
# │       └── length_X/
# ├── samples/                                  present records for each sample;
# │   └── type_X/
# │       └── length_X/
# └── precision_recall/                                         output measures;
#     └── type_X/
#         └── length_X/
#
workflow SV_Integration_UltralongAnalysis {
    input {
        String remote_workpackage_11_dir
        String remote_outdir
        
        String precision_recall_samples_csv
        Int sequence_similarity
        Int limit_to_dipcall_bed
        Int min_sv_length = 10000
        
        Array[String] sv_types = ["DEL", "INS"]
        Array[Int] sv_length_bins = [10000,50000,100000,500000,1000000]
        String sv_length_bins_str = "10000,50000,100000,500000,1000000"
        Int n_length_bins = 5
        
        File reference_fa
        File reference_fai
        File reference_agp
        File tandem_bed
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages:latest"
        Int preemptible_number = 4
    }
    parameter_meta {
        precision_recall_samples_csv: "Format: ID,dipcall_bed_uri,dipcall_vcf_uri"
        sequence_similarity: "0=OFF, 1=ON."
        limit_to_dipcall_bed: "0=NO, 1=YES."
        sv_length_bins_str: "Comma-separated"
    }
    
    # Preparing the VCFs to be benchmarked
    call ComplementBed {
        input:
            tandem_bed = tandem_bed,
            reference_fai = reference_fai,
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    call SplitBcfBySample {
        input:
            samples_csv = precision_recall_samples_csv,
            remote_workpackage_11_dir = remote_workpackage_11_dir,
            remote_outdir = remote_outdir + "/samples",
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    call CanonizeDipcall {
        input:
            samples_csv = precision_recall_samples_csv,
            min_sv_length = min_sv_length,
            reference_agp = reference_agp,
            limit_to_dipcall_bed = limit_to_dipcall_bed,
            remote_outdir = remote_outdir + "/dipcall",
            reference_fai = reference_fai,
            reference_agp = reference_agp,
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }    
    
    # For each sample: global analysis.
    call PrecisionRecallAnalysis {
        input:
            samples_csv = precision_recall_samples_csv,
            remote_indir_samples = remote_outdir + "/samples",
            remote_indir_dipcall = remote_outdir + "/dipcall",
            remote_outdir = remote_outdir + "/precision_recall",
        
            sequence_similarity = sequence_similarity,
            limit_to_dipcall_bed = limit_to_dipcall_bed,
            min_sv_length = min_sv_length,
        
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            tandem_bed = ComplementBed.sorted_bed,
            not_tandem_bed = ComplementBed.complement_bed,
            
            in_flag = [CanonizeDipcall.out_flag, SplitBcfBySample.out_flag],
            
            docker_image = docker_image,
            preemptible_number = preemptible_number
    }
    
    # For each sample: analysis by type and length.
    scatter (i in range(length(sv_types))) {
        call FilterByType as by_type {
            input:
                samples_csv = precision_recall_samples_csv,
                sv_type = sv_types[i],
                remote_indir = remote_outdir + "/samples",
                remote_outdir = remote_outdir + "/samples/type_" + sv_types[i],
                in_flag = [SplitBcfBySample.out_flag],
                docker_image = docker_image,
                preemptible_number = preemptible_number
        }
        call FilterByType as by_type_dipcall {
            input:
                samples_csv = precision_recall_samples_csv,
                sv_type = sv_types[i],
                remote_indir = remote_outdir + "/dipcall",
                remote_outdir = remote_outdir + "/dipcall/type_" + sv_types[i],
                in_flag = [CanonizeDipcall.out_flag],
                docker_image = docker_image,
                preemptible_number = preemptible_number
        }
        call PrecisionRecallAnalysis as pr_analysis_type {
            input:
                samples_csv = precision_recall_samples_csv,
                remote_indir_samples = remote_outdir + "/samples/type_" + sv_types[i],
                remote_indir_dipcall = remote_outdir + "/dipcall/type_" + sv_types[i],
                remote_outdir = remote_outdir + "/precision_recall/type_" + sv_types[i],
    
                sequence_similarity = sequence_similarity,
                limit_to_dipcall_bed = limit_to_dipcall_bed,
                min_sv_length = min_sv_length,
    
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                tandem_bed = ComplementBed.sorted_bed,
                not_tandem_bed = ComplementBed.complement_bed,
        
                in_flag = [by_type.out_flag, by_type_dipcall.out_flag],
        
                docker_image = docker_image,
                preemptible_number = preemptible_number
        }
        scatter (j in range(n_length_bins)) {
            call FilterByLength as by_length {
                input:
                    samples_csv = precision_recall_samples_csv,
                    sv_length_bins = sv_length_bins_str,
                    index = j,
                    remote_indir = remote_outdir + "/samples/type_" + sv_types[i],
                    remote_outdir = remote_outdir + "/samples/type_" + sv_types[i] + "/length_" + sv_length_bins[j],
                    in_flag = [by_type.out_flag],
                    docker_image = docker_image,
                    preemptible_number = preemptible_number
            }
            call FilterByLength as by_length_dipcall {
                input:
                    samples_csv = precision_recall_samples_csv,
                    sv_length_bins = sv_length_bins_str,
                    index = j,
                    remote_indir = remote_outdir + "/dipcall/type_" + sv_types[i],
                    remote_outdir = remote_outdir + "/dipcall/type_" + sv_types[i] + "/length_" + sv_length_bins[j],
                    in_flag = [by_type_dipcall.out_flag],
                    docker_image = docker_image,
                    preemptible_number = preemptible_number
            }
            call PrecisionRecallAnalysis as pr_analysis_type_length {
                input:
                    samples_csv = precision_recall_samples_csv,
                    remote_indir_samples = remote_outdir + "/samples/type_" + sv_types[i] + "/length_" + sv_length_bins[j],
                    remote_indir_dipcall = remote_outdir + "/dipcall/type_" + sv_types[i] + "/length_" + sv_length_bins[j],
                    remote_outdir = remote_outdir + "/precision_recall/type_" + sv_types[i] + "/length_" + sv_length_bins[j],
        
                    sequence_similarity = sequence_similarity,
                    limit_to_dipcall_bed = limit_to_dipcall_bed,
                    min_sv_length = min_sv_length,
        
                    reference_fa = reference_fa,
                    reference_fai = reference_fai,
                    tandem_bed = ComplementBed.sorted_bed,
                    not_tandem_bed = ComplementBed.complement_bed,
            
                    in_flag = [by_length.out_flag,by_length_dipcall.out_flag],
            
                    docker_image = docker_image,
                    preemptible_number = preemptible_number
            }
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
        
        String docker_image
        Int n_cpu = 1
        Int ram_size_gb = 4
        Int preemptible_number
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
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}


#
task CanonizeDipcall {
    input {
        File samples_csv
        Int min_sv_length
        File reference_agp
        Int limit_to_dipcall_bed
        
        String remote_outdir
        
        File reference_fai
        File reference_agp
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 50
        Int preemptible_number
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
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        
        
        
        
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
        
        
        # Puts in canonical form a raw VCF from dipcall. This is similar to
        # `SV_Integration_BuildTrainingResource.wdl`.
        #
        function CanonizeDipcallVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local INPUT_TBI=$3
            local DIPCALL_BED=$4
            local MIN_SV_LENGTH=$5
            local MAX_SV_LENGTH=$6
            local NOT_GAPS_BED=$7
            
            
            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_in.vcf.gz
            mv ${INPUT_TBI} ${SAMPLE_ID}_in.vcf.gz.tbi
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics - --output-type b ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.bcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.bcf ${SAMPLE_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_in.bcf
            
            # Removing SNVs, records with unresolved REF/ALT, records that are
            # not marked as present, and records with a FILTER. 
            ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || (GT!="alt" && GT!=".|1" && GT!="1|." && GT!="./1" && GT!="1/.") || (FILTER!="PASS" && FILTER!=".") || REF="*" || ALT="*"' --output-type b ${SAMPLE_ID}_in.bcf --output ${SAMPLE_ID}_out.bcf
            rm -f ${SAMPLE_ID}_in.bcf* ; mv ${SAMPLE_ID}_out.bcf ${SAMPLE_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_in.bcf
            
            # Removing records in reference gaps
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.bcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.bcf* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the dipcall BED, if needed.
            if [ ~{limit_to_dipcall_bed} -eq 1 ]; then
                ${TIME_COMMAND} bcftools filter --regions-file ${DIPCALL_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            fi
            
            # Making sure SVLEN and SVTYPE are consistently annotated
            truvari anno svinfo --minsize 1 ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the given length range
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type b ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.bcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.bcf ${SAMPLE_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_in.bcf
            
            mv ${SAMPLE_ID}_in.bcf ${SAMPLE_ID}.bcf
            mv ${SAMPLE_ID}_in.bcf.csi ${SAMPLE_ID}.bcf.csi
        }
        
        
        
        
        # --------------------------- Main program -----------------------------
        
        MAX_SV_LENGTH="3000000000"  # Arbitrary
        GetReferenceGaps ~{reference_agp} not_gaps.bed
        while read ROW; do
            SAMPLE_ID=$( echo ${ROW} | cut -d , -f 1 )
            REMOTE_BED=$( echo ${ROW} | cut -d , -f 2 )
            REMOTE_VCF=$( echo ${ROW} | cut -d , -f 3 )
            
            # Skipping the sample if it has already been canonized
            TEST=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            # Canonizing
            gcloud storage cp ${REMOTE_BED} ./${SAMPLE_ID}.bed
            gcloud storage cp ${REMOTE_VCF} ./${SAMPLE_ID}.vcf.gz
            bcftools index --threads ${N_THREADS} ${SAMPLE_ID}.vcf.gz
            CanonizeDipcallVcf ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz ${SAMPLE_ID}.vcf.gz.tbi ${SAMPLE_ID}.bed ~{min_sv_length} ${MAX_SV_LENGTH} not_gaps.bed
            
            # Uploading
            gcloud storage mv ${SAMPLE_ID}.'bcf*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            rm -f ${SAMPLE_ID}*
        done < ~{samples_csv}
        
        # Fake output
        echo "done" > out.txt
    >>>
    
    output {
        File out_flag = "out.txt"
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


# Writes to a separate file every sample column.
#
# Remark: we keep only records that are genotyped as present in each sample.
#
# Performance on 12'680 samples, whole genome, SSD, 10 samples:
#
# TOOL                              CPU     RAM      TIME
# bcftools +split ultralong        100%      2G     1h40m
# bcftools +split bnd              100%     85M       10s
#
task SplitBcfBySample {
    input {
        File samples_csv
        
        String remote_workpackage_11_dir
        String remote_outdir
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 50
        Int preemptible_number
    }
    parameter_meta {
        remote_outdir: "The result of the split is stored in this bucket location."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        
        gcloud storage cp ~{remote_workpackage_11_dir}/truvari_collapsed.'bcf*' .
        cut -d , -f 1 ~{samples_csv} | sort | uniq > samples.txt
        ${TIME_COMMAND} bcftools +split --samples-file samples.txt --output-type b --output . truvari_collapsed.bcf
        rm -f truvari_collapsed.bcf*
        for FILE in $(ls *.bcf); do
            ${TIME_COMMAND} bcftools filter --include 'GT="alt" || GT=".|1" || GT="1|." || GT="./1" || GT="1/."' --output-type b ${FILE} --output out.bcf
            rm -f ${FILE}* ; mv out.bcf ${FILE} ; bcftools index --threads ${N_THREADS} -f ${FILE}
        done
        ls -laht
        gcloud storage mv '*.bcf*' ~{remote_outdir}/        
        
        # Fake output
        echo "done" > out.txt
    >>>
    
    output {
        File out_flag = "out.txt"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}


#
task FilterByType {
    input {
        File samples_csv
        String sv_type
        
        String remote_indir
        String remote_outdir
        
        Array[File] in_flag
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 50
        Int preemptible_number
    }
    parameter_meta {
        sv_type: "INS or DEL only, since dipcall does not report DUP,INV."
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        
        cut -d , -f 1 ~{samples_csv} | sort | uniq > samples.txt
        while read SAMPLE_ID; do
            # Skipping the sample if it has already been filtered
            TEST=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            # Filtering
            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}'.bcf*' .
            if [ ~{sv_type} = "INS" ]; then
                ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'SVTYPE="INS" || SVTYPE="DUP"' --output-type b ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}_~{sv_type}.bcf
            else
                ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'SVTYPE="'~{sv_type}'"' --output-type b ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}_~{sv_type}.bcf
            fi
            
            rm -f ${SAMPLE_ID}.bcf* ; mv ${SAMPLE_ID}_~{sv_type}.bcf ${SAMPLE_ID}.bcf ; bcftools index --threads ${N_THREADS} ${SAMPLE_ID}.bcf
            
            # Uploading
            gcloud storage mv ${SAMPLE_ID}'.bcf*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
        done < samples.txt
        
        # Fake output
        echo "done" > out.txt
    >>>
    
    output {
        File out_flag = "out.txt"
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


#
task FilterByLength {
    input {
        File samples_csv
        String sv_length_bins
        Int index
        
        String remote_indir
        String remote_outdir
        
        Array[File] in_flag
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 50
        Int preemptible_number
    }
    parameter_meta {
        sv_length_bins: "Comma-separated, increasing."
        index: "Zero-based"
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        
        SVLEN_MAX=$( echo ~{sv_length_bins} | cut -d , -f $(( ~{index} + 1 )) )
        if [ ~{index} -eq 0 ]; then
            SVLEN_MIN="0"
        else
            SVLEN_MIN=$( echo ~{sv_length_bins} | cut -d , -f ~{index} )
        fi
        cut -d , -f 1 ~{samples_csv} | sort | uniq > samples.txt
        while read SAMPLE_ID; do
            # Skipping the sample if it has already been filtered
            TEST=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
        
            # Filtering
            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}.'bcf*' .    
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${SVLEN_MIN}' && ABS(SVLEN)<'${SVLEN_MAX} --output-type b ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}_${SVLEN_MAX}.bcf
            rm -f ${SAMPLE_ID}.bcf* ; mv ${SAMPLE_ID}_${SVLEN_MAX}.bcf ${SAMPLE_ID}.bcf ; bcftools index --threads ${N_THREADS} ${SAMPLE_ID}.bcf
            
            # Uploading
            gcloud storage mv ${SAMPLE_ID}'.bcf*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
        done < samples.txt
        
        # Fake output
        echo "done" > out.txt
    >>>
    
    output {
        File out_flag = "out.txt"
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


# Performance with 4 cores and 32GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# truvari bench             
# vcfdist
#
task PrecisionRecallAnalysis {
    input {
        File samples_csv
        String remote_indir_samples
        String remote_indir_dipcall
        String remote_outdir
        
        Int min_sv_length
        Int sequence_similarity
        Int limit_to_dipcall_bed
        
        File reference_fa
        File reference_fai
        File tandem_bed
        File not_tandem_bed
        
        Array[File] in_flag
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 4
        Int disk_size_gb = 20
        Int preemptible_number
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
        
        
        MAX_SV_LENGTH="3000000000"  # Arbitrary
        if [ ~{sequence_similarity} -eq 0 ]; then
            TRUVARI_MATCH_FLAGS="--pctseq 0 --dup-to-ins"
        elif [ ~{sequence_similarity} -eq 1 ]; then
            TRUVARI_MATCH_FLAGS="--reference ~{reference_fa} --max-resolve ${MAX_SV_LENGTH} --dup-to-ins"
        elif [ ~{sequence_similarity} -eq 2 ]; then
            TRUVARI_MATCH_FLAGS="--pctseq 0 --pick multi"
        fi
        if [ ~{limit_to_dipcall_bed} -eq 1 ]; then
            TRUVARI_BED_FLAG="--includebed ${DIPCALL_BED}"
        else
            TRUVARI_BED_FLAG=" "
        fi
        while read ROW; do
            SAMPLE_ID=$(echo ${ROW} | cut -d , -f 1)
            DIPCALL_BED_URI=$(echo ${ROW} | cut -d , -f 2)
            
            # Skipping the sample if it has already been processed
            TEST=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            # Downloading
            gcloud storage cp ${DIPCALL_BED_URI} ${SAMPLE_ID}_dipcall.bed
            gcloud storage cp ~{remote_indir_samples}/${SAMPLE_ID}.'bcf*' .
            gcloud storage cp ~{remote_indir_dipcall}/${SAMPLE_ID}.bcf ./${SAMPLE_ID}_dipcall.bcf
            gcloud storage cp ~{remote_indir_dipcall}/${SAMPLE_ID}.bcf.csi ./${SAMPLE_ID}_dipcall.bcf.csi
            
            # Extracting calls with POS inside and outside TRs
            ${TIME_COMMAND} bcftools view --output-type z ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}.vcf.gz &
            ${TIME_COMMAND} bcftools view --output-type z ${SAMPLE_ID}_dipcall.bcf --output ${SAMPLE_ID}_dipcall.vcf.gz &
            wait
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}.vcf.gz &
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_dipcall.vcf.gz &
            wait
            ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}_tr.vcf.gz &
            ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${SAMPLE_ID}.bcf --output ${SAMPLE_ID}_not_tr.vcf.gz &
            ${TIME_COMMAND} bcftools view --regions-file ~{tandem_bed} --regions-overlap pos --output-type z ${SAMPLE_ID}_dipcall.bcf --output ${SAMPLE_ID}_dipcall_tr.vcf.gz &
            ${TIME_COMMAND} bcftools view --regions-file ~{not_tandem_bed} --regions-overlap pos --output-type z ${SAMPLE_ID}_dipcall.bcf --output ${SAMPLE_ID}_dipcall_not_tr.vcf.gz &
            wait
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_tr.vcf.gz &
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_not_tr.vcf.gz &
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_dipcall_tr.vcf.gz &
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_dipcall_not_tr.vcf.gz &
            wait
            
            # Benchmarking
            rm -rf ./truvari_*
            ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_dipcall.vcf.gz        -c ${SAMPLE_ID}.vcf.gz        ${TRUVARI_MATCH_FLAGS} ${TRUVARI_BED_FLAG} --sizemin 1 --sizefilt 1 --sizemax ${MAX_SV_LENGTH} -o ./truvari_all/ &
            ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_dipcall_tr.vcf.gz     -c ${SAMPLE_ID}_tr.vcf.gz     ${TRUVARI_MATCH_FLAGS} ${TRUVARI_BED_FLAG} --sizemin 1 --sizefilt 1 --sizemax ${MAX_SV_LENGTH} -o ./truvari_tr/ &
            ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_dipcall_not_tr.vcf.gz -c ${SAMPLE_ID}_not_tr.vcf.gz ${TRUVARI_MATCH_FLAGS} ${TRUVARI_BED_FLAG} --sizemin 1 --sizefilt 1 --sizemax ${MAX_SV_LENGTH} -o ./truvari_not_tr/ &
            wait
            mv ./truvari_all/summary.json ./${SAMPLE_ID}_all.txt
            mv ./truvari_tr/summary.json ./${SAMPLE_ID}_tr.txt
            mv ./truvari_not_tr/summary.json ./${SAMPLE_ID}_not_tr.txt
        
            # Uploading
            gcloud storage mv ${SAMPLE_ID}_'*.txt' ~{remote_outdir}/
            echo "done" > ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            rm -f ${SAMPLE_ID}*
        done < ~{samples_csv}
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
