version 1.0


# Given a set of annotated, per-sample, per-SVTYPE query VCFs, and their 
# corresponding canonized svim-asm VCFs and dipcall BEDs, the program selects, 
# for each sample, the query records that are likely true according to the 
# assemblies.
#
# Remark: sequence similarity is not used to compute matches.
#
workflow SV_Integration_UltralongGetTrainingIntervals {
    input {
        File samples_tsv
        String suffix = "del"
        File reference_fai
        
        String remote_indir_annotated
        String remote_indir_truth
        String remote_outdir
        
        Int truvari_refdist = 500
        Float truvari_pctovl = 0
        Float truvari_pctsize = 0.9
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        samples_tsv: "Format: ID, DIPCALL_BED"
        remote_indir_annotated: "Without final slash. Contains per-sample annotated VCFs."
        remote_indir_truth: "Without final slash. Contains per-sample canonized and filtered dipcall VCFs."
        truvari_refdist: "For interval records we set it to a large value to disable the check."
        truvari_pctovl: "For interval records we enable this (it is disabled by default by truvari)."
        truvari_pctsize: "A non-stringent value"
    }
    
    call Impl {
        input:
            samples_tsv = samples_tsv,
            suffix = suffix,
            reference_fai = reference_fai,
            remote_indir_annotated = remote_indir_annotated,
            remote_indir_truth = remote_indir_truth,
            remote_outdir = remote_outdir,
            truvari_refdist = truvari_refdist,
            truvari_pctsize = truvari_pctsize,
            truvari_pctovl = truvari_pctovl,
            docker_image = docker_image
    }
    
    output {
    }
}







# Performance on a 4-core 8GB VM:
# TOOL                                  CPU%        RAM         TIME
#
task Impl {
    input {
        File samples_tsv
        String suffix
        File reference_fai
        
        String remote_indir_query
        String remote_indir_dipcall
        String remote_indir_svimasm
        String remote_outdir
        
        Int svimasm_slack_bp
        Float svimasm_length_similarity

        Int remap_max_length
        Float remap_cov_threshold
        File reference_fa
        File reference_fai

        Int truvari_refdist
        Float truvari_pctsize
        Float truvari_pctovl
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 50
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

        
        function GetTrainingIntervalsThread() {
            local CHUNK_CSV=$1

            local DIPCALL_SVTYPE=""
            local TRUVARI_EXTRA_FLAGS=""
            local SVIMASM_SVTYPE=~{suffix}
            local SVIMASM_SVTYPE=${SVIMASM_SVTYPE^^}


            
------> Allow these steps on complex dups to be toggled from the WDL, to see if adding these events decreases performance.            
            
            
            while read -u 3 LINE; do
                local SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
                local DIPCALL_BED=$(echo ${LINE} | cut -d , -f 2)
                
                # Skipping the sample if it has already been processed or if it
                # was not annotated.
                local TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "1" )
                if [ ${TEST} != "1" ]; then
                    continue
                fi
                local TEST=$( gsutil ls ~{remote_indir_query}/${SAMPLE_ID}_~{suffix}.vcf.gz || echo "1" )
                if [ ${TEST} == "1" ]; then
                    echo "WARNING: sample ${SAMPLE_ID} was not annotated?! Skipping." 1>&2
                    continue
                fi

                # 0. Extracting gaps from the dipcall BED
                # Remark: we do not handle gaps in the reference explicitly,
                # since we assume that calls in reference gaps have already been 
                # removed from the query VCFs upstream.
                gcloud storage cp ${DIPCALL_BED} ./${SAMPLE_ID}_dipcall.bed
                ${TIME_COMMAND} bedtools complement -L -i ${SAMPLE_ID}_dipcall.bed -g ~{reference_fai} > ${SAMPLE_ID}_gaps.bed
                rm -f ${SAMPLE_ID}_dipcall.bed
                local N_GAPS=$(wc -l < ${SAMPLE_ID}_gaps.bed)
                echo "The dipcall BED of ${SAMPLE_ID} has ${N_GAPS} gaps" 1>&2

                # 1. Downloading and filtering the truth VCF
                gcloud storage cp ~{remote_indir_svimasm}/${SAMPLE_ID}_canonized.vcf.gz ./${SAMPLE_ID}_svimasm.vcf.gz
                gcloud storage cp ~{remote_indir_svimasm}/${SAMPLE_ID}_canonized.vcf.gz.tbi ./${SAMPLE_ID}_svimasm.vcf.gz.tbi

                # 1.1 DEL
                bcftools filter --include 'SVTYPE=="DEL"' --output-type z ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_svimasm_del.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_del.vcf.gz

                # 1.2 DUP
                bcftools filter --include 'SVTYPE=="DUP" || SVTYPE=="DUP:TANDEM" || SVTYPE=="DUP:INT"' --output-type v ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_out.vcf
                java -cp ~{docker_dir} UltralongForceDup ${SAMPLE_ID}_out.vcf | bgzip > ${SAMPLE_ID}_svimasm_dup.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_dup.vcf.gz

                # 1.3 INS
                # Even though svim-asm explicitly emits DUPs, it might also 
                # emit INS records that are actually DUP. We reclassify such 
                # INS as DUP as follows.
                #
                # We reclassify as DUP every INS that is compatible with
                # a dipcall BED gap. Recall that a base is included in the 
                # dipcall BED iff: (1) it is covered by one alignment >=50kb 
                # mapQ>=5 in each haplotype; and (2) it is not covered by other 
                # >=10kb alignments in either parent. A simple DUP >10kb might 
                # satisfy (1) but it likely does not satisfy (2). This heuristic
                # is very fast.
                #
                # However, a complex DUP might satisfy (2), since its sub-
                # intervals might be covered by alignments <10kb, so the 
                # heuristic above would only mark simple DUPs as true. If 
                # svim-asm DUP records do not contain complex DUPs, then the 
                # latter would never be marked as true and the model would never
                # learn their BAM patterns.
                #
                # Thus, we use `truvari anno remap` on the remaining INS..................
                # Working only on the remaining INS reduces the runtime of this
                # slow step. However, using two methods might give rise to two
                # populations of DUPs with different breakpoint decisions. Since
                # these two populations are likely simple and complex DUPs, 
                # which are likely to have distinct BAM patterns anyway, we
                # accept this issue in exchange for the advantage in runtime.
                bcftools filter --include 'SVTYPE=="INS"' --output-type v ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_svimasm_ins.vcf
                local N_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_svimasm_ins.vcf | wc -l)
                if [ ${N_INS} -gt 0 ]; then
                    # 1.2.1 Filtering all INS with the dipcall BED
                    ${TIME_COMMAND} java -cp ~{docker_dir} UltralongSvimasmInsExtractDups ${SAMPLE_ID}_svimasm_ins.vcf ${SAMPLE_ID}_gaps.bed $((wc -l < ${SAMPLE_ID}_gaps.bed)) ~{svimasm_slack_bp} ~{svimasm_length_similarity} ${SAMPLE_ID}_svimasm_ins_dup.vcf ${SAMPLE_ID}_svimasm_ins_ins.vcf
                    rm -f ${SAMPLE_ID}_svimasm_ins.vcf ${SAMPLE_ID}_gaps.bed
                    local N_INS_DUP=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_svimasm_ins_dup.vcf | wc -l)
                    if [ ${N_INS_DUP} -gt 0 ]; then
                        ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_svimasm_ins_dup.vcf --output ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                        rm -f ${SAMPLE_ID}_svimasm_ins_dup.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                        ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type z ${SAMPLE_ID}_svimasm_dup.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
                        rm -f ${SAMPLE_ID}_svimasm_dup.vcf.gz* ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_svimasm_dup.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_dup.vcf.gz
                    fi
                    rm -f ${SAMPLE_ID}_svimasm_ins.vcf ; mv ${SAMPLE_ID}_svimasm_ins_ins.vcf ${SAMPLE_ID}_svimasm_ins.vcf

                    # 1.2.2 Filtering the remaining INS with truvari anno remap
                    local N_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_svimasm_ins.vcf | wc -l)
                    if [ ${N_INS} -gt 0 ]; then
                        bgzip --compress-level 1 ${SAMPLE_ID}_svimasm_ins.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins.vcf.gz
                        ${TIME_COMMAND} truvari anno remap --threads ${N_THREADS} --aligner minimap2 --min-length 1 --max-length ~{remap_max_length} --cov-threshold ~{remap_cov_threshold} -r ~{reference_fa} ${SAMPLE_ID}_svimasm_ins.vcf.gz -o ${SAMPLE_ID}_svimasm_ins_remap.vcf.gz
                        rm -f ${SAMPLE_ID}_svimasm_ins.vcf.gz* ; gunzip ${SAMPLE_ID}_svimasm_ins_remap.vcf.gz
                        ${TIME_COMMAND} java -cp ~{docker_dir} UltralongSvimasmInsExtractDupsPrime ${SAMPLE_ID}_svimasm_ins_remap.vcf ${SAMPLE_ID}_svimasm_ins_dup.vcf ${SAMPLE_ID}_svimasm_ins_ins.vcf
                        rm -f ${SAMPLE_ID}_svimasm_ins_remap.vcf
                        bgzip --compress-level 1 --stdout ${SAMPLE_ID}_svimasm_ins_ins.vcf > ${SAMPLE_ID}_svimasm_ins.vcf.gz
                        bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins.vcf.gz
                        ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_svimasm_ins_dup.vcf --output ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                        bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                        ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type z ${SAMPLE_ID}_svimasm_dup.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
                        rm -f ${SAMPLE_ID}_svimasm_dup.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_svimasm_dup.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_dup.vcf.gz
                    fi
                fi

                # 1.4 INV
                bcftools filter --include 'SVTYPE=="INV"' --output-type z ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_svimasm_inv.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_inv.vcf.gz

                rm -f ${SAMPLE_ID}_svimasm.vcf.gz*

                # 2. Downloading and filtering the query VCF
                gcloud storage cp ~{remote_indir_query}/${SAMPLE_ID}_~{suffix}.vcf.gz ./${SAMPLE_ID}_query.vcf.gz
                gcloud storage cp ~{remote_indir_query}/${SAMPLE_ID}_~{suffix}.vcf.gz.tbi ./${SAMPLE_ID}_query.vcf.gz.tbi
                bcftools filter --include 'SVTYPE=="DEL"' --output-type z ${SAMPLE_ID}_query.vcf.gz --output ${SAMPLE_ID}_del.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_del.vcf.gz
                bcftools filter --include 'SVTYPE=="DUP"' --output-type z ${SAMPLE_ID}_query.vcf.gz --output ${SAMPLE_ID}_dup.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_dup.vcf.gz
                bcftools filter --include 'SVTYPE=="INS"' --output-type z ${SAMPLE_ID}_query.vcf.gz --output ${SAMPLE_ID}_ins.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins.vcf.gz
                bcftools filter --include 'SVTYPE=="INV"' --output-type z ${SAMPLE_ID}_query.vcf.gz --output ${SAMPLE_ID}_inv.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_inv.vcf.gz

                # 3. Computing matches
                # Remark: for speed reasons, sequence similarity is never used.
                #
                # 3.1 DEL
                # Remark: svim-asm does not emit calls at some query DELs that 
                # correspond to gaps in the dipcall BED. DELs >10kb are likely
                # to create gaps in dipcall's BED by the definition of confident
                # BED, since such DELs induce split alignments. We choose to
                # mark as true every query DEL that matches a gap in dipcall's
                # BED.
                ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_svimasm_del.vcf.gz -c ${SAMPLE_ID}_del.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
                mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_del1.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_del1.vcf.gz
                rm -rf ${SAMPLE_ID}_truvari/
                bcftools view --header-only ${SAMPLE_ID}_del.vcf.gz > ${SAMPLE_ID}_gaps.vcf
                ${TIME_COMMAND} java -cp ~{docker_dir} UltralongBed2IntervalVcf ${SAMPLE_ID}_gaps.bed DEL >> ${SAMPLE_ID}_gaps.vcf
                bgzip -@ ${N_THREADS} ${SAMPLE_ID}_gaps.vcf
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_gaps.vcf.gz
                ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_gaps.vcf.gz -c ${SAMPLE_ID}_del.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
                mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_del2.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_del2.vcf.gz
                rm -rf ${SAMPLE_ID}_truvari/
                bcftools concat --allow-overlaps --remove-duplicates --output-type z ${SAMPLE_ID}_del1.vcf.gz ${SAMPLE_ID}_del2.vcf.gz --output ${SAMPLE_ID}_del_training.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_del_training.vcf.gz
                rm -f ${SAMPLE_ID}_del1.vcf.gz* ${SAMPLE_ID}_del2.vcf.gz* ${SAMPLE_ID}_gaps.vcf.gz*

                # 3.2 DUP
                # We observe that DUP records from svim-asm are generally 
                # comprehensive, and they do not need to be augmented with e.g.
                # gaps in dipcall's BED.
                ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_svimasm_dup.vcf.gz -c ${SAMPLE_ID}_dup.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
                mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_dup_training.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_dup_training.vcf.gz
                
                # 3.3 INS
                # Remark: we don't use `--dup-to-ins` in truvari, since we 
                # assume that DUP and INS calls are already accurately
                # classified in both the query and the svim-asm VCF, so we can
                # compare each SVTYPE in isolation.
                ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_svimasm_ins.vcf.gz -c ${SAMPLE_ID}_ins.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
                mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_ins_training.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins_training.vcf.gz

                # 3.4 INV
                # Remark: svim-asm does not emit calls at some query INVs that 
                # correspond to gaps in the dipcall BED. Since INV records are 
                # rare, we choose to mark as true such query INVs, at the risk
                # of polluting the training set with complex events that are not
                # just INVs.
                ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_svimasm_inv.vcf.gz -c ${SAMPLE_ID}_inv.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
                mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_inv1.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_inv1.vcf.gz
                rm -rf ${SAMPLE_ID}_truvari/
                bcftools view --header-only ${SAMPLE_ID}_inv.vcf.gz > ${SAMPLE_ID}_gaps.vcf
                ${TIME_COMMAND} java -cp ~{docker_dir} UltralongBed2IntervalVcf ${SAMPLE_ID}_gaps.bed INV >> ${SAMPLE_ID}_gaps.vcf
                bgzip -@ ${N_THREADS} ${SAMPLE_ID}_gaps.vcf
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_gaps.vcf.gz
                ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_gaps.vcf.gz -c ${SAMPLE_ID}_inv.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
                mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_inv2.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_inv2.vcf.gz
                rm -rf ${SAMPLE_ID}_truvari/
                bcftools concat --allow-overlaps --remove-duplicates --output-type z ${SAMPLE_ID}_inv1.vcf.gz ${SAMPLE_ID}_inv2.vcf.gz --output ${SAMPLE_ID}_inv_training.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_inv_training.vcf.gz
                rm -f ${SAMPLE_ID}_inv1.vcf.gz* ${SAMPLE_ID}_inv2.vcf.gz* ${SAMPLE_ID}_gaps.vcf.gz*
                
                # Uploading
                gcloud storage mv ${SAMPLE_ID}_'*_training.vcf.gz*' ~{remote_outdir}
            done 3< ${CHUNK_CSV}
        }



        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        samtools --version 1>&2
        bcftools --version 1>&2
        truvari --help 1>&2
        df -h 1>&2

        CHUNKSIZE_FLAG=""
        if [ ~{truvari_refdist} -gt 1000 ]; then
            # To avoid ERROR:root:--chunksize must be >= --refdist
            CHUNKSIZE_FLAG="--chunksize ~{truvari_refdist}"
        fi
        
        cat ~{samples_tsv} | tr '\t' ',' > samples.csv
        N_ROWS=$(wc -l < samples.csv)
        if [ ${N_ROWS} -gt ${N_THREADS} ]; then
            N_ROWS_PER_THREAD=$(( ${N_ROWS} / ${N_THREADS} ))
            split -l ${N_ROWS_PER_THREAD} -d -a 4 samples.csv chunk_
        else
            mv samples.csv chunk_0
        fi
        for FILE in $(ls chunk_*); do
            GetTrainingIntervalsThread ${FILE} &
        done
        wait
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




        # Old GetTrainingIntervals() function from UltralongAnnotate.wdl
        #
        # Given an interval-only input VCF, the procedure compresses it and
        # computes its records with a stringent match to the training resource.
        #
        # Remark: sequence similarity is not used to decide a match.
        #
#        function GetTrainingIntervals() {
#            local SAMPLE_ID=$1
#            local INPUT_VCF=$2
#            local TRUTH_VCF_GZ=$3
#            local SUFFIX=$4
#            
#            bcftools view --output-type z ${INPUT_VCF} --output ${SAMPLE_ID}_${SUFFIX}.vcf.gz
#            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${SUFFIX}.vcf.gz
#            rm -f ${INPUT_VCF}
#            
#            ${TIME_COMMAND} truvari bench -b ${TRUTH_VCF_GZ} -c ${SAMPLE_ID}_${SUFFIX}.vcf.gz --includebed ~{ultralong_training_resource_bed} --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0 --pick single -o ./truvari_${SAMPLE_ID}/
#            
#            mv truvari_${SAMPLE_ID}/tp-comp.vcf.gz ${SAMPLE_ID}_${SUFFIX}_training.vcf.gz
#            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${SUFFIX}_training.vcf.gz
#            rm -rf truvari_${SAMPLE_ID}/
#        }