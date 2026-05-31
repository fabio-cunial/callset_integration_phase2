version 1.0


# Given a set of annotated, per-sample, per-SVTYPE query VCFs, and their 
# corresponding canonized svim-asm VCFs and dipcall BEDs, the program selects, 
# for each sample, the query records that are likely true according to the 
# assemblies.
#
# Remark: sequence similarity is never used to compute matches.
#
workflow SV_Integration_UltralongGetTrainingIntervals {
    input {
        File samples_csv
        
        String remote_indir_query
        String remote_indir_svimasm
        String remote_outdir
        
        Int convert_svimasm_ins_to_dup = 1

        Int svimasm_ins_use_gaps = 1
        Int svimasm_ins_use_gaps_slack_bp = 200
        Float svimasm_ins_use_gaps_length_similarity = 0.9

        Int svimasm_ins_use_remap = 1
        Int svimasm_ins_remap_max_length = 2000000
        Float svimasm_ins_remap_cov_threshold = 0.8

        Int match_to_gaps = 1
        Int match_insdups_to_dups = 1

        File reference_fa
        File reference_fai
        
        Int truvari_refdist = 500
        Float truvari_pctsize = 0.9
        Float truvari_pctovl = 0
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong_remap:latest"
    }
    parameter_meta {
        samples_csv: "Format: ID, DIPCALL_BED"
        remote_indir_query: "Without final slash. Contains per-sample annotated VCFs created by `SV_Integration_UltralongAnnotate.wdl`."
        remote_indir_svimasm: "Without final slash. Contains per-sample canonized and filtered svim-asm VCFs."
        convert_svimasm_ins_to_dup: "1=SVIM-asm's INS records that correspond to duplications are rewritten as DUP records"
        truvari_refdist: "Should be set for stringent matches, since position is crucial for detecting BAM patterns."
        truvari_pctovl: "Should be set for stringent matches, since position is crucial for detecting BAM patterns."
        truvari_pctsize: "Should be set for stringent matches, since position is crucial for detecting BAM patterns."
    }
    
    call Impl {
        input:
            samples_csv = samples_csv,

            remote_indir_query = remote_indir_query,
            remote_indir_svimasm = remote_indir_svimasm,
            remote_outdir = remote_outdir,

            convert_svimasm_ins_to_dup = convert_svimasm_ins_to_dup,

            svimasm_ins_use_gaps = svimasm_ins_use_gaps,
            svimasm_ins_use_gaps_slack_bp = svimasm_ins_use_gaps_slack_bp,
            svimasm_ins_use_gaps_length_similarity = svimasm_ins_use_gaps_length_similarity,

            svimasm_ins_use_remap = svimasm_ins_use_remap,
            svimasm_ins_remap_max_length = svimasm_ins_remap_max_length,
            svimasm_ins_remap_cov_threshold = svimasm_ins_remap_cov_threshold,

            match_to_gaps = match_to_gaps,
            match_insdups_to_dups = match_insdups_to_dups,

            reference_fa = reference_fa,
            reference_fai = reference_fai,

            truvari_refdist = truvari_refdist,
            truvari_pctsize = truvari_pctsize,
            truvari_pctovl = truvari_pctovl,

            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 2-core 8GB VM:
#
# TOOL                                      CPU%        RAM         TIME
# UltralongSvimasmInsExtractDups            200%        70M           1s        
# bcftools sort                              70%        10M           1s
# bcftools concat                           100%        20M           1s
# truvari anno remap                        200%        20G          10m
# UltralongSvimasmInsExtractDupsPrime       200%        70M           1s
# truvari bench                             100%       100M           1s
# UltralongBed2IntervalVcf                  200%        50M           1s
#
task Impl {
    input {
        File samples_csv
        
        String remote_indir_query
        String remote_indir_svimasm
        String remote_outdir
        
        Int convert_svimasm_ins_to_dup

        Int svimasm_ins_use_gaps
        Int svimasm_ins_use_gaps_slack_bp
        Float svimasm_ins_use_gaps_length_similarity
        
        Int svimasm_ins_use_remap
        Int svimasm_ins_remap_max_length
        Float svimasm_ins_remap_cov_threshold

        Int match_to_gaps
        Int match_insdups_to_dups

        File reference_fa
        File reference_fai

        Int truvari_refdist
        Float truvari_pctsize
        Float truvari_pctovl
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 32
        Int disk_size_gb = 20
        Int preemptible_number = 3
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
        
        INFINITY="1000000000"
        samtools --version 1>&2
        bcftools --version 1>&2
        truvari --help 1>&2
        df -h 1>&2

        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            DIPCALL_BED=$(echo ${LINE} | cut -d , -f 2)
            
            # Skipping the sample if it has already been processed
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "1" )
            if [ ${TEST} != "1" ]; then
                continue
            fi


            # Extracting gaps from the dipcall BED
            # Remark: we do not handle gaps in the reference explicitly,
            # since we assume that calls in reference gaps have already been 
            # removed from the query VCFs upstream.
            gcloud storage cp ${DIPCALL_BED} ./${SAMPLE_ID}_dipcall.bed
            ${TIME_COMMAND} bedtools sort -i ${SAMPLE_ID}_dipcall.bed -g ~{reference_fai} > ${SAMPLE_ID}_dipcall_sorted.bed
            ${TIME_COMMAND} bedtools complement -L -i ${SAMPLE_ID}_dipcall_sorted.bed -g ~{reference_fai} > ${SAMPLE_ID}_gaps.bed
            rm -f ${SAMPLE_ID}_dipcall.bed ${SAMPLE_ID}_dipcall_sorted.bed
            N_GAPS=$(wc -l < ${SAMPLE_ID}_gaps.bed)
            echo "The dipcall BED of ${SAMPLE_ID} has ${N_GAPS} gaps" 1>&2


            # 1. Downloading and filtering the truth VCF -----------------------
            gcloud storage cp ~{remote_indir_svimasm}/${SAMPLE_ID}_canonized.vcf.gz ./${SAMPLE_ID}_svimasm.vcf.gz
            gcloud storage cp ~{remote_indir_svimasm}/${SAMPLE_ID}_canonized.vcf.gz.tbi ./${SAMPLE_ID}_svimasm.vcf.gz.tbi

            # 1.1 DEL
            bcftools filter --include 'SVTYPE=="DEL"' --output-type z ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_svimasm_del.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_del.vcf.gz

            # 1.2 DUP
            # We include the intervals of both DUP:TANDEM and of DUP:INT, even
            # though the latter (interspersed duplication) behaves like DUP just
            # in terms of depth, not in terms of breakpoints.
            bcftools filter --include 'SVTYPE=="DUP"' --output-type v ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_out.vcf
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
            # >=10kb alignments in either haplotype. A simple DUP >10kb might 
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
            # Thus, we use `truvari anno remap` on the remaining INS.
            # Working only on the remaining INS hopefully reduces the
            # runtime of this slow step. However, using two methods might 
            # give rise to two populations of DUPs with different breakpoint
            # decisions. Since these two populations are likely simple and
            # complex DUPs, which are likely to have distinct BAM patterns 
            # anyway, we accept this issue in exchange for the advantage in 
            # runtime.
            #
            # Remark: one could think of creating an INS for every DUP:INT 
            # record from svim-asm. However, we ran svim-asm with default 
            # params, which exclude `--interspersed_duplications_as_insertions`,
            # which means that an interspersed DUP is only represented as its 
            # source interval and the destination cannot be reconstructed. We 
            # should have run svim-asm twice, with and without the flag.
            bcftools filter --include 'SVTYPE=="INS"' --output-type v ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_svimasm_ins.vcf
            N_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_svimasm_ins.vcf | wc -l)
            if [ ${N_INS} -gt 0 ]; then
                # 1.3.1 Filtering all INS with the dipcall BED
                # Approx. 5% of all INS get marked as DUP in this step.
                if [ ~{svimasm_ins_use_gaps} -eq 1 -a ${N_GAPS} -gt 0 ]; then
                    ${TIME_COMMAND} java -cp ~{docker_dir} UltralongSvimasmInsExtractDups ${SAMPLE_ID}_svimasm_ins.vcf ${SAMPLE_ID}_gaps.bed $(wc -l < ${SAMPLE_ID}_gaps.bed) ~{svimasm_ins_use_gaps_slack_bp} ~{svimasm_ins_use_gaps_length_similarity} ${SAMPLE_ID}_svimasm_ins_dup.vcf ${SAMPLE_ID}_svimasm_ins_ins.vcf
                    rm -f ${SAMPLE_ID}_svimasm_ins.vcf ; mv ${SAMPLE_ID}_svimasm_ins_ins.vcf ${SAMPLE_ID}_svimasm_ins.vcf
                    N_INS_DUP=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_svimasm_ins_dup.vcf | wc -l)
                    if [ ${N_INS_DUP} -gt 0 ]; then
                        ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_svimasm_ins_dup.vcf --output ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                        rm -f ${SAMPLE_ID}_svimasm_ins_dup.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                    fi
                fi

                # 1.3.2 Filtering the remaining INS with truvari anno remap
                # Approx. 55% of the remaining INS get marked as DUP.
                N_INS=$(bcftools query --format '%ID\n' ${SAMPLE_ID}_svimasm_ins.vcf | wc -l)
                if [ ~{svimasm_ins_use_remap} -eq 1 -a ${N_INS} -gt 0 ]; then
                    bgzip --compress-level 1 ${SAMPLE_ID}_svimasm_ins.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins.vcf.gz
                    ${TIME_COMMAND} truvari anno remap --threads ${N_THREADS} --aligner minimap2 --min-length 1 --max-length ~{svimasm_ins_remap_max_length} --cov-threshold ~{svimasm_ins_remap_cov_threshold} -r ~{reference_fa} ${SAMPLE_ID}_svimasm_ins.vcf.gz -o ${SAMPLE_ID}_svimasm_ins_remap.vcf.gz
                    rm -f ${SAMPLE_ID}_svimasm_ins.vcf.gz*
                    ${TIME_COMMAND} java -cp ~{docker_dir} UltralongSvimasmInsExtractDupsPrime ${SAMPLE_ID}_svimasm_ins_remap.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf ${SAMPLE_ID}_svimasm_ins_ins.vcf
                    rm -f ${SAMPLE_ID}_svimasm_ins_remap.vcf.gz* ; mv ${SAMPLE_ID}_svimasm_ins_ins.vcf ${SAMPLE_ID}_svimasm_ins.vcf
                    ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf --output ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf.gz
                    rm -f ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf.gz
                    if [ -e ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz ]; then
                        ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type z ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
                        rm -f ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz* ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                    else
                        mv ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf.gz ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                        mv ${SAMPLE_ID}_svimasm_ins_dup_prime.vcf.gz.tbi ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz.tbi
                    fi
                fi
            fi
            bgzip --compress-level 1 ${SAMPLE_ID}_svimasm_ins.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_ins.vcf.gz
            if [ ! -e ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz ]; then
                bcftools view --header-only --output-type z ${SAMPLE_ID}_svimasm_ins.vcf.gz --output ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
                bcftools index -f -t ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz
            fi

            # 1.4 INV
            bcftools filter --include 'SVTYPE=="INV"' --output-type z ${SAMPLE_ID}_svimasm.vcf.gz --output ${SAMPLE_ID}_svimasm_inv.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_svimasm_inv.vcf.gz


            # 2. Computing matches of the query VCF ----------------------------
            # Remark: for speed reasons, sequence similarity is never used.
            gcloud storage cp ~{remote_indir_query}/${SAMPLE_ID}_'*.vcf.gz*' .

            # 3.1 DEL
            # Remark: svim-asm does not emit calls at some query DELs that 
            # correspond to gaps in the dipcall BED. DELs >10kb are likely
            # to create gaps in dipcall's BED by the definition of confident
            # BED, since such DELs induce split alignments (see above). We 
            # choose to mark as true every query DEL that matches a gap in
            # dipcall's BED.
            # Approx. 8% of all DEL get marked as true by a dipcall gap.
            ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_svimasm_del.vcf.gz -c ${SAMPLE_ID}_del.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
            mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_del1.vcf.gz
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_del1.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_del1.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_del1.vcf.gz
            rm -rf ${SAMPLE_ID}_truvari/
            if [ ~{match_to_gaps} -eq 1 -a ${N_GAPS} -gt 0 ]; then
                bcftools view --header-only ${SAMPLE_ID}_del.vcf.gz > ${SAMPLE_ID}_gaps.vcf
                ${TIME_COMMAND} java -cp ~{docker_dir} UltralongBed2IntervalVcf ${SAMPLE_ID}_gaps.bed DEL >> ${SAMPLE_ID}_gaps.vcf
                bgzip -@ ${N_THREADS} ${SAMPLE_ID}_gaps.vcf
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_gaps.vcf.gz
                ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_gaps.vcf.gz -c ${SAMPLE_ID}_del.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
                mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_del2.vcf.gz
                ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_del2.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
                mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_del2.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_del2.vcf.gz
                rm -rf ${SAMPLE_ID}_truvari/
                N_DEL1=$(bcftools index --nrecords ${SAMPLE_ID}_del1.vcf.gz)
                N_DEL2=$(bcftools index --nrecords ${SAMPLE_ID}_del2.vcf.gz)
                if [ ${N_DEL1} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_del1.vcf.gz*
                    mv ${SAMPLE_ID}_del2.vcf.gz ${SAMPLE_ID}_del_training.vcf.gz
                    mv ${SAMPLE_ID}_del2.vcf.gz.tbi ${SAMPLE_ID}_del_training.vcf.gz.tbi
                elif [ ${N_DEL2} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_del2.vcf.gz*
                    mv ${SAMPLE_ID}_del1.vcf.gz ${SAMPLE_ID}_del_training.vcf.gz
                    mv ${SAMPLE_ID}_del1.vcf.gz.tbi ${SAMPLE_ID}_del_training.vcf.gz.tbi
                else
                    ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type z ${SAMPLE_ID}_del1.vcf.gz ${SAMPLE_ID}_del2.vcf.gz --output ${SAMPLE_ID}_del_training.vcf.gz
                    bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_del_training.vcf.gz
                    rm -f ${SAMPLE_ID}_del1.vcf.gz* ${SAMPLE_ID}_del2.vcf.gz*
                fi
                rm -f ${SAMPLE_ID}_gaps.vcf.gz*
            else
                mv ${SAMPLE_ID}_del1.vcf.gz ${SAMPLE_ID}_del_training.vcf.gz
                mv ${SAMPLE_ID}_del1.vcf.gz.tbi ${SAMPLE_ID}_del_training.vcf.gz.tbi
            fi

            # 3.2 DUP
            # ------------> Fix comment once finalized!!!!!!!!!!!!!!!!!!
            # We compare the query DUPs to the entire DUP truth, which consists
            # of svim-asm's DUPs and svim-asm's INS->DUP.
            # Remark: DUP records from svim-asm seem to be generally 
            # comprehensive, and they do not need to be augmented with e.g. gaps
            # in dipcall's BED.
            ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_svimasm_dup.vcf.gz -c ${SAMPLE_ID}_dup.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
            mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_dup_training.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_dup_training.vcf.gz
            rm -rf ${SAMPLE_ID}_truvari/

            # 3.3 INSDUP
            # We compare the query INS->DUPs to both svim-asm's DUP and svim-
            # asm's INS->DUP.
            ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_svimasm_ins_dup.vcf.gz -c ${SAMPLE_ID}_insdup.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
            mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_insdup_training.vcf.gz
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_insdup_training.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_insdup_training.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_insdup_training.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_insdup_training.vcf.gz
            rm -rf ${SAMPLE_ID}_truvari/
            if [ ~{match_insdups_to_dups} -eq 1 ]; then
                ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_svimasm_dup.vcf.gz -c ${SAMPLE_ID}_insdup.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
                mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_insdup_training_prime.vcf.gz
                ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_insdup_training_prime.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_insdup_training_prime.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_insdup_training_prime.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_insdup_training_prime.vcf.gz
                rm -rf ${SAMPLE_ID}_truvari/
                N_INSDUP_TRAINING=$(bcftools index --nrecords ${SAMPLE_ID}_insdup_training.vcf.gz)
                N_INSDUP_TRAINING_PRIME=$(bcftools index --nrecords ${SAMPLE_ID}_insdup_training_prime.vcf.gz)
                if [ ${N_INSDUP_TRAINING} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_insdup_training.vcf.gz*
                    mv ${SAMPLE_ID}_insdup_training_prime.vcf.gz ${SAMPLE_ID}_insdup_training.vcf.gz
                    mv ${SAMPLE_ID}_insdup_training_prime.vcf.gz.tbi ${SAMPLE_ID}_insdup_training.vcf.gz.tbi
                elif [ ${N_INSDUP_TRAINING_PRIME} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_insdup_training_prime.vcf.gz*
                else
                    ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type z ${SAMPLE_ID}_insdup_training.vcf.gz ${SAMPLE_ID}_insdup_training_prime.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
                    rm -f ${SAMPLE_ID}_insdup_training.vcf.gz* ${SAMPLE_ID}_insdup_training_prime.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_insdup_training.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_insdup_training.vcf.gz
                fi
            fi
            
            # 3.4 INS
            # Remark: we don't use `--dup-to-ins` in truvari, since we 
            # assume that DUP and INS calls are already accurately
            # classified in both the query and the svim-asm VCF.
            ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_svimasm_ins.vcf.gz -c ${SAMPLE_ID}_ins.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
            mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_ins_training.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins_training.vcf.gz
            rm -rf ${SAMPLE_ID}_truvari/

            # 3.5 INV
            # Remark: svim-asm does not emit calls at some query INVs that 
            # correspond to gaps in the dipcall BED. We choose to mark as 
            # true such query INVs, at the risk of polluting the training 
            # set with events that are not simple INVs.
            # Approx. 15% of all INV get marked as true by a dipcall gap.
            ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_svimasm_inv.vcf.gz -c ${SAMPLE_ID}_inv.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
            mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_inv1.vcf.gz
            ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_inv1.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_inv1.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_inv1.vcf.gz
            rm -rf ${SAMPLE_ID}_truvari/
            if [ ~{match_to_gaps} -eq 1 -a ${N_GAPS} -gt 0 ]; then
                bcftools view --header-only ${SAMPLE_ID}_inv.vcf.gz > ${SAMPLE_ID}_gaps.vcf
                ${TIME_COMMAND} java -cp ~{docker_dir} UltralongBed2IntervalVcf ${SAMPLE_ID}_gaps.bed INV >> ${SAMPLE_ID}_gaps.vcf
                bgzip -@ ${N_THREADS} ${SAMPLE_ID}_gaps.vcf
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_gaps.vcf.gz
                ${TIME_COMMAND} truvari bench -b ${SAMPLE_ID}_gaps.vcf.gz -c ${SAMPLE_ID}_inv.vcf.gz --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --refdist ~{truvari_refdist} --pctseq 0 --pctsize ~{truvari_pctsize} --pctovl ~{truvari_pctovl} --pick single -o ./${SAMPLE_ID}_truvari/
                mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_inv2.vcf.gz
                ${TIME_COMMAND} bcftools sort --output-type z ${SAMPLE_ID}_inv2.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
                mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_inv2.vcf.gz
                bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_inv2.vcf.gz
                rm -rf ${SAMPLE_ID}_truvari/
                N_INV1=$(bcftools index --nrecords ${SAMPLE_ID}_inv1.vcf.gz)
                N_INV2=$(bcftools index --nrecords ${SAMPLE_ID}_inv2.vcf.gz)
                if [ ${N_INV1} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_inv1.vcf.gz*
                    mv ${SAMPLE_ID}_inv2.vcf.gz ${SAMPLE_ID}_inv_training.vcf.gz
                    mv ${SAMPLE_ID}_inv2.vcf.gz.tbi ${SAMPLE_ID}_inv_training.vcf.gz.tbi
                elif [ ${N_INV2} -eq 0 ]; then
                    rm -f ${SAMPLE_ID}_inv2.vcf.gz*
                    mv ${SAMPLE_ID}_inv1.vcf.gz ${SAMPLE_ID}_inv_training.vcf.gz
                    mv ${SAMPLE_ID}_inv1.vcf.gz.tbi ${SAMPLE_ID}_inv_training.vcf.gz.tbi
                else
                    ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type z ${SAMPLE_ID}_inv1.vcf.gz ${SAMPLE_ID}_inv2.vcf.gz --output ${SAMPLE_ID}_inv_training.vcf.gz
                    bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_inv_training.vcf.gz
                    rm -f ${SAMPLE_ID}_inv1.vcf.gz* ${SAMPLE_ID}_inv2.vcf.gz*
                fi
                rm -f ${SAMPLE_ID}_gaps.vcf.gz*
            else
                mv ${SAMPLE_ID}_inv1.vcf.gz ${SAMPLE_ID}_inv_training.vcf.gz
                mv ${SAMPLE_ID}_inv1.vcf.gz.tbi ${SAMPLE_ID}_inv_training.vcf.gz.tbi
            fi
            
            # Uploading and deallocating the sample
            gcloud storage mv ${SAMPLE_ID}_'*_training.vcf.gz*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            rm -rf ${SAMPLE_ID}_*
            ls -laht 1>&2
        done 3< ~{samples_csv}        
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
