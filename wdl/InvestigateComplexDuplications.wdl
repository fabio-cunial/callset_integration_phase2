version 1.0

#
workflow InvestigateComplexDuplications {
    input {
        File chunk_tsv
        String remote_outdir

        File reference_fa
        File reference_fai

        Float truvari_anno_remap_cov_threshold = 0.8
        String min_coverage_fractions = "0.9 0.8 0.7"

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong_remap:latest"
    }
    parameter_meta {
        chunk_tsv: "Format: sample_id, ultralong_bcf, hap1_fa, hap2_fa"
        min_coverage_fractions: "Space-separated"
    }

    call Impl {
        input:
            chunk_tsv = chunk_tsv,
            remote_outdir = remote_outdir,

            reference_fa = reference_fa,
            reference_fai = reference_fai,

            truvari_anno_remap_cov_threshold = truvari_anno_remap_cov_threshold,
            min_coverage_fractions = min_coverage_fractions,

            docker_image = docker_image
    }
}


# TOOL                                                CPU     RAM     TIME
# truvari anno remap                                 300%     20G       3m
# minimap2 haps                                      200%     13G       2m
# 
task Impl {
    input {
        File chunk_tsv
        String remote_outdir

        File reference_fa
        File reference_fai

        Float truvari_anno_remap_cov_threshold
        String min_coverage_fractions

        String docker_image
        Int n_cpu = 3
        Int mem_gb = 25
        Int disk_size_gb = 50
        Int preemptible_number = 0
    }

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        EFFECTIVE_RAM_MB=$(( ~{mem_gb}*1024 - 500 ))
        EFFECTIVE_RAM_MB_PER_THREAD=$(( ${EFFECTIVE_RAM_MB} / ${N_THREADS} ))
        



        # --------------------- Steps of the pipeline --------------------------

        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            local ULTRALONG_BCF=$(echo ${LINE} | cut -d , -f 2)
            local HAP1_FA=$(echo ${LINE} | cut -d , -f 3)
            local HAP2_FA=$(echo ${LINE} | cut -d , -f 4)

            gcloud storage cp ${ULTRALONG_BCF} ${SAMPLE_ID}.bcf
            gcloud storage cp ${HAP1_FA} ${SAMPLE_ID}_hap1.fa
            gcloud storage cp ${HAP2_FA} ${SAMPLE_ID}_hap2.fa
        }


        function DelocalizeSample() {
            local SAMPLE_ID=$1

            rm -f ${SAMPLE_ID}.bcf
            rm -f ${SAMPLE_ID}_hap1.fa
            rm -f ${SAMPLE_ID}_hap2.fa
        }


        # Stores in a multi-FASTA file every INS sequence that is annotated as 
        # tandem by `truvari anno remap`.
        #
        # Remark: we assume that INSDUPs have been converted back to INS.
        #
        function GetRemapSequences() {
            local SAMPLE_ID=$1
            local INPUT_BCF=$2
            local OUTPUT_FA=$3
            
            mv ${INPUT_BCF} ${SAMPLE_ID}_in.bcf
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'INFO/SVTYPE=="INS"' --output-type v ${SAMPLE_ID}_in.bcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.bcf; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            ${TIME_COMMAND} truvari anno remap --threads ${N_THREADS} --aligner minimap2 --min-length 1 --max-length 3000000000 --cov-threshold ~{truvari_anno_remap_cov_threshold} -r ./ref_files/${REF_FA_BASENAME} ${SAMPLE_ID}_in.vcf -o ${SAMPLE_ID}_out.vcf
            mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%SUPP_PBSV\t%SUPP_SNIFFLES\t%SUPP_PAV\t%remap_classification\t%remap_coords\t%remap_ori\t%remap_perc\t%remap_segments\t%remap_local_rescue\t%ALT\n' --include 'INFO/remap_classification=="tandem" || INFO/remap_classification=="tandem_complex" || INFO/remap_classification=="tandem_inverted"' ${SAMPLE_ID}_in.vcf | awk -F'\t' '{ printf ">%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n%s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13 }' > ${SAMPLE_ID}_out.fa
            rm -f ${SAMPLE_ID}_in.vcf

            mv ${SAMPLE_ID}_out.fa ${OUTPUT_FA}
        }


        # Returns the query sequences that are covered by >=1 alignment at
        # >=`MIN_COVERAGE_FRACTION`.
        #
        function FilterFullyAlignedSequences() {
            local SAMPLE_ID=$1
            local MIN_COVERAGE_FRACTION=$2
            local HAP_ID=$3

            awk -F'\t' -v min_coverage_fraction=${MIN_COVERAGE_FRACTION} '
                function qcov(cigar,    len, op, aligned, total) {
                    aligned = 0
                    total = 0
                    while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
                        len = substr(cigar, 1, RLENGTH - 1) + 0
                        op  = substr(cigar, RLENGTH, 1)
                        if (op ~ /[MI=X]/) {
                            aligned += len
                            total += len
                        } else if (op ~ /[SH]/) {
                            total += len
                        }
                        cigar = substr(cigar, RLENGTH + 1)
                    }
                    return (total > 0) ? (aligned / total) : 0
                }
                NR == FNR {
                    if ($0 !~ /^@/ && int($2 / 4) % 2 == 0 && $6 != "*" && qcov($6) >= min_coverage_fraction) {
                        keep[$1] = 1
                    }
                    next
                }
                END { for (q in keep) print q }
            ' ${SAMPLE_ID}_${HAP_ID}_remapped_alignments.sam > ${SAMPLE_ID}_${HAP_ID}_${MIN_COVERAGE_FRACTION}.txt
        }


        # Given a set of remapped sequences and a haplotype, returns the 
        # remapped sequences that are covered by >=1 alignment at >=
        # `MIN_COVERAGE_FRACTION`.
        #
        function GetFullyAlignedSequences() {
            local SAMPLE_ID=$1
            local HAP_FA=$2
            local HAP_ID=$3
            local REMAPPED_FA=$4

            ${TIME_COMMAND} minimap2 -t ${N_THREADS} -a -x map-hifi ${HAP_FA} ${REMAPPED_FA} > ${SAMPLE_ID}_${HAP_ID}_remapped_alignments.sam
            for MIN_COVERAGE_FRACTION in ~{min_coverage_fractions}; do
                FilterFullyAlignedSequences ${SAMPLE_ID} ${MIN_COVERAGE_FRACTION} ${HAP_ID}
            done
            ${TIME_COMMAND} samtools sort -@ ${N_THREADS} -m ${EFFECTIVE_RAM_MB_PER_THREAD}M --no-PG -o ${SAMPLE_ID}_${HAP_ID}_remapped_alignments.bam ${SAMPLE_ID}_${HAP_ID}_remapped_alignments.sam
            rm -f ${SAMPLE_ID}_${HAP_ID}_remapped_alignments.sam
        }




        # -------------------------- Main program ------------------------------

        # Preparing the reference for `truvari anno remap`.
        #
        # Remark: the reference should be minimap-indexed only once and reused
        # by every sample.
        mkdir ./ref_files
        mv ~{reference_fa} ./ref_files/
        mv ~{reference_fai} ./ref_files/
        REF_FA_BASENAME=$(basename ~{reference_fa})

        # Processing each sample in the chunk file
        cat ~{chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE || [ -n "${LINE}" ]; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            LocalizeSample ${SAMPLE_ID} ${LINE}
            GetRemapSequences ${SAMPLE_ID} ${SAMPLE_ID}.bcf ${SAMPLE_ID}_remapped.fa
            N_REMAPPED_SEQUENCES=$(wc -l < ${SAMPLE_ID}_remapped.fa)
            N_REMAPPED_SEQUENCES=$(( ${N_REMAPPED_SEQUENCES} / 2 ))
            echo "${SAMPLE_ID} has ${N_REMAPPED_SEQUENCES} remapped sequences" 2>&1
            if [ ${N_REMAPPED_SEQUENCES} -eq 0 ]; then
                DelocalizeSample ${SAMPLE_ID}
                continue
            fi
            GetFullyAlignedSequences ${SAMPLE_ID} ${SAMPLE_ID}_hap1.fa hap1 ${SAMPLE_ID}_remapped.fa
            GetFullyAlignedSequences ${SAMPLE_ID} ${SAMPLE_ID}_hap2.fa hap2 ${SAMPLE_ID}_remapped.fa
            for MIN_COVERAGE_FRACTION in ~{min_coverage_fractions}; do
                cat ${SAMPLE_ID}_hap1_${MIN_COVERAGE_FRACTION}.txt ${SAMPLE_ID}_hap2_${MIN_COVERAGE_FRACTION}.txt | sort | uniq > ${SAMPLE_ID}_${MIN_COVERAGE_FRACTION}_remapped_aligned.txt
                N_ALIGNED_SEQUENCES=$(wc -l < ${SAMPLE_ID}_${MIN_COVERAGE_FRACTION}_remapped_aligned.txt)
                PERCENT=$( bc <<< "scale=2; (${N_ALIGNED_SEQUENCES}*100)/${N_REMAPPED_SEQUENCES}" )
                echo "${SAMPLE_ID} has ${N_ALIGNED_SEQUENCES} remapped sequences covered >=${MIN_COVERAGE_FRACTION} by an alignment (${PERCENT}%)." 2>&1
            done
            rm -f ${SAMPLE_ID}_hap*_*.txt
            gcloud storage mv ${SAMPLE_ID}_remapped.fa ${SAMPLE_ID}'_*_'remapped_aligned.txt ${SAMPLE_ID}_hap1_remapped_alignments.bam ${SAMPLE_ID}_hap2_remapped_alignments.bam ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
        done 3< chunk.csv
    >>>
    
    output {
    }

    runtime {
        cpu: n_cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " +  disk_size_gb + " HDD"
        preemptible: preemptible_number
        docker: docker_image
    }
}
