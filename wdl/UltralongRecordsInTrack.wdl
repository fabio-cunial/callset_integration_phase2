version 1.0

# Adds basic TR and segdup annotations to every record in the input VCF.
#
workflow UltralongRecordsInTrack {
    input {
        File vcf_gz
        File vcf_tbi
        String svtype

        String output_prefix
        String remote_outdir

        File tr_bed
        File segdup_bed

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong"
    }

    call Impl {
        input:
            vcf_gz = vcf_gz,
            vcf_tbi = vcf_tbi,
            svtype = svtype,
            
            output_prefix = output_prefix,
            remote_outdir = remote_outdir,

            tr_bed = tr_bed,
            segdup_bed = segdup_bed,

            docker_image = docker_image
    }
}


#
task Impl {
    input {
        File vcf_gz
        File vcf_tbi
        String svtype

        String output_prefix
        String remote_outdir

        File tr_bed
        File segdup_bed

        String docker_image
        Int n_cpu = 4
        Int mem_gb = 8
        Int disk_size_gb = 50
    }

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"




        # ----------------------- Steps of the pipeline ------------------------

        function AnnotateStart() {
            local INPUT_BCF=$1
            local OUTPUT_BCF=$2
            local START_BED=$3
            local TRACK_BED=$4
            local TRACK_ID=$5

            ${TIME_COMMAND} bedtools intersect -wa -u -a ${START_BED} -b ${TRACK_BED} > start_track.bed
            awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t1\n",$1,$2,$4); }' start_track.bed > start_track.tsv
            ${TIME_COMMAND} bedtools intersect -v -wa -u -a ${START_BED} -b ${TRACK_BED} > start_track.bed
            awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t0\n",$1,$2,$4); }' start_track.bed >> start_track.tsv
            sort -k 1,1 -k 2,2n start_track.tsv | bgzip > start_track.tsv.gz
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 start_track.tsv.gz
            echo '##INFO=<ID=START_'${TRACK_ID}',Number=1,Type=Integer,Description="Left breakpoint is contained in a '${TRACK_ID}'">' > header.txt
            COLUMNS='CHROM,POS,~ID,INFO/START_'${TRACK_ID}
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations start_track.tsv.gz --header-lines header.txt --columns ${COLUMNS} --output-type b ${INPUT_BCF} --output ${OUTPUT_BCF}
            rm -f start_track* header.txt
        }


        function AnnotateEnd() {
            local INPUT_BCF=$1
            local OUTPUT_BCF=$2
            local END_BED=$3
            local TRACK_BED=$4
            local TRACK_ID=$5

            ${TIME_COMMAND} bedtools intersect -wa -u -a ${END_BED} -b ${TRACK_BED} > end_track.bed
            awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t1\n",$1,$2,$4); }' end_track.bed > end_track.tsv
            ${TIME_COMMAND} bedtools intersect -v -wa -u -a ${END_BED} -b ${TRACK_BED} > end_track.bed
            awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t0\n",$1,$2,$4); }' end_track.bed >> end_track.tsv
            sort -k 1,1 -k 2,2n end_track.tsv | bgzip > end_track.tsv.gz
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 end_track.tsv.gz
            echo '##INFO=<ID=END_'${TRACK_ID}',Number=1,Type=Integer,Description="Right breakpoint is contained in a '${TRACK_ID}'">' > header.txt
            COLUMNS='CHROM,POS,~ID,INFO/END_'${TRACK_ID}
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations end_track.tsv.gz --header-lines header.txt --columns ${COLUMNS} --output-type b ${INPUT_BCF} --output ${OUTPUT_BCF}
            rm -f end_track* header.txt
        }


        function AnnotateInterval() {
            local INPUT_BCF=$1
            local OUTPUT_BCF=$2
            local INTERVAL_BED=$3
            local TRACK_BED=$4
            local TRACK_ID=$5

            ${TIME_COMMAND} bedtools intersect -wa -u -f 1.0 -a ${INTERVAL_BED} -b ${TRACK_BED} > interval_track.bed
            awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t1\n",$1,$2,$4); }' interval_track.bed > interval_track.tsv
            ${TIME_COMMAND} bedtools intersect -v -wa -u -f 1.0 -a ${INTERVAL_BED} -b ${TRACK_BED} > interval_track.bed
            awk 'BEGIN { FS="\t"; OFS="\t"; } { printf("%s\t%d\t%s\t0\n",$1,$2,$4); }' interval_track.bed >> interval_track.tsv
            sort -k 1,1 -k 2,2n interval_track.tsv | bgzip > interval_track.tsv.gz
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 interval_track.tsv.gz
            echo '##INFO=<ID=INTERVAL_'${TRACK_ID}',Number=1,Type=Integer,Description="Interval is contained in a '${TRACK_ID}'">' > header.txt
            COLUMNS='CHROM,POS,~ID,INFO/INTERVAL_'${TRACK_ID}
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations interval_track.tsv.gz --header-lines header.txt --columns ${COLUMNS} --output-type b ${INPUT_BCF} --output ${OUTPUT_BCF}
            rm -f interval_track* header.txt
        }

        


        # ---------------------------- Main program ----------------------------

        bcftools view --output-type b ~{vcf_gz} --output input.bcf

        # Creating start, end, and interval BED.
        ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%INFO/SVLEN\t%ID\n' input.bcf > matrix.tsv
        ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } {
            printf("%s\t%d\t%d\t%s\n",$1,$2,$2,$4); \
        }' matrix.tsv | sort -k1,1 -k2,2n > start.bed
        if [ ~{svtype} != "INS" ]; then
            ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } { 
                printf("%s\t%d\t%d\t%s\n",$1, $2 + $3 - 1, $2 + $3 - 1, $4); \
            }' matrix.tsv | sort -k1,1 -k2,2n > end.bed
            ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } { 
                printf("%s\t%d\t%d\t%s\n",$1,$2, $2 + $3 - 1, $4); \
            }' matrix.tsv | sort -k1,1 -k2,2n > interval.bed
        fi
        rm -f matrix.tsv

        # Annotating with the TR track
        AnnotateStart input.bcf output.bcf start.bed ~{tr_bed} TR
        rm -f input.bcf ; mv output.bcf input.bcf
        if [ ~{svtype} != "INS" ]; then
            AnnotateEnd input.bcf output.bcf end.bed ~{tr_bed} TR
            rm -f input.bcf ; mv output.bcf input.bcf
            AnnotateInterval input.bcf output.bcf interval.bed ~{tr_bed} TR
            rm -f input.bcf ; mv output.bcf input.bcf
        fi

        # Annotating with the segdup track
        AnnotateStart input.bcf output.bcf start.bed ~{segdup_bed} SD
        rm -f input.bcf ; mv output.bcf input.bcf
        if [ ~{svtype} != "INS" ]; then
            AnnotateEnd input.bcf output.bcf end.bed ~{segdup_bed} SD
            rm -f input.bcf ; mv output.bcf input.bcf
            AnnotateInterval input.bcf output.bcf interval.bed ~{segdup_bed} SD
            rm -f input.bcf ; mv output.bcf input.bcf
        fi
        
        # Outputting
        bcftools view --output-type z --threads ${N_THREADS} input.bcf --output ${output_prefix}.vcf.gz
        bcftools index --threads ${N_THREADS} -f -t ${output_prefix}.vcf.gz
        gcloud storage cp ${output_prefix}.vcf.'gz*' ~{remote_outdir}
    >>>
    
    output {
    }

    runtime {
        cpu: n_cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " +  disk_size_gb + " HDD"
        preemptible: 0
        docker: docker_image
    }
}