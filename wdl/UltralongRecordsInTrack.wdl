version 1.0

#
workflow UltralongRecordsInTrack {
    input {
        File vcf_gz
        File vcf_tbi

        File tr_bed
        File segdup_bed
        File segdup_gt10kb_bed

        String docker_image
    }

    call Impl {
        input:
            vcf_gz = vcf_gz,
            vcf_tbi = vcf_tbi,

            tr_bed = tr_bed,
            segdup_bed = segdup_bed,
            segdup_gt10kb_bed = segdup_gt10kb_bed, 

            docker_image = docker_image
    }
}


#
task Impl {
    input {
        File vcf_gz
        File vcf_tbi

        File tr_bed
        File segdup_bed
        File segdup_gt10kb_bed

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


        ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%INFO/SVLEN\n' ~{vcf_gz} > matrix.tsv
        ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } { 
            $3 = $2; \
            print $0 \
        }' matrix.tsv | sort -k1,1 -k2,2n > start.bed
        ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } { 
            $2 = $2 + $3 - 1; \
            $3 = $2; \
            print $0 \
        }' matrix.tsv | sort -k1,1 -k2,2n > end.bed
        ${TIME_COMMAND} awk 'BEGIN { FS="\t"; OFS="\t"; } { 
            $3 = $2 + $3 - 1; \
            print $0 \
        }' matrix.tsv | sort -k1,1 -k2,2n > interval.bed
        rm -f matrix.tsv
        
        ${TIME_COMMAND} bedtools intersect -wa -u -a start.bed -b ~{tr_bed} > start_tr.bed &
        ${TIME_COMMAND} bedtools intersect -wa -u -a end.bed -b ~{tr_bed} > end_tr.bed &
        ${TIME_COMMAND} bedtools intersect -wa -u -f 1.0 -a interval.bed -b ~{tr_bed} > interval_tr.bed &
        wait
        ${TIME_COMMAND} bedtools intersect -wa -u -a start.bed -b ~{segdup_bed} > start_segdup.bed &
        ${TIME_COMMAND} bedtools intersect -wa -u -a end.bed -b ~{segdup_bed} > end_segdup.bed &
        ${TIME_COMMAND} bedtools intersect -wa -u -f 1.0 -a interval.bed -b ~{segdup_bed} > interval_segdup.bed &
        wait
        ${TIME_COMMAND} bedtools intersect -wa -u -a start.bed -b ~{segdup_gt10kb_bed} > start_segdup_gt10kb.bed &
        ${TIME_COMMAND} bedtools intersect -wa -u -a end.bed -b ~{segdup_gt10kb_bed} > end_segdup_gt10kb.bed &
        ${TIME_COMMAND} bedtools intersect -wa -u -f 1.0 -a interval.bed -b ~{segdup_gt10kb_bed} > interval_segdup_gt10kb.bed &
        wait
        ls -laht 1>&2
    >>>
    
    output {
        File start_tr_bed = "start_tr.bed"
        File end_tr_bed = "end_tr.bed"
        File interval_tr_bed = "interval_tr.bed"

        File start_segdup_bed = "start_segdup.bed"
        File end_segdup_bed = "end_segdup.bed"
        File interval_segdup_bed = "interval_segdup.bed"

        File start_segdup_gt10kb_bed = "start_segdup_gt10kb.bed"
        File end_segdup_gt10kb_bed = "end_segdup_gt10kb.bed"
        File interval_segdup_gt10kb_bed = "interval_segdup_gt10kb.bed"
    }

    runtime {
        cpu: n_cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " +  disk_size_gb + " HDD"
        preemptible: 0
        docker: docker_image
    }
}