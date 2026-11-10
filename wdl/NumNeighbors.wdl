version 1.0


# For each long call, computes the number of long neighbors located within a max
# distance.
#
workflow NumNeighbors {
    input {
        String chromosome_id

        File input_bcf
        File input_csi

        File tr_bed
        Int tr_bed_slack_bp = 0
        File reference_fai

        Int min_long_length = 50
        Int max_distance_bp = 1000
        Float min_overlap_fraction = 0.9

        String remote_outdir

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        input_bcf: "Assumed to contain all chromosomes"
        min_long_length: "Only records with length >=this are considered"
        max_distance_bp: "Max bp distance to compare two records"
        min_overlap_fraction: "Records whose interval is covered by the TR track by >=this fraction are considered as being inside TRs."
        tr_bed_slack_bp: "Number of bps to add as slack to the TR BED. This is useful for INS."
        remote_outdir: "Without final slash"
    }


    call Impl {
        input:
            chromosome_id = chromosome_id,

            input_bcf = input_bcf,
            input_csi = input_csi,

            tr_bed = tr_bed,
            reference_fai = reference_fai,

            min_long_length = min_long_length,
            max_distance_bp = max_distance_bp,
            tr_bed_slack_bp = tr_bed_slack_bp,
            min_overlap_fraction = min_overlap_fraction,

            remote_outdir = remote_outdir,

            docker_image = docker_image
    }
}


# Performance on chr22:
#
# TOOL                                            CPU        RAM         TIME
# 
# bcftools view                                  100%       300M          26s
# truvari anno numneigh | bedtools coverage       50%       200M          13s
# awk | bcftools query                                                     5s
#
task Impl {
    input {
        String chromosome_id

        File input_bcf
        File input_csi

        File tr_bed
        Int tr_bed_slack_bp
        File reference_fai

        Int min_long_length
        Int max_distance_bp
        Float min_overlap_fraction

        String remote_outdir

        String docker_image
        Int n_cpu = 2
        Int mem_gb = 4
        Int disk_size_gb = 100
    }
    parameter_meta {
    }

    String docker_dir = "/callset_integration"

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_MB=$(( ~{mem_gb}*1024 - 512 ))
        TIME_COMMAND="/usr/bin/time --verbose"

        RUN_ID="~{chromosome_id}_~{min_long_length}_~{max_distance_bp}_~{min_overlap_fraction}"

        if [ ~{tr_bed_slack_bp} -gt 0 ]; then
            ${TIME_COMMAND} bedtools slop -b ~{tr_bed_slack_bp} -i ~{tr_bed} -g ~{reference_fai} > tr.bed
        else
            mv ~{tr_bed} tr.bed
        fi
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --drop-genotypes --include 'ABS(SVLEN)>='~{min_long_length} --output-type z ~{input_bcf} ~{chromosome_id} --output chromosome.vcf.gz
        ${TIME_COMMAND} bcftools index -t chromosome.vcf.gz
        date 1>&2
        ${TIME_COMMAND} truvari anno numneigh --sizemin 0 --refdist ~{max_distance_bp} chromosome.vcf.gz | bedtools coverage -header -a - -b tr.bed > coverage.tsv
        date 1>&2
        awk -F'\t' -v OFS='\t' '/^#/ { print; next } $NF >='~{min_overlap_fraction}' { NF -= 4; print }' coverage.tsv | bcftools query --format '%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/NumNeighbors\n' > ${RUN_ID}_tr.tsv
        date 1>&2
        head -n 10 ${RUN_ID}_tr.tsv 1>&2 || true
        date 1>&2
        awk -F'\t' -v OFS='\t' '/^#/ { print; next } $NF <'~{min_overlap_fraction}' { NF -= 4; print }' coverage.tsv | bcftools query --format '%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/NumNeighbors\n' > ${RUN_ID}_not_tr.tsv
        date 1>&2
        head -n 10 ${RUN_ID}_not_tr.tsv 1>&2 || true
        gcloud storage cp ${RUN_ID}'*.tsv' ~{remote_outdir}/
    >>>
    
    output {
    }

    runtime {
        cpu: n_cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " +  disk_size_gb + " SSD"
        preemptible: 0
        docker: docker_image
    }
}