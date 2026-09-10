version 1.0

#
workflow Rsquare {
    input {
        String chr_id
        File input_bcf
        File input_csi
        File input_bed
        File reference_fai
        String remote_outdir

        Int min_long_length = 50
        Int max_short_length = 49
        Int max_distance_bp = 1000000

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        input_bcf: "Every record must belong to a single chromosome and be biallelic."
        input_bed: "Whole-genome BED, not necessarily limited to the current chromosome."
        max_short_length: "Must be <min_long_length"
    }


    call SubsetByBed {
        input:
            chr_id = chr_id,
            input_bcf = input_bcf,
            input_csi = input_csi,
            input_bed = input_bed,
            reference_fai = reference_fai,
            docker_image = docker_image
    }
    call Rsquare as all {
        input:
            id = chr_id + "_all",
            input_bcf = input_bcf,
            input_csi = input_csi,
            remote_outdir = remote_outdir,

            min_long_length = min_long_length,
            max_short_length = max_short_length,
            max_distance_bp = max_distance_bp,

            docker_image = docker_image
    }
    call Rsquare as in_bed {
        input:
            id = chr_id + "_in_bed",
            input_bcf = SubsetByBed.in_bed_bcf,
            input_csi = SubsetByBed.in_bed_csi,
            remote_outdir = remote_outdir,

            min_long_length = min_long_length,
            max_short_length = max_short_length,
            max_distance_bp = max_distance_bp,

            docker_image = docker_image
    }
    call Rsquare as not_in_bed {
        input:
            id = chr_id + "_not_in_bed",
            input_bcf = SubsetByBed.not_in_bed_bcf,
            input_csi = SubsetByBed.not_in_bed_csi,
            remote_outdir = remote_outdir,

            min_long_length = min_long_length,
            max_short_length = max_short_length,
            max_distance_bp = max_distance_bp,

            docker_image = docker_image
    }
}


task SubsetByBed {
    input {
        String chr_id
        File input_bcf
        File input_csi
        File input_bed
        File reference_fai

        String docker_image
        Int n_cpu = 4
        Int mem_gb = 4
        Int preemptible_number = 0
    }

    Int disk_size_gb = 10 + 4*( ceil(size(input_bcf,"GB")) )
    String docker_dir = "/callset_integration"

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"

        ${TIME_COMMAND} bedtools sort -i ~{input_bed} -faidx ~{reference_fai} > sorted.bed
        ${TIME_COMMAND} bedtools complement -i sorted.bed -g ~{reference_fai} > complement.bed
        ${TIME_COMMAND} awk -F'\t' -v c=~{chr_id} '$1==c' sorted.bed > sorted_chr.bed
        ${TIME_COMMAND} awk -F'\t' -v c=~{chr_id} '$1==c' complement.bed > complement_chr.bed
        if [ -s sorted_chr.bed ]; then
            ${TIME_COMMAND} bcftools view --threads $(( ${N_THREADS} / 2 )) --output-type b --regions-file sorted_chr.bed     --regions-overlap pos --write-index ~{input_bcf} --output in_bed.bcf & PID1=$!
        else
            bcftools view --header-only --output-type b --write-index ~{input_bcf} --output in_bed.bcf & PID1=$!
        fi
        if [ -s complement_chr.bed ]; then
            ${TIME_COMMAND} bcftools view --threads $(( ${N_THREADS} / 2 )) --output-type b --targets-file complement_chr.bed --targets-overlap pos --write-index ~{input_bcf} --output not_in_bed.bcf & PID2=$!
        else
            bcftools view --header-only --output-type b --write-index ~{input_bcf} --output not_in_bed.bcf & PID2=$!
        fi
        wait ${PID1} ; wait ${PID2}
    >>>
    
    output {
        File in_bed_bcf = "in_bed.bcf"
        File in_bed_csi = "in_bed.bcf.csi"
        File not_in_bed_bcf = "not_in_bed.bcf"
        File not_in_bed_csi = "not_in_bed.bcf.csi"
    }
    runtime {
        cpu: n_cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " +  disk_size_gb + " SSD"
        preemptible: preemptible_number
        docker: docker_image
    }
}


# TOOL                                                CPU     RAM     TIME
# 
# bcftools query | awk | gzip
# Rsquare.java
#
task Rsquare {
    input {
        String id
        File input_bcf
        File input_csi
        String remote_outdir

        Int min_long_length
        Int max_short_length
        Int max_distance_bp

        String docker_image
        Int n_cpu = 4
        Int mem_gb = 16
        Int preemptible_number = 0
    }

    Int disk_size_gb = 10 + 10*( ceil(size(input_bcf,"GB")) )
    String docker_dir = "/callset_integration"

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_MB=$(( ~{mem_gb}*1024 - 512 ))
        TIME_COMMAND="/usr/bin/time --verbose"

        # Building a sparse projection of the VCF
        bcftools view --header-only ~{input_bcf} | tail -n 1 | tr '\t' '\n' | tail -n +10 > samples.txt
        N_SAMPLES=$(wc -l < samples.txt)
        date 1>&2
        bcftools query --format '%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' --include 'GT="alt" | GT="mis"' ~{input_bcf} | awk -F'\t' -v OFS='\t' '
        BEGIN {
            n_samples=0
            while ((getline sample_id < "samples.txt") > 0) { n_samples++; sample_index[sample_id]=n_samples }
        }
        {
            n_printed = 0
            for (i=4; i<=NF; i++) {
                p = index($i, "="); sid = substr($i, 1, p-1)
                n_alleles = split(substr($i, p+1), alleles, "[/|]")
                alt_count = 0
                for (j = 1; j <= n_alleles; j++) if (alleles[j] + 0 > 0) alt_count++
                if (alt_count > 0) {
                    if (n_printed == 0) printf "%s\t%d\t%d", $1, length($2), length($3)
                    printf "\t%d=%d", sample_index[sid], alt_count
                    n_printed++
                }
            }
            if (n_printed > 0) printf "\n"
        }' | gzip -1 -c > matrix.tsv.gz
        date 1>&2

        # Computing R^2
        bcftools index --nrecords ~{input_csi} 1>&2
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M Rsquare matrix.tsv.gz ${N_SAMPLES} ~{min_long_length} ~{max_short_length} ~{max_distance_bp} ~{id}.tsv
        gcloud storage mv ~{id}.tsv ~{remote_outdir}/
    >>>
    
    output {
    }

    runtime {
        cpu: n_cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " +  disk_size_gb + " SSD"
        preemptible: preemptible_number
        docker: docker_image
    }
}