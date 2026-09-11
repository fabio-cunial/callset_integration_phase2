version 1.0


# For each long call, computes the max R^2 coefficient with the short calls 
# located within a max distance from the long's endpoints. See `Rsquare.java` 
# for details.
#
workflow Rsquare {
    input {
        String chromosome_id
        String samples_id

        File input_bcf
        File input_csi
        File samples_txt

        Int matrix_only = 0

        Int min_long_length = 50
        Int max_short_length = 49
        Int max_distance_bp = 1000000
        Float min_af = 0.001

        String remote_outdir

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        input_bcf: "Every record must belong to a single autosome and be biallelic."
        samples_txt: "One sample ID per line, in any order. Samples that do not occur in the VCF are discarded, and the remaining ones are reordered as in `input_bcf`."
        max_short_length: "Must be <min_long_length"
        max_distance_bp: "Max bp distance to compare a long to a short"
        remote_outdir: "Without final slash"
        matrix_only: "1=the program only computes and saves the sparse projection of the VCF."
    }


    call Rsquare {
        input:
            chromosome_id = chromosome_id,
            samples_id = samples_id,

            input_bcf = input_bcf,
            input_csi = input_csi,
            samples_txt = samples_txt,

            matrix_only = matrix_only,

            min_long_length = min_long_length,
            max_short_length = max_short_length,
            max_distance_bp = max_distance_bp,
            min_af = min_af,

            remote_outdir = remote_outdir,

            docker_image = docker_image
    }
}


# Performance on chr22 (1,908,914 records in VCF; 1,905,492 records in sparse 
# matrix; 36,787 long; 1,868,705 short; 1Mbp max distance):
#
# TOOL                                       CPU        RAM         TIME
# 
# bcftools query | awk | gzip -1                                      1h
# Rsquare.java                              100%        4.3G       1h10m
#
task Rsquare {
    input {
        String chromosome_id
        String samples_id

        File input_bcf
        File input_csi
        File samples_txt

        Int matrix_only

        Int min_long_length
        Int max_short_length
        Int max_distance_bp
        Float min_af

        String remote_outdir

        String docker_image
        Int n_cpu = 3
        Int mem_gb = 8
        Int preemptible_number = 0
    }
    parameter_meta {
        n_cpu: "3 because there can be 3 concurrent processes in a pipe."
    }

    Int disk_size_gb = 10 + 10*( ceil(size(input_bcf,"GB")) )
    String docker_dir = "/callset_integration"

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_MB=$(( ~{mem_gb}*1024 - 512 ))
        TIME_COMMAND="/usr/bin/time --verbose"

        MATRIX_FILENAME="~{chromosome_id}_~{samples_id}.tsv.gz"
        RUN_ID="~{chromosome_id}_~{samples_id}_~{min_long_length}_~{max_short_length}_~{max_distance_bp}_~{min_af}"

        # Making sure that `samples_txt` is a subset of the samples in the VCF
        # and in the same relative order.
        bcftools query --list-samples ~{input_bcf} > vcf_samples.txt
        N_SAMPLES_IN_VCF=$(wc -l < vcf_samples.txt)
        awk '
        NR==FNR { requested[$0]=1; next }
        ($0 in requested)
        ' ~{samples_txt} vcf_samples.txt > samples.txt
        rm -f vcf_samples.txt
        N_SAMPLES_REQUESTED=$(wc -l < ~{samples_txt})
        N_SAMPLES_KEPT=$(wc -l < samples.txt)
        if [ ${N_SAMPLES_KEPT} -eq 0 ]; then
            echo "No requested samples found in the VCF." 1>&2
            exit 1
        fi
        echo "Number of samples in the VCF: ${N_SAMPLES_IN_VCF}" > ${RUN_ID}.log
        echo "Number of samples requested: ${N_SAMPLES_REQUESTED}" >> ${RUN_ID}.log
        echo "Number of samples kept for R^2 computation: ${N_SAMPLES_KEPT}" >> ${RUN_ID}.log

        # Building and caching a sparse projection of the VCF, that contains
        # only the requested samples and the information needed for R^2.
        #
        # Remark: the sparse matrix could be made more compressible by storing
        # all the sample indexes first, and then all their GT counts. We leave
        # this to the future.
        TEST=$(gcloud storage ls ~{remote_outdir}/matrices/${MATRIX_FILENAME} || echo "0")
        if [ "$TEST" != "0" ]; then
            gcloud storage cp ~{remote_outdir}/matrices/${MATRIX_FILENAME} .
        else
            date 1>&2
            bcftools query --format '%POS\t%REF\t%ALT\t%ID[\t%SAMPLE=%GT]\n' --include 'GT="alt" | GT="mis"' ~{input_bcf} | awk -F'\t' -v OFS='\t' '
            BEGIN {
                n_samples=0
                while ((getline sample_id < "samples.txt") > 0) { n_samples++; sample_index[sample_id]=n_samples }
            }
            {
                n_printed = 0
                for (i=4; i<=NF; i++) {
                    p = index($i, "="); sid = substr($i, 1, p-1)
                    if (!(sid in sample_index)) continue
                    n_alleles = split(substr($i, p+1), alleles, "[/|]")
                    alt_count = 0
                    for (j = 1; j <= n_alleles; j++) if (alleles[j] + 0 > 0) alt_count++
                    if (alt_count > 0) {
                        if (n_printed == 0) printf "%s\t%d\t%d\t%s", $1, length($2), length($3), $4
                        printf "\t%d=%d", sample_index[sid], alt_count
                        n_printed++
                    }
                }
                if (n_printed > 0) printf "\n"
            }' | gzip -4 -c > ${MATRIX_FILENAME}
            date 1>&2
            gcloud storage cp ${MATRIX_FILENAME} ${RUN_ID}.log ~{remote_outdir}/matrices/
            zcat ${MATRIX_FILENAME} | head -n 10 1>&2 || true
            ls -laht 1>&2
        fi

        # Computing R^2
        if [ ~{matrix_only} -eq 0 ]; then
            N_RECORDS_IN_VCF=$(bcftools index --nrecords ~{input_csi})
            echo "Number of records in the VCF: ${N_RECORDS_IN_VCF}" >> ${RUN_ID}.log
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M Rsquare ${MATRIX_FILENAME} ~{chromosome_id} ${N_SAMPLES_KEPT} ~{min_long_length} ~{max_short_length} ~{max_distance_bp} ~{min_af} ${RUN_ID}.bed >> ${RUN_ID}.log
            gcloud storage mv ${RUN_ID}.bed ${RUN_ID}.log ~{remote_outdir}/
        fi
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