version 1.0


# For each long call, computes the smallest distance from every other long
# call located within a max distance from its endpoints. See
# `SmallestDistance.java` for details. This is essentially identical to
# `Rsquare.wdl`.
#
workflow SmallestDistance {
    input {
        String chromosome_id

        File input_bcf
        File input_csi

        Int matrix_only = 0

        Int min_long_length = 50
        Int max_distance_bp = 1000000
        Float min_af = 0.001
        Int distance_type = 0

        String remote_outdir

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        input_bcf: "Assumed to contain all chromosomes"
        min_long_length: "Minimum length for a long call"
        max_distance_bp: "Max bp distance to compare two longs"
        remote_outdir: "Without final slash"
        matrix_only: "1=the program only computes and saves the sparse projection of the VCF."
        distance_type: "0=first, 1=last, 2=gap."
    }


    call SmallestDistance {
        input:
            chromosome_id = chromosome_id,

            input_bcf = input_bcf,
            input_csi = input_csi,

            matrix_only = matrix_only,

            min_long_length = min_long_length,
            max_distance_bp = max_distance_bp,
            min_af = min_af,
            distance_type = distance_type,

            remote_outdir = remote_outdir,

            docker_image = docker_image
    }
}


# Performance on chr22 (1,908,914 records in VCF; 1,905,492 records in sparse 
# matrix; 36,787 long; 1,868,705 short; 1Mbp max distance):
#
# TOOL                                       CPU        RAM         TIME
# 
# bcftools query | awk | gzip -1                                      
# SmallestDistance.java                     
#
task SmallestDistance {
    input {
        String chromosome_id

        File input_bcf
        File input_csi

        Int matrix_only

        Int min_long_length
        Int max_distance_bp
        Float min_af
        Int distance_type

        String remote_outdir

        String docker_image
        Int n_cpu = 4
        Int mem_gb = 8
        Int disk_size_gb = 100
    }
    parameter_meta {
        n_cpu: "There can be 4 concurrent processes in a pipe."
    }

    String docker_dir = "/callset_integration"

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_MB=$(( ~{mem_gb}*1024 - 512 ))
        TIME_COMMAND="/usr/bin/time --verbose"

        MATRIX_FILENAME="~{chromosome_id}.tsv.gz"
        RUN_ID="~{chromosome_id}_~{min_long_length}_~{min_af}_~{max_distance_bp}_~{distance_type}"

        # Building and caching a sparse projection of the VCF, that contains
        # only the information needed for smallest distance computation.
        TEST=$(gcloud storage ls ~{remote_outdir}/matrices/${MATRIX_FILENAME} || echo "0")
        if [ "$TEST" != "0" ]; then
            gcloud storage cp ~{remote_outdir}/matrices/${MATRIX_FILENAME} .
        else
            date 1>&2
            bcftools view --output-type u ~{input_bcf} ~{chromosome_id} | bcftools +fill-tags --output-type u - -- -t AF | bcftools query --format '%POS\t%INFO/SVLEN\t%INFO/AF\t%ID\t%INFO/SVTYPE\n' | gzip -4 -c > ${MATRIX_FILENAME}
            date 1>&2
            gcloud storage cp ${MATRIX_FILENAME} ~{remote_outdir}/matrices/
            zcat ${MATRIX_FILENAME} | head -n 10 1>&2 || true
            ls -laht 1>&2
        fi

        # Computing smallest distance
        if [ ~{matrix_only} -eq 0 ]; then
            N_RECORDS_IN_VCF=$(bcftools index --nrecords ~{input_csi})
            echo "Number of records in the VCF: ${N_RECORDS_IN_VCF}" >> ${RUN_ID}.log
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M SmallestDistance ${MATRIX_FILENAME} ~{chromosome_id} ~{min_long_length} ~{min_af} ~{max_distance_bp} ~{distance_type} ${RUN_ID}.bed >> ${RUN_ID}.log
            gcloud storage mv ${RUN_ID}.bed ${RUN_ID}.log ~{remote_outdir}/
        fi
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