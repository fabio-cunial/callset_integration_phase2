version 1.0


# 
#
workflow SameChromPos {
    input {
        String chunk_id
        File sv_integration_chunk_tsv
        String remote_indir
        Int main_or_ultralong

        File tr_bed
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        main_or_ultralong: "0=main, 1=ultralong"
    }
    
    call Impl {
        input:
            chunk_id = chunk_id,
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            remote_indir = remote_indir,
            main_or_ultralong = main_or_ultralong,

            tr_bed = tr_bed,

            docker_image = docker_image
    }
    
    output {
        File stats_lt_50 = Impl.stats_lt_50
        File stats_ge_50 = Impl.stats_ge_50
        File cluster_sizes_lt_50 = Impl.cluster_sizes_lt_50
        File cluster_sizes_ge_50 = Impl.cluster_sizes_ge_50
        File deltas_lt_50 = Impl.deltas_lt_50
        File deltas_ge_50 = Impl.deltas_ge_50
        File different_annotations_lt_50 = Impl.different_annotations_lt_50
        File different_annotations_ge_50 = Impl.different_annotations_ge_50
        File n_tr_clusters_lt_50 = Impl.n_tr_clusters_lt_50
        File n_tr_clusters_ge_50 = Impl.n_tr_clusters_ge_50
    }
}


task Impl {
    input {
        String chunk_id
        File sv_integration_chunk_tsv
        String remote_indir
        Int main_or_ultralong

        File tr_bed
        
        String docker_image
        Int n_cpu = 1
        Int ram_size_gb = 8
        Int disk_size_gb = 20
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

        function ProcessVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF=$2
            local OUTPUT_CLUSTER_SIZE_TXT=$3  # Write
            local OUTPUT_DELTAS_CSV=$4  # Write
            local OUTPUT_COUNTS_CSV=$5  # Append
            local OUTPUT_ANNOTATIONS_CSV=$6  # Append
            local OUTPUT_TR_CSV=$7  # Append

            # Distinct pairs
            N_DISTINCT_PAIRS=$(bcftools query --format '%CHROM\t%POS\n' ${INPUT_VCF} | sort -k1,1 -k2,2n | uniq | wc -l)

            # Detailed counts
            java -cp ~{docker_dir} SameChromPos ${INPUT_VCF} ${OUTPUT_CLUSTER_SIZE_TXT} ${OUTPUT_DELTAS_CSV} ${SAMPLE_ID}_sites.vcf 1>> ${OUTPUT_COUNTS_CSV} 2>> ${OUTPUT_ANNOTATIONS_CSV}

            # Tandem repeats
            bgzip ${SAMPLE_ID}_sites.vcf ; tabix -f ${SAMPLE_ID}_sites.vcf.gz
            N_PROBLEMATIC_PAIRS=$(bcftools index --nrecords ${SAMPLE_ID}_sites.vcf.gz)
            N_PROBLEMATIC_PAIRS_IN_TRS=$(bcftools view --no-header --regions-file ~{tr_bed} --regions-overlap pos ${SAMPLE_ID}_sites.vcf.gz | wc -l)
            echo "${N_DISTINCT_PAIRS},${N_PROBLEMATIC_PAIRS},${N_PROBLEMATIC_PAIRS_IN_TRS}" >> ${OUTPUT_TR_CSV}
            rm -f ${SAMPLE_ID}_sites.vcf.gz*
        }




        # --------------------------- Main program -----------------------------
        
        touch ~{chunk_id}_stats_lt_50.csv 
        touch ~{chunk_id}_stats_ge_50.csv
        touch ~{chunk_id}_different_annotations_lt_50.csv
        touch ~{chunk_id}_different_annotations_ge_50.csv
        touch ~{chunk_id}_trs_lt_50.csv
        touch ~{chunk_id}_trs_ge_50.csv
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            if [ ~{main_or_ultralong} -eq 0 ]; then
                SUFFIX="_kanpig.vcf.gz"
            else
                SUFFIX="_ultralong.bcf"
            fi
            gcloud storage cp ~{remote_indir}/${SAMPLE_ID}${SUFFIX}'*' .

            if [ ~{main_or_ultralong} -eq 0 ]; then
                bcftools filter --include 'ABS(SVLEN)<50' --output-type v ${SAMPLE_ID}${SUFFIX} --output ${SAMPLE_ID}_lt_50.vcf
                ProcessVcf ${SAMPLE_ID} ${SAMPLE_ID}_lt_50.vcf ${SAMPLE_ID}_cluster_sizes_lt_50.txt ${SAMPLE_ID}_deltas_lt_50.txt ~{chunk_id}_stats_lt_50.csv ~{chunk_id}_different_annotations_lt_50.csv ~{chunk_id}_trs_lt_50.csv
            fi

            bcftools filter --include 'ABS(SVLEN)>=50' --output-type v ${SAMPLE_ID}${SUFFIX} --output ${SAMPLE_ID}_ge_50.vcf
            ProcessVcf ${SAMPLE_ID} ${SAMPLE_ID}_ge_50.vcf ${SAMPLE_ID}_cluster_sizes_ge_50.txt ${SAMPLE_ID}_deltas_ge_50.txt ~{chunk_id}_stats_ge_50.csv ~{chunk_id}_different_annotations_ge_50.csv ~{chunk_id}_trs_ge_50.csv

            rm -f ${SAMPLE_ID}_*.vcf* ${SAMPLE_ID}_*.bcf*
        done 3< chunk.csv
        if [ ~{main_or_ultralong} -eq 0 ]; then
            cat *_cluster_sizes_lt_50.txt > ~{chunk_id}_cluster_sizes_lt_50.txt
            cat *_deltas_lt_50.txt > ~{chunk_id}_deltas_lt_50.txt
        else
            touch ~{chunk_id}_cluster_sizes_lt_50.txt
            touch ~{chunk_id}_deltas_lt_50.txt
        fi
        cat *_cluster_sizes_ge_50.txt > ~{chunk_id}_cluster_sizes_ge_50.txt
        cat *_deltas_ge_50.txt > ~{chunk_id}_deltas_ge_50.txt
    >>>
    
    output {
        File stats_lt_50 = chunk_id + "_stats_lt_50.csv"
        File stats_ge_50 = chunk_id + "_stats_ge_50.csv"
        File cluster_sizes_lt_50 = chunk_id + "_cluster_sizes_lt_50.txt"
        File cluster_sizes_ge_50 = chunk_id + "_cluster_sizes_ge_50.txt"
        File deltas_lt_50 = chunk_id + "_deltas_lt_50.txt"
        File deltas_ge_50 = chunk_id + "_deltas_ge_50.txt"
        File different_annotations_lt_50 = chunk_id + "_different_annotations_lt_50.csv"
        File different_annotations_ge_50 = chunk_id + "_different_annotations_ge_50.csv"
        File n_tr_clusters_lt_50 = chunk_id + "_trs_lt_50.csv"
        File n_tr_clusters_ge_50 = chunk_id + "_trs_ge_50.csv"
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
