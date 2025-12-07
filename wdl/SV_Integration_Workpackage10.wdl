version 1.0


# Merging by ID the re-genotyped personalized BCFs of all samples, in a given 
# chromosome chunk.
#
workflow SV_Integration_Workpackage10 {
    input {
        Int chunk_id
        File sample_ids
        Array[String] bi_samples_to_prefer_over_ha
        
        String remote_indir_bi
        String remote_indir_ha
        String remote_indir_bcm
        String remote_indir_uw
        String remote_indir_controls_15x
        String remote_indir_controls_30x
        
        Int n_expected_samples_bi
        Int n_expected_samples_ha
        Int n_expected_samples_bcm
        Int n_expected_samples_uw
        Int n_expected_samples_controls_15x
        Int n_expected_samples_controls_30x
        
        String remote_outdir
    }
    parameter_meta {
        sample_ids: "Speficies the order of the samples used by bcftools merge."
        bi_samples_to_prefer_over_ha: "If a sample occurs both in BI and HA, the HA one is usually preferred, except for these samples."
        remote_indir_bi: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            chunk_id = chunk_id,
            sample_ids = sample_ids,
            bi_samples_to_prefer_over_ha = bi_samples_to_prefer_over_ha,
            remote_indir_bi = remote_indir_bi,
            remote_indir_ha = remote_indir_ha,
            remote_indir_bcm = remote_indir_bcm,
            remote_indir_uw = remote_indir_uw,
            remote_indir_controls_15x = remote_indir_controls_15x,
            remote_indir_controls_30x = remote_indir_controls_30x,
            n_expected_samples_bi = n_expected_samples_bi,
            n_expected_samples_ha = n_expected_samples_ha,
            n_expected_samples_bcm = n_expected_samples_bcm,
            n_expected_samples_uw = n_expected_samples_uw,
            n_expected_samples_controls_15x = n_expected_samples_controls_15x,
            n_expected_samples_controls_30x = n_expected_samples_controls_30x,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, one 10 MB chunk of chr6:
#
# TOOL               CPU     RAM     TIME
# bcftools merge     120%    23G     30m
#
# Peak disk usage: ??????
#
task Impl {
    input {
        Int chunk_id
        File sample_ids
        Array[String] bi_samples_to_prefer_over_ha
        
        String remote_indir_bi
        String remote_indir_ha
        String remote_indir_bcm
        String remote_indir_uw
        String remote_indir_controls_15x
        String remote_indir_controls_30x
        
        Int n_expected_samples_bi
        Int n_expected_samples_ha
        Int n_expected_samples_bcm
        Int n_expected_samples_uw
        Int n_expected_samples_controls_15x
        Int n_expected_samples_controls_30x
        
        String remote_outdir
        
        Int n_cpu = 16
        Int ram_size_gb = 128
        Int disk_size_gb = 50
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"

        
        # Checking that every input dataset has the expected number of samples
        # in the chunk.
        N_FILES=$(gsutil ls ~{remote_indir_bi}/'*_chunk_'~{chunk_id}.bcf | wc -l)
        if [ ${N_FILES} -ne ~{n_expected_samples_bi} ]; then
            echo "ERROR: BI has ${N_FILES} files != ~{n_expected_samples_bi}"
            exit
        fi
        N_FILES=$(gsutil ls ~{remote_indir_ha}/'*_chunk_'~{chunk_id}.bcf | wc -l)
        if [ ${N_FILES} -ne ~{n_expected_samples_ha} ]; then
            echo "ERROR: HA has ${N_FILES} files != ~{n_expected_samples_ha}"
            exit
        fi
        N_FILES=$(gsutil ls ~{remote_indir_bcm}/'*_chunk_'~{chunk_id}.bcf | wc -l)
        if [ ${N_FILES} -ne ~{n_expected_samples_bcm} ]; then
            echo "ERROR: BCM has ${N_FILES} files != ~{n_expected_samples_bcm}"
            exit
        fi
        N_FILES=$(gsutil ls ~{remote_indir_uw}/'*_chunk_'~{chunk_id}.bcf | wc -l)
        if [ ${N_FILES} -ne ~{n_expected_samples_uw} ]; then
            echo "ERROR: UW has ${N_FILES} files != ~{n_expected_samples_uw}"
            exit
        fi
        N_FILES=$(gsutil ls ~{remote_indir_controls_15x}/'*_chunk_'~{chunk_id}.bcf | wc -l)
        if [ ${N_FILES} -ne ~{n_expected_samples_controls_15x} ]; then
            echo "ERROR: CONTROLS_15X has ${N_FILES} files != ~{n_expected_samples_controls_15x}"
            exit
        fi
        N_FILES=$(gsutil ls ~{remote_indir_controls_30x}/'*_chunk_'~{chunk_id}.bcf | wc -l)
        if [ ${N_FILES} -ne ~{n_expected_samples_controls_30x} ]; then
            echo "ERROR: CONTROLS_30X has ${N_FILES} files != ~{n_expected_samples_controls_30x}"
            exit
        fi
        
        # Localizing all the samples for the given chunk, and handling samples
        # that occur in multiple input datasets.
        mkdir ./input_bcfs/
        gsutil -m cp ~{remote_indir_bi}/'*'_chunk_~{chunk_id}.'bcf*' ./input_bcfs/
        gsutil -m cp ~{remote_indir_ha}/'*'_chunk_~{chunk_id}.'bcf*' ./input_bcfs/
        for SAMPLE_ID in ~{bi_samples_to_prefer_over_ha}; do
            gsutil -m cp ~{remote_indir_bi}/${SAMPLE_ID}_chunk_~{chunk_id}.'bcf*' ./input_bcfs/
        done
        gsutil -m cp ~{remote_indir_uw}/'*.bcf*' ./input_bcfs/
        gsutil -m cp ~{remote_indir_bcm}/'*.bcf*' ./input_bcfs/
        gsutil -m cp ~{remote_indir_controls_15x}/'*.bcf*' ./input_bcfs/
        gsutil -m cp ~{remote_indir_controls_30x}/'*.bcf*' ./input_bcfs/
        N_DOWNLOADED_SAMPLES=$(ls ./input_bcfs/*_chunk_~{chunk_id}.bcf | wc -l)
        N_SAMPLES=$(cat ~{sample_ids} | wc -l)
        if [ ${N_DOWNLOADED_SAMPLES} -ne ${N_SAMPLES} ]; then
            echo "Error: the number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is different from the number of samples specified (${N_SAMPLES})."
            exit 1
        fi
        df -h
        
        # Merging records by ID, since the records in every BCF originate from
        # the same cohort BCF, which had distinct IDs.
        while read SAMPLE_ID; do
            echo ./input_bcfs/${SAMPLE_ID}_chunk_~{chunk_id}.bcf >> list.txt
        done < ~{sample_ids}
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge id --output-type b --file-list list.txt > chunk_~{chunk_id}.bcf
        bcftools index chunk_~{chunk_id}.bcf
        ls -laht
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp chunk_~{chunk_id}.'bcf*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading merged BCF. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
