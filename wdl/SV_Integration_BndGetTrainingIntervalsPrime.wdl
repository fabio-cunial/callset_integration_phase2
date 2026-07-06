version 1.0


# Like `SV_Integration_BndGetTrainingIntervals.wdl`, but uses SAM-derived 
# breakpoints rather than SVIM-asm's BND calls as a source of truth.
#
workflow SV_Integration_BndGetTrainingIntervalsPrime {
    input {
        File samples_tsv

        Int sam_min_alignment_distance = 1000
        Int sam_min_clip_length = 1000

        Int breakpoint_max_distance = 500
        Int breakpoint_filter_mode = 2
        File reference_agp
        
        String remote_indir_query
        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        samples_tsv: "Format: ID, HAP1_BAM, HAP2_BAM"
        sam_min_alignment_distance: "Min distance between an assembly-to-ref alignment and the closest one to call a breakpoint"
        sam_min_clip_length: "Min clip length of an assembly-to-ref alignment to call a breakpoint"
        breakpoint_max_distance: "Max distance between a BND and an assembly breakpoint to mark TPs"
        breakpoint_filter_mode: "1=at least one side of the BND must be close to an assembly breakpoint; 2=both sides of the BND must be close to an assembly breakpoint."
        reference_agp: "Reference AGP file."
        remote_indir_query: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            samples_tsv = samples_tsv,

            sam_min_alignment_distance = sam_min_alignment_distance,
            sam_min_clip_length = sam_min_clip_length,

            breakpoint_max_distance = breakpoint_max_distance,
            breakpoint_filter_mode = breakpoint_filter_mode,
            reference_agp = reference_agp,

            remote_indir_query = remote_indir_query,
            remote_outdir = remote_outdir,

            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 1-core, 4GB VM:
#
# TOOL                                      CPU%        RAM         TIME
# AssemblySam2Breakpoints
# BndFilterWithAssemblyBreakpoints
#
task Impl {
    input {
        File samples_tsv

        Int sam_min_alignment_distance
        Int sam_min_clip_length

        Int breakpoint_max_distance
        Int breakpoint_filter_mode
        File reference_agp
        
        String remote_indir_query
        String remote_outdir
        
        String docker_image
        Int n_cpu = 1
        Int ram_size_gb = 4
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
        EFFECTIVE_RAM_MB=$(( ~{ram_size_gb} * 1024 - 500 ))



        # ----------------------- Steps of the pipeline ------------------------

        function BuildAssemblyBreakpoints() {
            SAMPLE_ID=$1
            REMOTE_HAP1_BAM=$2
            REMOTE_HAP2_BAM=$3
            RAM_SIZE_MB=$4

            rm -f ${SAMPLE_ID}_breakpoints.csv

            ${TIME_COMMAND} gcloud storage cp ${REMOTE_HAP1_BAM} ./hap1.bam
            samtools view hap1.bam > hap1.sam
            rm -f hap1.bam
            java -cp ~{docker_dir} -Xmx${RAM_SIZE_MB}M AssemblySam2Breakpoints hap1.sam ~{sam_min_alignment_distance} ~{sam_min_clip_length} > ${SAMPLE_ID}_breakpoints.csv
            rm -f hap1.sam

            ${TIME_COMMAND} gcloud storage cp ${REMOTE_HAP2_BAM} ./hap2.bam
            samtools view hap2.bam > hap2.sam
            rm -f hap2.bam
            java -cp ~{docker_dir} -Xmx${RAM_SIZE_MB}M AssemblySam2Breakpoints hap2.sam ~{sam_min_alignment_distance} ~{sam_min_clip_length} >> ${SAMPLE_ID}_breakpoints.csv
            rm -f hap2.sam
        }




        # --------------------------- Main program -----------------------------

        samtools --version 1>&2
        bcftools --version 1>&2
        df -h 1>&2

        cat ~{samples_tsv} | tr '\t' ',' > samples.csv
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            REMOTE_HAP1_BAM=$(echo ${LINE} | cut -d , -f 2)
            REMOTE_HAP2_BAM=$(echo ${LINE} | cut -d , -f 3)
            
            # Skipping the sample if it has already been processed
            TEST=$( gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "1" )
            if [ "${TEST}" != "1" ]; then
                continue
            fi

            # Filtering
            BuildAssemblyBreakpoints ${SAMPLE_ID} ${REMOTE_HAP1_BAM} ${REMOTE_HAP2_BAM} ${EFFECTIVE_RAM_MB}
            gcloud storage cp ~{remote_indir_query}/${SAMPLE_ID}_bnd.vcf.'gz*' .
            java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M BndCanonize ${SAMPLE_ID}_bnd.vcf.gz > ${SAMPLE_ID}_bnd_canonized.vcf
            rm -f ${SAMPLE_ID}_bnd.vcf.gz*
            java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_MB}M BndFilterWithAssemblyBreakpoints ${SAMPLE_ID}_bnd_canonized.vcf ${SAMPLE_ID}_breakpoints.csv ~{breakpoint_max_distance} ~{reference_agp} ~{breakpoint_filter_mode} | bcftools sort - --output-type z > ${SAMPLE_ID}_bnd_training.vcf.gz
            rm -f ${SAMPLE_ID}_bnd_canonized.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_bnd_training.vcf.gz
            
            # Uploading and deallocating the sample
            gcloud storage mv ${SAMPLE_ID}_bnd_training.vcf.'gz*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
            rm -rf ${SAMPLE_ID}_*
            ls -laht 1>&2
        done 3< samples.csv
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
