version 1.0


# Builds a BND truthset using diploid-assembly-to-reference BAMs
#
workflow SV_Integration_BndBuildTruth {
    input {
        String sample_id
        File hap1_bam
        File hap2_bam

        Int max_adjacency_distance = 1000
        Int min_violation_distance = 50000
        Int output_mode = 0
        File? header_vcf

        String remote_outdir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
        max_adjacency_distance: "Max distance (on an assembled contig) between two alignments for them to be considered adjacent. 1kbp seems a good value based on a histogram of nearest-neighbor distances."
        min_violation_distance: "Min distance (on the same reference chr) between two alignments (that are adjacent on some contig) for them to be considered a colinearity violation. Any setting will capture some ultralong DELs."
        output_mode: "0=a CSV of points; 1=a BND VCF"
        header_vcf: "Header to be given to the output VCF, if any."
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            sample_id = sample_id,
            hap1_bam = hap1_bam,
            hap2_bam = hap2_bam,

            max_adjacency_distance = max_adjacency_distance,
            min_violation_distance = min_violation_distance,
            output_mode = output_mode,
            header_vcf = header_vcf,

            remote_outdir = remote_outdir,

            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 2-core, 8GB VM:
#
# TOOL                                      CPU%        RAM         TIME
# samtools sort
# AssemblySam2Breakpoints2
#
task Impl {
    input {
        String sample_id
        File hap1_bam
        File hap2_bam

        Int max_adjacency_distance
        Int min_violation_distance
        Int output_mode
        File? header_vcf

        String remote_outdir
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 20
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
        EFFECTIVE_RAM_MB=$(( (~{ram_size_gb} - 1) * 1024 ))
        RAM_PER_THREAD_MB=$(( ${EFFECTIVE_RAM_MB} / ${N_THREADS} ))


        samtools --version 1>&2
        df -h 1>&2


        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} -n -O SAM -o hap1.sam ~{hap1_bam}
        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} -n -O SAM -o hap2.sam ~{hap2_bam}
        if [ ~{output_mode} -eq 0 ]; then
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${RAM_PER_THREAD_MB}M AssemblySam2Breakpoints2 hap1.sam ~{max_adjacency_distance} ~{min_violation_distance} 0 > ${SAMPLE_ID}_breakpoints1.csv &
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${RAM_PER_THREAD_MB}M AssemblySam2Breakpoints2 hap2.sam ~{max_adjacency_distance} ~{min_violation_distance} 0 > ${SAMPLE_ID}_breakpoints2.csv &
            wait
            cat ${SAMPLE_ID}_breakpoints1.csv ${SAMPLE_ID}_breakpoints2.csv | sort -t , -k1,1 -k2,2n > ${SAMPLE_ID}_breakpoints.csv
            gcloud storage mv ${SAMPLE_ID}_breakpoints.csv ~{remote_outdir}/
        else
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${RAM_PER_THREAD_MB}M AssemblySam2Breakpoints2 hap1.sam ~{max_adjacency_distance} ~{min_violation_distance} 1 > ${SAMPLE_ID}_breakpoints1.vcf &
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${RAM_PER_THREAD_MB}M AssemblySam2Breakpoints2 hap2.sam ~{max_adjacency_distance} ~{min_violation_distance} 1 > ${SAMPLE_ID}_breakpoints2.vcf &
            wait
            cat ~{header_vcf} ${SAMPLE_ID}_breakpoints1.vcf ${SAMPLE_ID}_breakpoints2.vcf > ${SAMPLE_ID}_breakpoints.vcf
            gcloud storage mv ${SAMPLE_ID}_breakpoints.vcf ~{remote_outdir}/
        fi
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
