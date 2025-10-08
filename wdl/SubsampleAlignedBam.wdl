version 1.0


#
workflow SubsampleAlignedBam {
    input {
        String sample_id
        
        File aligned_bam
        File aligned_bai
        Float aligned_bam_coverage
        
        String target_coverages = "15,30"
        String remote_output_dir
    }
    parameter_meta {
        target_coverages: "Comma-separated. Each coverage is sampled independently, so there is no guarantee that a bigger-coverage sample contains the alignments in a smaller-coverage sample."
        n_cores: ">=1 per coverage"
    }
    
    call SubsampleImpl {
        input:
            sample_id = sample_id,
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            aligned_bam_coverage = aligned_bam_coverage,
            target_coverages = target_coverages,
            remote_output_dir = remote_output_dir
    }
    
    output {
    }
}


# Performance with 16 cores and 32GB of RAM on a 66x, 158GB BAM:
#
# COMMAND           CPU     RAM     TIME
# samtools view     400%    30M     20m       // regardless of sampling factor
#
task SubsampleImpl {
    input {
        String sample_id
        
        File aligned_bam
        File aligned_bai
        Float aligned_bam_coverage
        
        String target_coverages
        String remote_output_dir

        Int n_cores = 8
        Int mem_gb = 8
    }
    parameter_meta {
        target_coverages: "Comma-separated. Each coverage is sampled independently, so there is no guarantee that a bigger-coverage sample contains the alignments in a smaller-coverage sample."
        n_cores: ">=1 per coverage"
    }
    
    Int disk_size_gb = 10*( ceil(size(aligned_bam,"GB")) )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        function subsample() {
            local COVERAGE=$1
            local PREFIX=$2
            local NT=$3
            
            FACTOR=$(echo "scale=2; ${COVERAGE}/~{aligned_bam_coverage} " | bc)
            ${TIME_COMMAND} samtools view -@ ${NT} --subsample ${FACTOR} --bam ~{aligned_bam} > ${PREFIX}_${COVERAGE}x.bam
            samtools index ${PREFIX}_${COVERAGE}x.bam
        }
        
        
        # Subsampling in parallel
        echo ~{target_coverages} | tr ',' '\n' > coverages.txt
        N_COVERAGES=$(wc -l < coverages.txt)
        N_THREADS_PER_COVERAGE=$(( ${N_THREADS} / ${N_COVERAGES} ))
        date
        while read COVERAGE; do
            subsample ${COVERAGE} ~{sample_id} ${N_THREADS_PER_COVERAGE} &
        done < coverages.txt
        wait
        date
        ls -laht
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{sample_id}_'*x.bam*' ~{remote_output_dir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading the BAMs. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 2
    }
}
