version 1.0


#
workflow GetMapqDistribution {
    input {
        File alignments_bam
        File alignments_bai
        
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File pav_vcf_gz
        File pav_vcf_gz_tbi
        
        Int clustering_distance_bp = 500
        
        Int n_cpu
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    call MapqImpl {
        input:
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            
            pbsv_vcf_gz = pbsv_vcf_gz,
            pbsv_vcf_gz_tbi = pbsv_vcf_gz_tbi,
            sniffles_vcf_gz = sniffles_vcf_gz,
            sniffles_vcf_gz_tbi = sniffles_vcf_gz_tbi,
            pav_vcf_gz = pav_vcf_gz,
            pav_vcf_gz_tbi = pav_vcf_gz_tbi,
            
            clustering_distance_bp = clustering_distance_bp,
            
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb
    }
    
    output {
    	File mapq_txt = MapqImpl.mapq_txt
    }
}


#
task MapqImpl {
    input {
        File alignments_bam
        File alignments_bai
        
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File pav_vcf_gz
        File pav_vcf_gz_tbi
        
        Int clustering_distance_bp
        
        Int n_cpu
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = ceil(size(alignments_bam,"GB")) + 50
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} bedtools merge -i ~{pbsv_vcf_gz} -i ~{sniffles_vcf_gz} -i ~{pav_vcf_gz} -d ~{clustering_distance_bp} > merged.bed
        ${TIME_COMMAND} samtools view -@ ${N_THREADS} --regions-file merged.bed ~{alignments_bam} | awk '{print $5}' | awk 'BEGIN {OFS="\t"} { count[$1]++ } END { for (mapq in count) print mapq, count[mapq] }' > mapq.txt
    >>>
    
    output {
    	File mapq_txt = "mapq.txt"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
