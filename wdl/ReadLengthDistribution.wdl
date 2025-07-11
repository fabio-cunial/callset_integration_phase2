version 1.0


#
# Performance on HPRC's HG002 8x. 16 physical cores, 32 GB RAM, 500 GB HDD. 
# Total time: 3h 30m
# Total cost: ????
#
# STEP                  TIME            CPU            RAM
# samtools fastq        15 m             33 %         16 MB
# split                 30 m              3 %          2 MB
# samtools fqidx        10 m              1 %         20 MB
#
workflow ReadLengthDistribution {
    input {
        String bam_address
        String billing_project = "broad-firecloud-dsde-methods"
        
        Int n_cores = 4
        Int mem_gb = 32
        Int disk_size_gb = 1000
    }
    parameter_meta {
        bam_address: "Can be .bam, .fastq, .fastq.gz"
    }
    
    call DistributionImpl {
        input:
            bam_address = bam_address,
            billing_project = billing_project,
            n_cores = n_cores,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
        File out_fai = DistributionImpl.out_lengths
    }
}


task DistributionImpl {
    input {
        String bam_address
        String billing_project
        
        Int n_cores
        Int mem_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        function getLengths() {
            local FASTQ_FILE=$1
            
            ${TIME_COMMAND} samtools fqidx --fai-idx ${FASTQ_FILE}.fai ${FASTQ_FILE}
            cut -f 2 ${FASTQ_FILE}.fai > ${FASTQ_FILE}.lengths
            # sort --parallel 1 --mmap --numeric-sort
            rm -f ${FASTQ_FILE}.fai
        }
        
        
        # Main program
        
        # Downloading and converting to FASTQ
        if [[ ~{bam_address} == gs://* ]]; then
            gsutil -u ~{billing_project} -m cp ~{bam_address} .
        elif [[ ~{bam_address} == s3://* ]]; then
            aws s3 --no-sign-request cp ~{bam_address} .
        else
            wget ~{bam_address}
        fi
        df -h
        FILE_NAME=$(basename ~{bam_address})
        if [[ ${FILE_NAME} == *.bam ]]; then
            ${TIME_COMMAND} samtools fastq -@ ${N_THREADS} -n ${FILE_NAME} > reads.fastq 
        elif [[ ${FILE_NAME} == *.fastq.gz ]]; then
            ${TIME_COMMAND} gunzip -c ${FILE_NAME} > reads.fastq
        elif [[ ${FILE_NAME} == *.fastq ]]; then
            mv ${FILE_NAME} reads.fastq
        fi
        df -h
        rm -f ${FILE_NAME}
        
        # Splitting into chunks
        N_LINES=$(wc -l < reads.fastq)
        N_READS=$(( ${N_LINES} / 4 ))
        N_LINES=$(( (${N_READS} / ${N_THREADS})*4 ))
        ${TIME_COMMAND} split -l ${N_LINES} -a 2 -d reads.fastq chunk_
        df -h
        rm -f reads.fastq
        for FILE in $(ls chunk_*); do
            mv ${FILE} ${FILE}.fastq
            getLengths ${FILE}.fastq &
        done
        wait
        date
        cat *.lengths > lengths.txt
        # sort --merge --parallel ${N_THREADS} --mmap --numeric-sort --batch-size=$(( ${N_THREADS} + 1 )) *.lengths > lengths.txt
        date
    >>>
    
    output {
        File out_lengths = work_dir + "/lengths.txt"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
