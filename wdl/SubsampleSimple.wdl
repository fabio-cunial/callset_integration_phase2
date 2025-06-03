version 1.0


# Simpler variant with just one BAM in input.
#
workflow SubsampleSimple {
    input {
        String sample_id
        String bam_address
        String coverages = "4,8,16,32"
        String remote_dir
        String billing_project = "broad-firecloud-dsde-methods"
        Int haploid_genome_length_gb = 3
        Int n_cores = 8
        Int mem_gb = 32
        Int disk_size_gb = 500
    }
    parameter_meta {
        bam_address: "Can be .bam, .fastq, .fastq.gz"
        coverages: "Comma-separated"
        remote_dir: "Output directory in a remote bucket"
        n_cores: "At least one per coverage"
    }
    
    call SubsampleImpl {
        input:
            sample_id = sample_id,
            bam_address = bam_address,
            coverages = coverages,
            remote_dir = remote_dir,
            billing_project = billing_project,
            haploid_genome_length_gb = haploid_genome_length_gb,
            n_cores = n_cores,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}


task SubsampleImpl {
    input {
        String sample_id
        String bam_address
        String coverages
        String remote_dir
        String billing_project
        Int haploid_genome_length_gb
        Int n_cores
        Int mem_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TOTAL_RAM_KB=$(cat /proc/meminfo | grep MemTotal | awk '{print $2}')
        TOTAL_RAM_GB=$(( ${TOTAL_RAM_KB} / 1000000 ))
        RAM_PER_THREAD_BYTES=$(( (1000000000*( ${TOTAL_RAM_GB} - 10)) / ${N_THREADS} ))
        HAPLOID_GENOME_LENGTH_GB=$(( ~{haploid_genome_length_gb} * 1000000000 ))
        df -h
        
        MAX_COVERAGE=~{coverages}
        MAX_COVERAGE=${MAX_COVERAGE##*,}
        
        # 1. Downloading the BAM
        if [[ ~{bam_address} == gs://* ]]; then
            SUCCESS=$(gsutil -u ~{billing_project} -m cp ~{bam_address} . && echo 1 || echo 0)
        else
            SUCCESS=$(wget ~{bam_address} && echo 1 || echo 0)
        fi
        if [[ ${SUCCESS} -eq 0 ]]; then
            echo "Error downloading file"
            return 1
        fi
        FILE_NAME=$(basename ~{bam_address})
        if [[ ${FILE_NAME} == *.bam ]]; then
            ${TIME_COMMAND} samtools fastq -@ ${N_THREADS} -n ${FILE_NAME} | pigz --processes ${N_THREADS} --fast --to-stdout > tmp1.fastq.gz
        elif [[ ${FILE_NAME} == *.fastq.gz ]]; then
            mv ${FILE_NAME} tmp1.fastq.gz
        elif [[ ${FILE_NAME} == *.fastq ]]; then
            ${TIME_COMMAND} pigz --processes ${N_THREADS} --fast --to-stdout ${FILE_NAME} > tmp1.fastq.gz
        fi
        df -h
        rm -f ${FILE_NAME}
        ${TIME_COMMAND} ~{docker_dir}/seqkit stats --threads ${N_THREADS} --tabular tmp1.fastq.gz > stats.txt
        TOTAL_N_CHARS=$(cut -f 5 stats.txt | tail -n 1)
        TOTAL_N_READS=$(cut -f 4 stats.txt | tail -n 1)
        
        # 2. Subsampling
        function sample() {
            local INPUT_FASTQ_GZ=$1
            local INPUT_FASTQ_GZ_N_CHARS=$2
            local COVERAGE=$3
            
            TARGET_N_CHARS=$(( ${COVERAGE} * ${HAPLOID_GENOME_LENGTH_GB} ))
            if [ ${TARGET_N_CHARS} -ge ${INPUT_FASTQ_GZ_N_CHARS} ]; then
                echo "WARNING: the remote files contain just ${INPUT_FASTQ_GZ_N_CHARS} bps in total, but we need ${TARGET_N_CHARS} bps to achieve coverage ${COVERAGE}."
                cp ${INPUT_FASTQ_GZ} ~{sample_id}_${COVERAGE}.fastq.gz
            else
                PROPORTION=$( bc <<< "scale=2; ${TARGET_N_CHARS}/${INPUT_FASTQ_GZ_N_CHARS}" )
                ${TIME_COMMAND} ~{docker_dir}/seqkit sample --proportion ${PROPORTION} -o ~{sample_id}_${COVERAGE}.fastq.gz ${INPUT_FASTQ_GZ}
            fi
        }
        echo ~{coverages} | tr ',' '\n' > coverages.txt
        while read COVERAGE; do
            sample tmp1.fastq.gz ${TOTAL_N_CHARS} ${COVERAGE} &
        done < coverages.txt
        wait
        df -h
        rm -f tmp1.fastq.gz
        
        # 3. Uploading
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}_'*.fastq.gz' ~{remote_dir} && echo 0 || echo 1)
            if [[ ${TEST} -eq 1 ]]; then
                echo "Error uploading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/hapestry_experiments"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
