version 1.0


# 
#
workflow BuildKanpigHapVcf {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        Int min_sv_length = 20
        
        File windows_bed
        File reference_fa
        String remote_dir
    }
    parameter_meta {
        truvari_collapsed_vcf_gz: "A cohort-level truvari-collapsed VCF."
        min_sv_length: "Keep only calls >=X from the input VCF."
    }
    
    call BuildWindowVcfs {
        input:
            truvari_collapsed_vcf_gz = truvari_collapsed_vcf_gz,
            truvari_collapsed_tbi = truvari_collapsed_tbi,
            min_sv_length = min_sv_length,
            windows_bed = windows_bed,
            reference_fa = reference_fa,
            remote_output_dir = remote_dir+"/chunks"
    }
    call ConcatWindowVcfs {
        input:
            remote_input_dir = remote_dir+"/chunks",
            remote_output_dir = remote_dir+"/concatenated",
            input_flag = BuildWindowVcfs.out_flag
    }
    
    output {
    }
}


# Creates a VCF for each window in a BED file.
#
task BuildWindowVcfs {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        Int min_sv_length
        
        File windows_bed
        File reference_fa
        String remote_output_dir
        
        Int n_cpu = 4
        Int ram_size_gb = 64
        Int disk_size_gb = 100
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
        EFFECTIVE_RAM_GB_PER_THREAD=$(( ( ~{ram_size_gb} - 2 ) / ${N_THREADS} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        
        
        # Remark: the procedure uploads the window BCFs at regular intervals,
        # to use minimal disk.
        #
        function ChunkThread() {
            local THREAD_ID=$1
            local CHUNK_FILE=$2
            local INPUT_BCF=$3
            local FIRST_WINDOW=$4
            
            mkdir -p ./${THREAD_ID}_input_vcfs/
            mkdir -p ./${THREAD_ID}_output_vcfs/
            WINDOW_ID=${FIRST_WINDOW}; COUNTER="0"
            while read WINDOW; do
                CHROM=$(echo ${WINDOW} | cut -d , -f 1)
                START=$(echo ${WINDOW} | cut -d , -f 2)
                END=$(echo ${WINDOW} | cut -d , -f 3)
                echo -e "${CHROM}\t${START}\t${END}" > ${THREAD_ID}.bed
                bcftools view --no-header --regions-file ${THREAD_ID}.bed --regions-overlap pos ${INPUT_BCF} > ./${THREAD_ID}_input_vcfs/chunk_${WINDOW_ID}.vcf
                COUNTER=$(( ${COUNTER} + 1 ))
                if [ ${COUNTER} -eq ${N_FILES_FOR_HAP_VCF} ]; then
                    ls ./${THREAD_ID}_input_vcfs/chunk_*.vcf | sort -V > ${THREAD_ID}_list.txt
                    java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB_PER_THREAD}G BuildKanpigHapVcf ${THREAD_ID}_list.txt reference.fa ${N_SAMPLES} ./${THREAD_ID}_output_vcfs 0 0
                    for VCF_FILE in $(ls ./${THREAD_ID}_output_vcfs/*_haps.vcf); do
                        cat header.txt ${VCF_FILE} | bcftools view --output-type z > ${VCF_FILE}.gz
                        tabix -f ${VCF_FILE}.gz
                        rm -f ${VCF_FILE}
                    done
                    gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ./${THREAD_ID}_output_vcfs/'*_haps.vcf.gz*' ~{remote_output_dir}/
                    gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ./${THREAD_ID}_output_vcfs/'*_records.vcf' ~{remote_output_dir}/
                    gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ./${THREAD_ID}_output_vcfs/'*_map.csv' ~{remote_output_dir}/
                    rm -f ./${THREAD_ID}_input_vcfs/*.vcf
                    COUNTER="0"
                fi
                WINDOW_ID=$(( ${WINDOW_ID} + 1 ))
            done < ${CHUNK_FILE}
            N_CHUNKS=$(ls ./${THREAD_ID}_input_vcfs/*.vcf | wc -l || echo 0)
            if [ ${N_CHUNKS} -gt 0 ]; then
                ls ./${THREAD_ID}_input_vcfs/chunk_*.vcf | sort -V > ${THREAD_ID}_list.txt
                java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB_PER_THREAD}G BuildKanpigHapVcf ${THREAD_ID}_list.txt reference.fa ${N_SAMPLES} ./${THREAD_ID}_output_vcfs 0
                for VCF_FILE in $(ls ./${THREAD_ID}_output_vcfs/*_haps.vcf); do
                    cat header.txt ${VCF_FILE} | bcftools view --output-type z > ${VCF_FILE}.gz
                    tabix -f ${VCF_FILE}.gz
                    rm -f ${VCF_FILE}
                done
                gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ./${THREAD_ID}_output_vcfs/'*_haps.vcf.gz*' ~{remote_output_dir}/
                gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ./${THREAD_ID}_output_vcfs/'*_records.vcf' ~{remote_output_dir}/
                gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ./${THREAD_ID}_output_vcfs/'*_map.csv' ~{remote_output_dir}/
                rm -f ./${THREAD_ID}_input_vcfs/*.vcf
            fi
        }
        
        
        
        
        # -------------------------- Main program ------------------------------
        N_FILES_FOR_HAP_VCF="50"
        
        mv ~{truvari_collapsed_vcf_gz} in.vcf.gz
        mv ~{truvari_collapsed_tbi} in.vcf.gz.tbi
        
        # Formatting the header
        bcftools view --header-only in.vcf.gz > tmp.txt
        N_ROWS=$(wc -l < tmp.txt)
        head -n $(( ${N_ROWS} - 1 )) tmp.txt > header.txt
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> header.txt
        N_SAMPLES=$(tail -n 1 tmp.txt | tr '\t' '\n' | wc -l)
        N_SAMPLES=$(( ${N_SAMPLES} - 9 ))
        
        # Keeping only records in the given length range
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'ABS(SVLEN)>='~{min_sv_length} --output-type b in.vcf.gz > out.bcf
        rm -f in.vcf.gz* ; mv out.bcf in.bcf ; bcftools index in.bcf
        
        # Building window VCFs in parallel
        mv ~{reference_fa} reference.fa
        samtools faidx reference.fa
        bedtools slop -b 1 -i ~{windows_bed} -g reference.fa.fai | tr '\t' ',' > windows.csv
        N_ROWS=$(wc -l < windows.csv)
        N_ROWS_PER_THREAD=$(( ${N_ROWS} / ${N_THREADS} ))
        split -l ${N_ROWS_PER_THREAD} -d -a 4 windows.csv chunk_
        cat chunk_0000
        THREAD_ID="0"; FIRST_WINDOW="0"
        for CHUNK_FILE in $(ls chunk_* | sort -V); do
            ChunkThread ${THREAD_ID} ${CHUNK_FILE} in.bcf ${FIRST_WINDOW} &
            THREAD_ID=$(( ${THREAD_ID} + 1 ))
            N_WINDOWS=$(wc -l < ${CHUNK_FILE})
            FIRST_WINDOW=$(( ${FIRST_WINDOW} + ${N_WINDOWS} ))
        done
        wait
        
        echo "1" > out.txt
    >>>
    
    output {
        File out_flag = "out.txt"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


#
task ConcatWindowVcfs {
    input {
        String remote_input_dir
        String remote_output_dir
        File input_flag
        
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 200
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        date
        gsutil -m cp ~{remote_input_dir}/'*_haps.vcf.gz*' .
        date
        ls *_haps.vcf.gz | sort -V > list.txt
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --output-type z --file-list list.txt > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
        
        gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv out.vcf.'gz*' ~{remote_output_dir}
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
