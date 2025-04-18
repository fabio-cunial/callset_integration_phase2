version 1.0


#
workflow InterCenterBench {
    input {
        Array[File] center_vcf_gz
        Array[File] center_tbi
    }
    parameter_meta {
    }

    scatter (chr in chromosomes) {
        call Merge {
            input:
                chromosome = chr,
                center_vcf_gz = center_vcf_gz,
                center_tbi = center_tbi
        }
    }
    call ConcatenateChromosomes {
        input:
            chromosomes_vcf_gz = Merge.merged_vcf_gz,
            chromosomes_tbi = Merge.merged_tbi
    }
    
    output {
        File vcf_gz = ConcatenateChromosomes.vcf_gz
        File vcf_gz_tbi = ConcatenateChromosomes.vcf_gz_tbi
    }
}


#
task Merge {
    input {
        String chromosome
        
        Array[File] center_vcf_gz
        Array[File] center_tbi
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Filtering by chromosome
        INPUT_FILES=~{sep=',' center_vcf_gz}
        echo ${INPUT_FILES} | tr ',' '\n' | sort > list.txt
        rm -f list_filtered.txt
        i="0"
        while read FILE; do
            i=$(( ${i} + 1 ))
            bcftools view --header-only ${FILE} > header.txt
            N_ROWS=$(wc -l < header.txt)
            head -n $(( ${N_ROWS} - 1 )) header.txt > ${i}_~{chromosome}.vcf
            tail -n 1 header.txt | awk '{ \
                printf("%s",$1); \
                for (i=2; i<=9; i++) printf("\t%s",$i); \
                printf("\tSAMPLE\n"); \
            }' >> ${i}_~{chromosome}.vcf
            bcftools view --threads ${N_THREADS} --no-header --output-type v ${FILE} ~{chromosome} | awk '{ \
                printf("%s",$1); \
                for (i=2; i<=8; i++) printf("\t%s",$i); \
                printf("\tGT\t0/1\n"); \
            }' >> ${i}_~{chromosome}.vcf
            bgzip -@ ${N_THREADS} --compress-level 1 ${i}_~{chromosome}.vcf
            tabix -f ${i}_~{chromosome}.vcf.gz
            echo ${i}_~{chromosome}.vcf.gz >> list_filtered.txt
            rm -f ${FILE}*
        done < list.txt
        rm -f list.txt
        mv list_filtered.txt list.txt
        
        # BCFTOOLS
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type z > ~{chromosome}.merged.vcf.gz
        tabix -f ~{chromosome}.merged.vcf.gz
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type z ~{chromosome}.merged.vcf.gz > ~{chromosome}.normed.vcf.gz
        tabix -f ~{chromosome}.normed.vcf.gz
        rm -f ~{chromosome}.merged.vcf.gz
        
        # TRUVARI
        ${TIME_COMMAND} truvari collapse --input ~{chromosome}.normed.vcf.gz --sizemin 0 --sizemax 1000000 --keep common --gt off --output tmp.vcf
        ${TIME_COMMAND} bcftools sort --max-mem $(( ~{ram_size_gb} - 4 ))G --output-type z tmp.vcf > ~{chromosome}.collapsed.vcf.gz
        tabix -f ~{chromosome}.collapsed.vcf.gz
    >>>
    
    output {
        File merged_vcf_gz = work_dir + "/" + chromosome + ".collapsed.vcf.gz"
        File merged_tbi = work_dir + "/" + chromosome + ".collapsed.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}


#
task ConcatenateChromosomes {
    input {
        Array[File] chromosomes_vcf_gz
        Array[File] chromosomes_tbi
        
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 256  # Arbitrary

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        INPUT_FILES=~{sep=',' chromosomes_vcf_gz}
        echo ${INPUT_FILES} | tr ',' '\n' | sort > list.txt
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --file-list list.txt --output-type z > concat.vcf.gz
        tabix -f concat.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/concat.vcf.gz"
        File vcf_gz_tbi = work_dir + "/concat.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
