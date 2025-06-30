version 1.0


#
workflow BuildTrainingResourceIntersample {
    input {
        Array[File] intrasample_vcf_gz
    }
    parameter_meta {
    }
    
    call BuildTrainingResourceIntersampleImpl {
        input:
            intrasample_vcf_gz = intrasample_vcf_gz
    }
    
    output {
        File merged_vcf = BuildTrainingResourceIntersampleImpl.merged_vcf
        File merged_tbi = BuildTrainingResourceIntersampleImpl.merged_tbi
    }
}


# Simply runs `bcftools merge` and forces the output to have a single sample
# with all calls `0|0`.
#
task BuildTrainingResourceIntersampleImpl {
    input {
        Array[File] intrasample_vcf_gz
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/mnt/disks/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        INPUT_FILES=~{sep=',' intrasample_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            tabix -f ${INPUT_FILE}
            echo ${INPUT_FILE} >> list.txt
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list.txt --output-type v > tmp1.vcf
        bcftools view --header-only tmp1.vcf > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > merged.vcf
        tail -n 1 header.txt | awk '{ \
            printf("%s",$1); \
            for (i=2; i<=9; i++) printf("\t%s",$i); \
            printf("\tDUMMY_SAMPLE_NAME\n"); \
        }' >> merged.vcf
        rm -f header.txt
        bcftools view --no-header tmp1.vcf | awk '{ \
            printf("%s",$1); \
            for (i=2; i<=8; i++) printf("\t%s",$i); \
            printf("\tGT\t0|0\n"); \
        }' >> merged.vcf
        rm -f tmp1.vcf
        bgzip --threads ${N_THREADS} --compress-level 1 merged.vcf
        tabix -f merged.vcf.gz
    >>>

    output {
        File merged_vcf = work_dir + "/merged.vcf.gz"
        File merged_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: 16
        memory: "32GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
