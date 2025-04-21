version 1.0


#
workflow InterCenterBench {
    input {
        File center1_vcf_gz
        File center1_tbi
        File center2_vcf_gz
        File center2_tbi
        File center3_vcf_gz
        File center3_tbi
    }
    parameter_meta {
    }

    call Bench as Bench1 {
        input:
            center1_vcf_gz = center1_vcf_gz,
            center1_tbi = center1_tbi,
            center2_vcf_gz = center2_vcf_gz,
            center2_tbi = center2_tbi
    }
    call Bench as Bench2 {
        input:
            center1_vcf_gz = center1_vcf_gz,
            center1_tbi = center1_tbi,
            center2_vcf_gz = center3_vcf_gz,
            center2_tbi = center3_tbi
    }
    call Bench as Bench3 {
        input:
            center1_vcf_gz = center2_vcf_gz,
            center1_tbi = center2_tbi,
            center2_vcf_gz = center3_vcf_gz,
            center2_tbi = center3_tbi
    }
    
    output {
    }
}


#
task Bench {
    input {
        File center1_vcf_gz
        File center1_tbi
        File center2_vcf_gz
        File center2_tbi
        
        Int n_cpu = 4
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
        
        # Forces exactly one sample with all GTs equal to 0/1
        function filter() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF=$2
            
            bcftools view --header-only ${INPUT_VCF_GZ} > header.txt
            N_ROWS=$(wc -l < header.txt)
            head -n $(( ${N_ROWS} - 1 )) header.txt > ${OUTPUT_VCF}
            tail -n 1 header.txt | awk '{ \
                printf("%s",$1); \
                for (i=2; i<=9; i++) printf("\t%s",$i); \
                printf("\tSAMPLE\n"); \
            }' >> ${OUTPUT_VCF}
            bcftools view --threads ${N_THREADS} --no-header --output-type v ${INPUT_VCF_GZ} | awk '{ \
                printf("%s",$1); \
                for (i=2; i<=8; i++) printf("\t%s",$i); \
                printf("\tGT\t0/1\n"); \
            }' >> ${OUTPUT_VCF}
            bgzip -@ ${N_THREADS} --compress-level 1 ${OUTPUT_VCF}
            tabix -f ${OUTPUT_VCF}.gz
        }
        
        # Main program
        filter ~{center1_vcf_gz} center1.vcf
        filter ~{center2_vcf_gz} center2.vcf
        ${TIME_COMMAND} truvari bench -b center1.vcf.gz -c center2.vcf.gz -o ./truvari/
        mv ./truvari/summary.json .
    >>>
    
    output {
        File summary = work_dir + "/summary.json"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}
