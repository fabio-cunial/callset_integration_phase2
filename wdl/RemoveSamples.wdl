version 1.0


# Replaces all sample columns with a single sample column where all calls have
# 0/1 GT. Assigns a unique integer ID to every call, moving the original ID to
# the INFO field (this is necessary for kanpig, which otherwise complains
# about duplicated IDs, which arise naturally from the previous steps of the 
# pipeline).
#
workflow RemoveSamples {
    input {
        File intersample_vcf_gz
        File intersample_tbi
    }
    parameter_meta {
    }
    
    call RemoveSamplesImpl {
        input:
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi
    }
    
    output {
        File cleaned_vcf_gz = RemoveSamplesImpl.cleaned_vcf_gz
        File cleaned_tbi = RemoveSamplesImpl.cleaned_tbi
    }
}


task RemoveSamplesImpl {
    input {
        File intersample_vcf_gz
        File intersample_tbi
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 40*ceil(size(intersample_vcf_gz, "GB"))
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        ${TIME_COMMAND} bcftools view --header-only ~{intersample_vcf_gz} > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > cleaned.vcf
        echo '##INFO=<ID=ORIGINAL_ID,Number=1,Type=String,Description="Original ID from truvari collapse">' >> cleaned.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> cleaned.vcf
        bcftools view --no-header ~{intersample_vcf_gz} | awk 'BEGIN { i=0; } { gsub(/;/,"_",$3); printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s;ORIGINAL_ID=%s\tGT\t0/1\n",$1,$2,++i,$4,$5,$6,$7,$8,$3); }' >> cleaned.vcf
        rm -f ~{intersample_vcf_gz}
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} cleaned.vcf
        tabix -f cleaned.vcf.gz
    >>>

    output {
        File cleaned_vcf_gz = work_dir + "/cleaned.vcf.gz"
        File cleaned_tbi = work_dir + "/cleaned.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: 2
        memory: "8GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}