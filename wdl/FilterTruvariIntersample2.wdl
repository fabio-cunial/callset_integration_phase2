version 1.0


# 
#
workflow FilterTruvariIntersample2 {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        
        Array[Int] min_n_samples
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            truvari_collapsed_vcf_gz = truvari_collapsed_vcf_gz,
            truvari_collapsed_tbi = truvari_collapsed_tbi,
            min_n_samples = min_n_samples
    }
    
    output {
        Array[File] filtered_vcf_gz = Impl.filtered_vcf_gz
        Array[File] filtered_tbi = Impl.filtered_tbi
    }
}


#
task Impl {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        
        Array[Int] min_n_samples
    }
    parameter_meta {
    }
    
    Int n_cpus = 4*length(min_n_samples)
    Int ram_size_gb = 8*length(min_n_samples)
    Int disk_size_gb = (length(min_n_samples)+1)*10*ceil(size(truvari_collapsed_vcf_gz, "GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        function filter_thread() {
            local N_SAMPLES=$1
            
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS_PER_VALUE} --include 'COUNT(GT="alt")>='${N_SAMPLES} --output-type z ~{truvari_collapsed_vcf_gz} > ${N_SAMPLES}_tmp.vcf.gz
            tabix -f ${N_SAMPLES}_tmp.vcf.gz
            bcftools view --header-only ${N_SAMPLES}_tmp.vcf.gz > ${N_SAMPLES}_header.txt
            N_ROWS=$(wc -l < ${N_SAMPLES}_header.txt)
            head -n $(( ${N_ROWS} - 1 )) ${N_SAMPLES}_header.txt > ${N_SAMPLES}_filtered.vcf
            echo '##INFO=<ID=ORIGINAL_ID,Number=1,Type=String,Description="Original ID from truvari collapse">' >> ${N_SAMPLES}_filtered.vcf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> ${N_SAMPLES}_filtered.vcf
            bcftools view --threads ${N_THREADS_PER_VALUE} --no-header ${N_SAMPLES}_tmp.vcf.gz | awk 'BEGIN { i=0; } { gsub(/;/,"_",$3); printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s;ORIGINAL_ID=%s\tGT\t0/1\n",$1,$2,++i,$4,$5,$6,$7,$8,$3); }' >> ${N_SAMPLES}_filtered.vcf
            rm -f ${N_SAMPLES}_tmp.vcf.gz*
            ${TIME_COMMAND} bgzip -@ ${N_THREADS_PER_VALUE} ${N_SAMPLES}_filtered.vcf
            ${TIME_COMMAND} tabix -f ${N_SAMPLES}_filtered.vcf.gz
        }
        
        
        # Main program
        echo ~{sep="," min_n_samples} | tr ',' '\n' > min_n_samples.txt
        N_VALUES=$(wc -l < min_n_samples.txt)
        N_THREADS_PER_VALUE=$(( ${N_THREADS} / ${N_VALUES} ))
        while read N_SAMPLES; do
            filter_thread ${N_SAMPLES} &
        done < min_n_samples.txt
        wait
        ls -laht
    >>>
    
    output {
        Array[File] filtered_vcf_gz = glob("*_filtered.vcf.gz")
        Array[File] filtered_tbi = glob("*_filtered.vcf.gz.tbi")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
