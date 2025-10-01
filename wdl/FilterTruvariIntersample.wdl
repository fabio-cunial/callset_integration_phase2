version 1.0


# 
#
workflow FilterTruvariIntersample {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        
        Array[String] ids
        Array[Int] min_depth
        Array[Int] max_depth
        Array[Int] min_alt_reads
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            truvari_collapsed_vcf_gz = truvari_collapsed_vcf_gz,
            truvari_collapsed_tbi = truvari_collapsed_tbi,
            ids = ids,
            min_depth = min_depth,
            max_depth = max_depth,
            min_alt_reads = min_alt_reads
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
        
        Array[String] ids
        Array[Int] min_depth
        Array[Int] max_depth
        Array[Int] min_alt_reads
    }
    parameter_meta {
    }
    
    Int n_cpus = 4*length(ids)
    Int ram_size_gb = 8*length(ids)
    Int disk_size_gb = (length(ids)+1)*10*ceil(size(truvari_collapsed_vcf_gz, "GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        function filter_thread() {
            local ID=$1
            local MIN_DEPTH=$2
            local MAX_DEPTH=$3
            local MIN_ALT_READS=$4
            
            echo "'GT=\"alt\" && (DP < ${MIN_DEPTH} || DP > ${MAX_DEPTH} || AD[*:1] < ${MIN_ALT_READS})'" > ${ID}_filter.txt
            while read FILTER; do
                ${TIME_COMMAND} bcftools filter --exclude ${FILTER} --set-GTs . --output-type z ~{truvari_collapsed_vcf_gz} > ${ID}_tmp1.vcf.gz
            done < ${ID}_filter.txt
            tabix -f ${ID}_tmp1.vcf.gz
            ${TIME_COMMAND} bcftools filter --threads 1 --include 'COUNT(GT="alt")>0' --output-type z ${ID}_tmp1.vcf.gz > ${ID}_tmp2.vcf.gz
            tabix -f ${ID}_tmp2.vcf.gz
            rm -f ${ID}_tmp1.vcf.gz*
            bcftools view --header-only ${ID}_tmp2.vcf.gz > ${ID}_header.txt
            N_ROWS=$(wc -l < ${ID}_header.txt)
            head -n $(( ${N_ROWS} - 1 )) ${ID}_header.txt > ${ID}_filtered.vcf
            echo '##INFO=<ID=ORIGINAL_ID,Number=1,Type=String,Description="Original ID from truvari collapse">' >> ${ID}_filtered.vcf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> ${ID}_filtered.vcf
            bcftools view --no-header ${ID}_tmp2.vcf.gz | awk 'BEGIN { i=0; } { gsub(/;/,"_",$3); printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s;ORIGINAL_ID=%s\tGT\t0/1\n",$1,$2,++i,$4,$5,$6,$7,$8,$3); }' >> ${ID}_filtered.vcf
            rm -f ${ID}_tmp2.vcf.gz*
            ${TIME_COMMAND} bgzip -@ 4 ${ID}_filtered.vcf
            ${TIME_COMMAND} tabix -f ${ID}_filtered.vcf.gz
        }
        
        
        # Main program
        echo ~{sep="," ids} | tr ',' '\n' > ids.txt
        echo ~{sep="," min_depth} | tr ',' '\n' > min_depth.txt
        echo ~{sep="," max_depth} | tr ',' '\n' > max_depth.txt
        echo ~{sep="," min_alt_reads} | tr ',' '\n' > min_alt_reads.txt
        paste -d , ids.txt min_depth.txt max_depth.txt min_alt_reads.txt > tasks.tsv
        while read ROW; do
            ID=$(echo ${ROW} | cut -d , -f 1)
            MIN_DEPTH=$(echo ${ROW} | cut -d , -f 2)
            MAX_DEPTH=$(echo ${ROW} | cut -d , -f 3)
            MIN_ALT_READS=$(echo ${ROW} | cut -d , -f 4)
            filter_thread ${ID} ${MIN_DEPTH} ${MAX_DEPTH} ${MIN_ALT_READS} &
        done < tasks.tsv
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
