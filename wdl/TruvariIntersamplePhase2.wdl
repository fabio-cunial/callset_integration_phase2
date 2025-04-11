version 1.0


# Runs a truvari inter-sample merge in parallel on every chromosome.
#
workflow TruvariIntersamplePhase2 {
    input {
        String source_dir
        String filter_string = "none"
        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        File samples_file
        
        File reference_fai
        File density_counter_py
        Int max_records_per_chunk = 10000
    }
    parameter_meta {
        source_dir: "Contains per-chromosome files built by `Split.wdl`."
        filter_string: "Apply this filter to every VCF before merging. 'none'=no filter."
        max_records_per_chunk: "Discards chunks that contain more than this many records. Setting it to 10k keeps 99.9% of all chunks in AoU Phase 1 (1027 samples) on CHM13."
    }

    scatter (chr in chromosomes) {
        call BcftoolsMerge {
            input:
                source_dir = source_dir,
                chromosome = chr,
                filter_string = filter_string
        }
        call GetBed {
            input:
                chromosome = chr,
                bcftools_merged_vcf_gz = BcftoolsMerge.bcftools_merged_vcf_gz,
                bcftools_merged_tbi = BcftoolsMerge.bcftools_merged_tbi,
                reference_fai = reference_fai,
                density_counter_py = density_counter_py,
                max_records_per_chunk = max_records_per_chunk
        }
        call Truvari {
            input:
                chromosome = chr,
                bcftools_merged_vcf_gz = BcftoolsMerge.bcftools_merged_vcf_gz,
                bcftools_merged_tbi = BcftoolsMerge.bcftools_merged_tbi,
                include_bed = GetBed.include_bed
        }
    }
    call ConcatenateChromosomes {
        input:
            chromosomes_vcf_gz = Truvari.collapsed_vcf_gz,
            chromosomes_tbi = Truvari.collapsed_tbi,
            samples_file = samples_file
    }
    
    output {
        File vcf_gz = ConcatenateChromosomes.vcf_gz
        File vcf_gz_tbi = ConcatenateChromosomes.vcf_gz_tbi
    }
}


# Resource usage on 449 samples on chr1:
#
# SETTING     CPU       RAM     TIME
# hg38<=0.7   150%      4G      1h30m
# hg38<=0.9   160%      5G      2h
# chm13<=0.7  130%      2G      1h10m
# chm13<=0.9  140%      3G      1h40m
#
task BcftoolsMerge {
    input {
        String source_dir
        String chromosome
        String filter_string
        
        Int n_cpu = 2
        Int ram_size_gb = 8
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
        
        # Downloading
        while : ; do
            TEST=$(gsutil -m cp ~{source_dir}'/*_'~{chromosome}_split.vcf.gz'*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        find . -maxdepth 1 -name '*.vcf.gz' > list.txt
        
        # Filtering, if needed.
        FILTER_STRING="~{filter_string}"
        if [ ${FILTER_STRING} != none ]; then
            INCLUDE_STR="--include ${FILTER_STRING}"
            rm -f list_filtered.txt
            while read FILE; do
                ID=$(basename ${FILE} .vcf.gz)
                bcftools filter --threads ${N_THREADS} ${INCLUDE_STR} --output-type z ${FILE} > ${ID}_filtered.vcf.gz
                tabix -f ${ID}_filtered.vcf.gz
                echo ${ID}_filtered.vcf.gz >> list_filtered.txt
                rm -f ${FILE}*
            done < list.txt
            rm -f list.txt
            mv list_filtered.txt list.txt
        fi
        
        # BCFTOOLS
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type z > ~{chromosome}.merged.vcf.gz
        tabix -f ~{chromosome}.merged.vcf.gz
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type z ~{chromosome}.merged.vcf.gz > ~{chromosome}.normed.vcf.gz
        tabix -f ~{chromosome}.normed.vcf.gz
    >>>
    
    output {
        File bcftools_merged_vcf_gz = work_dir + "/" + chromosome + ".normed.vcf.gz"
        File bcftools_merged_tbi = work_dir + "/" + chromosome + ".normed.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}


# Discards chunks with too many records
#
task GetBed {
    input {
        String chromosome
        File bcftools_merged_vcf_gz
        File bcftools_merged_tbi
        Int max_records_per_chunk
        File density_counter_py
        File reference_fai
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
        
        python ~{density_counter_py} ~{bcftools_merged_vcf_gz} > ~{chromosome}.chunks.bed
        awk '$4 >= ~{max_records_per_chunk}' ~{chromosome}.chunks.bed > ~{chromosome}.excluded.bed
        bedtools complement -i ~{chromosome}.excluded.bed -g ~{reference_fai} > ~{chromosome}.included.bed
        bedtools sort -faidx ~{reference_fai} -i ~{chromosome}.included.bed > ~{chromosome}.included.sorted.bed
    >>>
    
    output {
        File include_bed = work_dir + "/" + chromosome + ".included.sorted.bed"
        File chunks_bed = work_dir + "/" + chromosome + ".chunks.bed"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_resolve"
        cpu: 1
        memory: "8GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}


# Resource usage on 449 samples on chr1:
#
# SETTING     CPU       RAM     TIME
# hg38<=0.9   
# chm13<=0.9  
#
task Truvari {
    input {
        String chromosome
        File bcftools_merged_vcf_gz
        File bcftools_merged_tbi
        File include_bed
        Int ram_size_gb = 64
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
        
        ${TIME_COMMAND} truvari collapse --input ~{bcftools_merged_vcf_gz} --sizemin 0 --sizemax 1000000 --keep common --bed ~{include_bed} --gt all --output tmp.vcf
        ${TIME_COMMAND} bcftools sort --max-mem $(( ~{ram_size_gb} - 4 ))G --output-type z tmp.vcf > ~{chromosome}.collapsed.vcf.gz
        tabix -f ~{chromosome}.collapsed.vcf.gz
    >>>
    
    output {
        File collapsed_vcf_gz = work_dir + "/" + chromosome + ".collapsed.vcf.gz"
        File collapsed_tbi = work_dir + "/" + chromosome + ".collapsed.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: 8
        memory: ram_size_gb + "GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}


# Resource usage on 449 samples on chr1:
#
# SETTING     CPU       RAM     TIME
# hg38<=0.7   270%      4G      4m
# chm13<=0.7  200%      2G      4m
#
task ConcatenateChromosomes {
    input {
        Array[File] chromosomes_vcf_gz
        Array[File] chromosomes_tbi
        File samples_file
        
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
        
        # Ensuring that samples have the same order in all chromosome files
        rm -f list_filtered.txt
        while read FILE; do
            ID=$(basename ${FILE} .vcf.gz)
            bcftools view --samples-file ~{samples_file} --output-type z ${FILE} > ${ID}_filtered.vcf.gz
            tabix -f ${ID}_filtered.vcf.gz
            echo ${ID}_filtered.vcf.gz >> list_filtered.txt
            rm -f ${FILE}
        done < list.txt
        
        # Concatenating
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --file-list list_filtered.txt --output-type z > concat.vcf.gz
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
