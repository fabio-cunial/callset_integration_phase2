version 1.0


# 
#
workflow TestBetaBinomial {
    input {
        Array[String] sample_ids
        
        File intersample_vcf_gz
        File intersample_tbi
        
        String betabinomial_params = ""
    }
    parameter_meta {
        sample_ids: "Only these samples in `intersample_vcf_gz` will be re-genotyped."
    }
    
    scatter (i in range(length(sample_ids))) {
        call BetaBinomial {
            input:
                sample_id = sample_ids[i],
                intersample_vcf_gz = intersample_vcf_gz,
                intersample_tbi = intersample_tbi,
                betabinomial_params = betabinomial_params
        }
    }
    call Merge {
        input:
            sample_vcf_gz = BetaBinomial.regenotyped_vcf_gz,
            sample_tbi = BetaBinomial.regenotyped_tbi
    }
    
    output {
        File merged_vcf_gz = Merge.merged_vcf_gz
        File merged_tbi = Merge.merged_tbi
    }
}


#
task BetaBinomial {
    input {
        String sample_id
        
        File intersample_vcf_gz
        File intersample_tbi
        
        String betabinomial_params
        
        Int n_cpu = 4
        Int ram_size_gb = 32
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(intersample_vcf_gz,"GB")) ) + 50
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        source activate pyro-kanpig
        ${TIME_COMMAND} bcftools view --samples ~{sample_id} --output-type z ~{intersample_vcf_gz} > kanpig.vcf.gz
        tabix -f kanpig.vcf.gz
        ${TIME_COMMAND} python3 ~{docker_dir}/genotype-beta-binomial-mixture.py --kanpig-vcf kanpig.vcf.gz --output-prefix out
        ls -laht
        df -h
        ${TIME_COMMAND} bgzip out.delta.tsv
        bcftools view --header-only kanpig.vcf.gz > annotations.vcf
        cat out.annot.tsv >> annotations.vcf
        rm -f out.annot.tsv
        ${TIME_COMMAND} bgzip annotations.vcf
        tabix -f annotations.vcf.gz
        ls -laht
        df -h
        ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --columns CHROM,POS,REF,ALT,FORMAT/GT,FORMAT/GQ,FORMAT/SQ --annotations annotations.vcf.gz --output-type z kanpig.vcf.gz > ~{sample_id}_kanpig_betabinomial.vcf.gz
        ls -laht
        df -h
        rm -f kanpig.vcf.gz*
        tabix -f ~{sample_id}_kanpig_betabinomial.vcf.gz
        mv annotations.vcf.gz ~{sample_id}_annotations.vcf.gz
        mv out.delta.tsv.gz ~{sample_id}_delta.tsv.gz
        ls -laht
        df -h
    >>>
    
    output {
        File regenotyped_vcf_gz = sample_id + "_kanpig_betabinomial.vcf.gz"
        File regenotyped_tbi = sample_id + "_kanpig_betabinomial.vcf.gz.tbi"
        File annotations_vcf_gz = sample_id + "_annotations.vcf.gz"
        File deltas_tsv_gz = sample_id + "_delta.tsv.gz"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_squish"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


#
task Merge {
    input {
        Array[File] sample_vcf_gz
        Array[File] sample_tbi
        
        Int n_cpu = 4
        Int ram_size_gb = 32
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(sample_vcf_gz,"GB")) ) + 50
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
        
        SAMPLES=~{sep="," sample_vcf_gz}
        SAMPLES=$(echo ${SAMPLES} | tr ',' ' ')
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --output-type z ${SAMPLES} > merged.vcf.gz
        ${TIME_COMMAND} tabix -f merged.vcf.gz
    >>>
    
    output {
        File merged_vcf_gz = "merged.vcf.gz"
        File merged_tbi = "merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_squish"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
