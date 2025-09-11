version 1.0


# 
#
workflow TestBetaBinomial2 {
    input {
        String sample_id
        
        File kanpig_vcf_gz
        File kanpig_tbi
        
        String betabinomial_params = ""
    }
    parameter_meta {
    }
    
    call BetaBinomial {
        input:
            sample_id = sample_id,
            kanpig_vcf_gz = kanpig_vcf_gz,
            kanpig_tbi = kanpig_tbi,
            betabinomial_params = betabinomial_params
    }
    
    output {
        File out_vcf_gz = BetaBinomial.regenotyped_vcf_gz
        File out_tbi = BetaBinomial.regenotyped_tbi
    }
}


#
task BetaBinomial {
    input {
        String sample_id
        
        File kanpig_vcf_gz
        File kanpig_tbi
        
        String betabinomial_params
        
        Int n_cpu = 4
        Int ram_size_gb = 32
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(kanpig_vcf_gz,"GB")) ) + 50
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
        ${TIME_COMMAND} python3 ~{docker_dir}/genotype-beta-binomial-mixture.py --kanpig-vcf ~{kanpig_vcf_gz} --output-prefix out
        ls -laht
        df -h
        ${TIME_COMMAND} bgzip out.delta.tsv
        bcftools view --header-only ~{kanpig_vcf_gz} > annotations.vcf
        cat out.annot.tsv >> annotations.vcf
        rm -f out.annot.tsv
        ${TIME_COMMAND} bgzip annotations.vcf
        tabix -f annotations.vcf.gz
        ls -laht
        df -h
        ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --columns CHROM,POS,REF,ALT,FORMAT/GT,FORMAT/GQ,FORMAT/SQ --annotations annotations.vcf.gz --output-type z ~{kanpig_vcf_gz} > ~{sample_id}_kanpig_betabinomial.vcf.gz
        ls -laht
        df -h
        rm -f ~{kanpig_vcf_gz}
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
