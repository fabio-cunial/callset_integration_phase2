version 1.0


# Given a single-sample VCF re-genotyped with kanpig, the program recomputes its
# genotypes with beta binomial.
#
workflow TestBetaBinomial2 {
    input {
        String sample_id
        
        File kanpig_vcf_gz
        File kanpig_tbi
        
        String betabinomial_params = ""
        String remote_output_dir
    }
    parameter_meta {
    }
    
    call BetaBinomial {
        input:
            sample_id = sample_id,
            kanpig_vcf_gz = kanpig_vcf_gz,
            kanpig_tbi = kanpig_tbi,
            betabinomial_params = betabinomial_params,
            remote_output_dir = remote_output_dir
    }
    
    output {
    }
}


#
task BetaBinomial {
    input {
        String sample_id
        
        File kanpig_vcf_gz
        File kanpig_tbi
        
        String betabinomial_params
        String remote_output_dir
        
        Int n_cpu = 4
        Int ram_size_gb = 32
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(kanpig_vcf_gz,"GB")) ) + 50
    String docker_dir = "/callset_integration"
    
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
        
        # Computing new genotypes
        source activate pyro-kanpig
        ${TIME_COMMAND} python3 ~{docker_dir}/genotype-beta-binomial-mixture.py --kanpig-vcf ~{kanpig_vcf_gz} --output-prefix out
        
        # Annotating the input VCF
        bcftools view --header-only ~{kanpig_vcf_gz} > annotations.vcf
        cat out.annot.tsv >> annotations.vcf
        rm -f out.annot.tsv
        ${TIME_COMMAND} bgzip annotations.vcf
        tabix -f annotations.vcf.gz
        ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --columns CHROM,POS,REF,ALT,FORMAT/GT,FORMAT/GQ,FORMAT/SQ --annotations annotations.vcf.gz --output-type z ~{kanpig_vcf_gz} > ~{sample_id}_kanpig_betabinomial.vcf.gz
        tabix -f ~{sample_id}_kanpig_betabinomial.vcf.gz        
        
        # Outputting
        mv annotations.vcf.gz ~{sample_id}_annotations.vcf.gz
        ${TIME_COMMAND} bgzip out.delta.tsv
        mv out.delta.tsv.gz ~{sample_id}_delta.tsv.gz
        for FILE in $(ls *.png); do
            mv ${FILE} ~{sample_id}_${FILE}
        done
        while : ; do
            TEST=$(gsutil -m cp ~{sample_id}'_*' ~{remote_output_dir} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_squish"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
