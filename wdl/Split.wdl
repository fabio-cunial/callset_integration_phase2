version 1.0


# Splits a VCF into a given set of chromosomes.
#
workflow Split {
    input {
        String sample_id
        File sample_vcf_gz
        File sample_vcf_gz_tbi
        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        String destination_dir
    }
    parameter_meta {
        sample_vcf_gz: "Assumed to be already sorted"
        destination_dir: "The filtered and split files are stored in this remote directory."
    }

    call SplitImpl {
        input:
            sample_id = sample_id,
            sample_vcf_gz = sample_vcf_gz,
            sample_vcf_gz_tbi = sample_vcf_gz_tbi,
            chromosomes = chromosomes,
            destination_dir = destination_dir
    }
    
    output {
    }
}


task SplitImpl {
    input {
        String sample_id
        File sample_vcf_gz
        File sample_vcf_gz_tbi
        Array[String] chromosomes
        String destination_dir
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 10*(ceil(size(sample_vcf_gz,"GB")))  # Arbitrary

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Transfering annotations to the FORMAT field, so that they are
        # preserved by the inter-sample merge later.
        echo '##FORMAT=<ID=CALIBRATION_SENSITIVITY,Number=1,Type=Float,Description="Calibration sensitivity according to the model applied by ScoreVariantAnnotations">' > header.txt
        echo '##FORMAT=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> header.txt
        echo '##FORMAT=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by sniffles">' >> header.txt
        echo '##FORMAT=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by pav">' >> header.txt
        bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CALIBRATION_SENSITIVITY\t%SUPP_PBSV\t%SUPP_SNIFFLES\t%SUPP_PAV\n' ~{sample_vcf_gz} | bgzip -c > annotations.tsv.gz
        tabix -s1 -b2 -e2 annotations.tsv.gz
        ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations annotations.tsv.gz --header-lines header.txt --columns CHROM,POS,ID,REF,ALT,FORMAT/CALIBRATION_SENSITIVITY,FORMAT/SUPP_PBSV,FORMAT/SUPP_SNIFFLES,FORMAT/SUPP_PAV ~{sample_vcf_gz} --output-type z > formatted.vcf.gz
        tabix -f formatted.vcf.gz
        
        # Splitting
        CHROMOSOMES=~{sep='-' chromosomes}
        CHROMOSOMES=$(echo ${CHROMOSOMES} | tr '-' ' ')
        for CHROMOSOME in ${CHROMOSOMES}; do
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --output-type z formatted.vcf.gz ${CHROMOSOME} > ~{sample_id}_${CHROMOSOME}_split.vcf.gz
            tabix -f ~{sample_id}_${CHROMOSOME}_split.vcf.gz
        done
        
        # Uploading
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp '*_split.vcf.gz*' ~{destination_dir} && echo 0 || echo 1)
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
        docker: "fcunial/callset_integration_phase2"
        cpu: 1
        memory: "8GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
