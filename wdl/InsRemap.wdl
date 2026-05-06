version 1.0

#
workflow InsRemap {
    input {
        File vcf_gz
        File vcf_tbi

        File ref_fa
        File ref_fai

        Int max_length = 1000000
        Float cov_threshold = 0.8

        String docker_image = "quay.io/ymostovoy/lr-remap:latest"
    }

    call Impl {
        input:
            vcf_gz = vcf_gz,
            vcf_tbi = vcf_tbi,

            ref_fa = ref_fa,
            ref_fai = ref_fai,

            max_length = max_length,
            cov_threshold = cov_threshold,

            docker_image = docker_image
    }
}


#
task Impl {
    input {
        File vcf_gz
        File vcf_tbi

        File ref_fa
        File ref_fai

        Int max_length
        Float cov_threshold

        String docker_image
        Int n_cpu = 8
        Int mem_gb = 32
        Int disk_size_gb = 50
    }

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        mkdir ./ref_files
        mv ~{ref_fa} ./ref_files/
        mv ~{ref_fai} ./ref_files/
        REF_FA_BASENAME=$(basename ~{ref_fa})
        time truvari anno remap --threads ${N_THREADS} --aligner minimap2 --min-length 1 --max-length ~{max_length} --cov-threshold ~{cov_threshold} -r ./ref_files/${REF_FA_BASENAME} ~{vcf_gz} -o remapped.vcf.gz     
        time bcftools index --threads ${N_THREADS} --tbi remapped.vcf.gz
        bcftools query --format '%INFO/SUPP_SNIFFLES,%INFO/SUPP_PBSV,%INFO/SUPP_PAV,%INFO/SVLEN,%INFO/remap_classification,%INFO/remap_perc\n' remapped.vcf.gz > matrix.csv
    >>>
    
    output {
        File out_vcf_gz = "remapped.vcf.gz"
        File out_tbi = "remapped.vcf.gz.tbi"
        File matrix_csv = "matrix.csv"
    }

    runtime {
        cpu: n_cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " +  disk_size_gb + " HDD"
        preemptible: 0
        docker: docker_image
    }
}