version 1.0

#
workflow InsRemap {
    input {
        String sample_id
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
            sample_id = sample_id,
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
        String sample_id
        File vcf_gz
        File vcf_tbi

        File ref_fa
        File ref_fai

        Int max_length
        Float cov_threshold

        String docker_image
        Int n_cpu = 8
        Int mem_gb = 16
    }
   
    Int disk_size = ceil(size(vcf_gz, "GB") + size(ref_fa, "GB")) * 3 + 50

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        

        mkdir ./ref_files
        mv ~{ref_fa} ./ref_files/
        mv ~{ref_fai} ./ref_files/
        REF_FA_BASENAME=$(basename ~{ref_fa})
        time truvari anno remap --threads ${N_THREADS} --aligner minimap2 --min-length 1 --max-length ~{max_length} --cov-threshold ~{cov_threshold} -r ./ref_files/${REF_FA_BASENAME} ~{vcf_gz} -o ~{sample_id}_ins_remapped.vcf.gz     
        time bcftools index --threads ${N_THREADS} --tbi ~{sample_id}_ins_remapped.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = "~{sample_id}_ins_remapped.vcf.gz"
        File out_tbi = "~{sample_id}_ins_remapped.vcf.gz.tbi"
    }

    runtime {
        cpu: n_cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " +  disk_size + " HDD"
        preemptible: 0
        docker: docker_image
    }
}