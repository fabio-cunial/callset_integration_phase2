version 1.0


# 
#
workflow UltralongMerge {
    input {
        String remote_indir
        String remote_outdir
        
        String svtype
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong"
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            svtype = svtype,
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 4-core, 4GB VM:
#
# TOOL                                                CPU     RAM     TIME
# bcftools merge
# bcftools norm
#
task Impl {
    input {
        String remote_indir
        String remote_outdir
        
        String svtype
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 4
        Int disk_size_gb = 50
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Merging
        df -h 1>&2
        ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/'*_'~{svtype}'.vcf.gz*' ~{remote_indir}/'*_'~{svtype}'_training.vcf.gz*' .
        df -h 1>&2
        ls *_~{svtype}.vcf.gz > list1.txt
        ls *_~{svtype}_training.vcf.gz > list2.txt
        ${TIME_COMMAND} bcftools merge --threads 2 --force-samples --merge none --file-list list1.txt --output-type z --output ~{svtype}_merged.vcf.gz &
        ${TIME_COMMAND} bcftools merge --threads 2 --force-samples --merge none --file-list list2.txt --output-type z --output ~{svtype}_training_merged.vcf.gz &
        wait
        ${TIME_COMMAND} bcftools index --threads 2 -f ~{svtype}_merged.vcf.gz &
        ${TIME_COMMAND} bcftools index --threads 2 -f ~{svtype}_training_merged.vcf.gz &
        wait
        ls -laht 1>&2
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type z ~{svtype}_merged.vcf.gz --output ~{svtype}_normed.vcf.gz
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type z ~{svtype}_training_merged.vcf.gz --output ~{svtype}_training_normed.vcf.gz
        ${TIME_COMMAND} bcftools index --threads 2 -f ~{svtype}_normed.vcf.gz &
        ${TIME_COMMAND} bcftools index --threads 2 -f ~{svtype}_training_normed.vcf.gz &
        wait
        ls -laht 1>&2
        
        # Uploading
        ${TIME_COMMAND} gcloud storage mv ~{svtype}_normed.vcf.gz ~{remote_outdir}/~{svtype}_merged.vcf.gz
        ${TIME_COMMAND} gcloud storage mv ~{svtype}_normed.vcf.gz.tbi ~{remote_outdir}/~{svtype}_merged.vcf.gz.tbi
        ${TIME_COMMAND} gcloud storage mv ~{svtype}_training_normed.vcf.gz ~{remote_outdir}/~{svtype}_training_merged.vcf.gz
        ${TIME_COMMAND} gcloud storage mv ~{svtype}_training_normed.vcf.gz.tbi ~{remote_outdir}/~{svtype}_training_merged.vcf.gz.tbi
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
