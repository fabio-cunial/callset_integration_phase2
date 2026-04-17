version 1.0


# 
#
workflow UltralongMerge {
    input {
        String remote_indir
        String remote_outdir
        
        String svtype
        String suffix
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong"
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            svtype = svtype,
            suffix = suffix,
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 4-core, 4GB VM:
#
# TOOL                                                CPU     RAM     TIME
# bcftools concat
#
task Impl {
    input {
        String remote_indir
        String remote_outdir
        
        String svtype
        String suffix
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
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
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        cat << 'END' > fix_sample.sh
#!/bin/bash
INPUT_VCF_GZ=$1

bcftools reheader --samples-list SAMPLE ${INPUT_VCF_GZ} --output ${INPUT_VCF_GZ}.reheader
rm -f ${INPUT_VCF_GZ} ${INPUT_VCF_GZ}.tbi
mv ${INPUT_VCF_GZ}.reheader ${INPUT_VCF_GZ}
bcftools index -f -t ${INPUT_VCF_GZ}
END
        chmod +x fix_sample.sh
        
        
        # ---------------------------- Main program ----------------------------
        
        # Simple concatenation, with only exact duplicate removal. In the
        # future we may run truvari collapse to remove approximate duplicates.
        df -h 1>&2
        ls -laht 1>&2
        ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/'*_'~{svtype}~{suffix}'.vcf.gz*' .
        df -h 1>&2
        ls *.vcf.gz > list.txt
        ${TIME_COMMAND} xargs --arg-file=list.txt --max-lines=1 --max-procs=${N_THREADS} ./fix_sample.sh
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --remove-duplicates --file-list list.txt --output-type z --output ~{svtype}~{suffix}_merged.vcf.gz
        ${TIME_COMMAND} bcftools index --threads 2 -f -t ~{svtype}~{suffix}_merged.vcf.gz
        ls -laht 1>&2
        
        # Uploading
        ${TIME_COMMAND} gcloud storage mv ~{svtype}~{suffix}'*_merged.vcf.gz*' ~{remote_outdir}/
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
