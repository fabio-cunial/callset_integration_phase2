version 1.0


# 
#
workflow SV_Integration_Workpackage5_Ultralong {
    input {
        File sample_ids
        
        String remote_indir
        String remote_outdir
    }
    parameter_meta {
        sample_ids: "Speficies the order of the samples used by bcftools merge."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            sample_ids = sample_ids,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}



#
task BcftoolsMerge {
    input {
        File sample_ids
        
        String remote_indir
        
        Int n_cpu = 4
        Int ram_size_gb = 36
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        # Localizing all samples
        N_SAMPLES=$(cat ~{sample_ids} | wc -l)
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir}/'*'_ultralong.vcf.'gz*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading some files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        N_DOWNLOADED_SAMPLES=$(ls *_ultralong.vcf.gz | wc -l)
        if [ ${N_DOWNLOADED_SAMPLES} -ne ${N_SAMPLES} ]; then
            echo "Error: the number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is different from the number of samples specified (${N_SAMPLES})."
            exit 1
        fi
        
        # Merging
        while read SAMPLE_ID; do
            echo ${SAMPLE_ID}_ultralong.vcf.gz >> list.txt
        done < ~{sample_ids}
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type z > merged.vcf.gz
        tabix -f merged.vcf.gz
        ls -laht merged.vcf.gz
        df -h
        while read FILE; do
            rm -f ${FILE}
        done < list.txt
        
        # Splitting multiallelic records into biallelic records
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics - --output-type z merged.vcf.gz > normed.vcf.gz
        tabix -f normed.vcf.gz
    >>>
    
    output {
        File merged_normed_vcf_gz = "normed.vcf.gz"
        File merged_normed_tbi = "normed.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}



#
task TruvariCollapse {
    input {
        File merged_vcf_gz
        File merged_tbi
        
        String remote_outdir
        
        Int n_cpu = 4
        Int ram_size_gb = 36
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"

        
        # Collapsing
-----------> Should follow the corresponding workpackage...
----------> Also remember that these have all been compressed to symbolic by Workpackage1, to save space.....
        
        
        
        while read SAMPLE_ID; do
            echo ${SAMPLE_ID}_ultralong.vcf.gz >> list.txt
        done < ~{sample_ids}
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type z > merged.vcf.gz
        tabix -f merged.vcf.gz
        ls -laht merged.vcf.gz
        df -h
        while read FILE; do
            rm -f ${FILE}
        done < list.txt
        
        # Splitting multiallelic records into biallelic records
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics - --output-type z merged.vcf.gz > normed.vcf.gz
        tabix -f normed.vcf.gz
    >>>
    
    output {
        File merged_normed_vcf_gz = "normed.vcf.gz"
        File merged_normed_tbi = "normed.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
