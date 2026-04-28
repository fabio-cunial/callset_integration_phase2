version 1.0


# 
#
workflow SvimAsm {
    input {
        String sample_id
        
        File hap1_bam
        File hap1_bai
        File hap2_bam
        File hap2_bai

        File reference_fa
        File reference_fai
        
        String remote_output_dir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            sample_id = sample_id,
            hap1_bam = hap1_bam,
            hap1_bai = hap1_bai,
            hap2_bam = hap2_bam,
            hap2_bai = hap2_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            remote_output_dir = remote_output_dir,
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on a 4-core 8GB VM:
# TOOL                                  CPU%        RAM         TIME
# svim-asm                              100%        1.5GB       2m                       
#
task Impl {
    input {
        String sample_id
        
        File hap1_bam
        File hap1_bai
        File hap2_bam
        File hap2_bai

        File reference_fa
        File reference_fai
        
        String remote_output_dir
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 32
        Int disk_size_gb = 20
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        
        svim-asm --version 2>&1 || echo "1"

        ${TIME_COMMAND} svim-asm diploid ./svim/ ~{hap1_bam} ~{hap2_bam} ~{reference_fa}
        ${TIME_COMMAND} bcftools view --output-type z ./svim/variants.vcf --output ./~{sample_id}_svim.vcf.gz
        bcftools index --threads ${N_THREADS} -f -t ~{sample_id}_svim.vcf.gz

        gcloud storage cp ~{sample_id}_svim.vcf.'gz*' ~{remote_output_dir}/


        # # Splitting multiallelic records into biallelic records
        # ${TIME_COMMAND} bcftools norm --multiallelics -any --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
        # rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz

        # # Removing SNVs, replacement records, records that are not marked
        # # as ALT, records with a FILTER, and records with unresolved
        # # REF/ALT.
        # ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || (STRLEN(REF)>1 && STRLEN(ALT)>1) || GT!="alt" || (FILTER!="PASS" && FILTER!=".") || REF="*" || ALT="*"' --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
        # rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz

        # # Making sure SVLEN and SVTYPE are consistently annotated        
        # ${TIME_COMMAND} java -cp ~{docker_dir} AddSvtypeSvlen ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
        # rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz

        # # Keeping only INS and DEL in the given length range
        # ${TIME_COMMAND} bcftools filter --include '(SVTYPE="INS" || SVTYPE="DEL") && ABS(SVLEN)>='${MIN_SV_LENGTH} --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
        # rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz

        # mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_canonized.vcf.gz
        # mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_canonized.vcf.gz.tbi

    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
