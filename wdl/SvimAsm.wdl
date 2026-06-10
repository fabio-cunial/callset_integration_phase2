version 1.0


# 
#
workflow SvimAsm {
    input {
        String sample_id
        Int max_sv_length = 10000000
        
        File hap1_bam
        File hap1_bai
        File hap2_bam
        File hap2_bai

        File reference_fa
        File reference_fai
        
        String remote_output_dir
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:latest"
        Int preemptible_number = 0
    }
    parameter_meta {
        max_sv_length: "SVIM-asm's default is 100kb."
    }
    
    call Impl {
        input:
            sample_id = sample_id,
            max_sv_length = max_sv_length,
            hap1_bam = hap1_bam,
            hap1_bai = hap1_bai,
            hap2_bam = hap2_bam,
            hap2_bai = hap2_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            remote_output_dir = remote_output_dir,
            docker_image = docker_image,
            preemptible_number = preemptible_number
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
        Int max_sv_length
        
        File hap1_bam
        File hap1_bai
        File hap2_bam
        File hap2_bai

        File reference_fa
        File reference_fai
        
        String remote_output_dir
        
        String docker_image
        Int n_cpu = 1
        Int ram_size_gb = 3
        Int disk_size_gb = 20
        Int preemptible_number
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        svim-asm --version 2>&1 || echo "1"

        # Capturing all SVs, and the destination of interspersed duplications.
        ${TIME_COMMAND} svim-asm diploid --max_sv_size ~{max_sv_length} --interspersed_duplications_as_insertions ./svim/ ~{hap1_bam} ~{hap2_bam} ~{reference_fa}
        bgzip -c ./svim/variants.vcf > ~{sample_id}_in1.vcf.gz
        bcftools index -f -t ~{sample_id}_in1.vcf.gz
        rm -rf ./svim/

        # Capturing only the source of interspersed duplications
        ${TIME_COMMAND} svim-asm diploid --max_sv_size ~{max_sv_length} ./svim/ ~{hap1_bam} ~{hap2_bam} ~{reference_fa}
        bcftools filter --include 'SVTYPE="DUP:INT"' --output-type z ./svim/variants.vcf --output ~{sample_id}_in2.vcf.gz
        bcftools index -f -t ~{sample_id}_in2.vcf.gz
        rm -rf ./svim/
        
        # Concatenating
        ${TIME_COMMAND} bcftools concat --allow-overlaps --remove-duplicates --output-type v ~{sample_id}_in1.vcf.gz ~{sample_id}_in2.vcf.gz --output ~{sample_id}_in.vcf
        rm -f ~{sample_id}_in1.vcf.gz* ~{sample_id}_in2.vcf.gz*

        # Splitting multiallelic records into biallelic records, if any.
        ${TIME_COMMAND} bcftools norm --multiallelics -any --output-type v ~{sample_id}_in.vcf --output ~{sample_id}_out.vcf
        rm -f ~{sample_id}_in.vcf ; mv ~{sample_id}_out.vcf ~{sample_id}_in.vcf

        # Removing records that are not marked as ALT and records with 
        # unresolved REF/ALT, if any.
        # Remark: we keep records with a FILTER, since the only filters are
        # incomplete_inversion ("Only one inversion breakpoint is supported")
        # and not_fully_covered ("Tandem duplication is not fully covered by a 
        # contig").
        ${TIME_COMMAND} bcftools filter --exclude 'GT!="alt" || REF="*" || ALT="*"' --output-type v ~{sample_id}_in.vcf --output ~{sample_id}_out.vcf
        rm -f ~{sample_id}_in.vcf ; mv ~{sample_id}_out.vcf ~{sample_id}_in.vcf

        # Making sure SVLEN and SVTYPE are consistently annotated        
        ${TIME_COMMAND} java -cp ~{docker_dir} AddSvtypeSvlen ~{sample_id}_in.vcf > ~{sample_id}_out.vcf
        rm -f ~{sample_id}_in.vcf ; mv ~{sample_id}_out.vcf ~{sample_id}_in.vcf

        # Uploading
        bgzip -c ~{sample_id}_in.vcf > ~{sample_id}_canonized.vcf.gz
        bcftools index -f -t ~{sample_id}_canonized.vcf.gz
        gcloud storage cp ~{sample_id}_canonized.vcf.'gz*' ~{remote_output_dir}/
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
    }
}
