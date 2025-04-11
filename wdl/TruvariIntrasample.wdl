# Truvari intra-merge for AoU SV
version 1.0


# Workflow for intra-sample merging for AoU.
#
workflow TruvariIntrasample {
    input {
        String sample_id
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File pav_vcf_gz
        File pav_vcf_gz_tbi
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    call TruvariIntrasampleImpl {
        input:
            sample_id = sample_id,
            pbsv_vcf_gz = pbsv_vcf_gz,
            pbsv_vcf_gz_tbi = pbsv_vcf_gz_tbi,
            sniffles_vcf_gz = sniffles_vcf_gz,
            sniffles_vcf_gz_tbi = sniffles_vcf_gz_tbi,
            pav_vcf_gz = pav_vcf_gz,
            pav_vcf_gz_tbi = pav_vcf_gz_tbi,
            ram_size_gb = ram_size_gb
    }
    
    output {
    	File truvari_collapsed = TruvariIntrasampleImpl.truvari_collapsed
    	File truvari_collapsed_idx = TruvariIntrasampleImpl.truvari_collapsed_idx
    	File bcftools_merged = TruvariIntrasampleImpl.bcftools_merged
    	File bcftools_merged_idx = TruvariIntrasampleImpl.bcftools_merged_idx
    }
}


# Other intermediate files created, but likely aren't useful for production are:
# - preprocessed/ directory for each caller's cleaned result
# - ~{sample_id}.removed.vcf.gz variant representations removed during
# collapsing
#
task TruvariIntrasampleImpl {
    input {
        String sample_id
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File pav_vcf_gz
        File pav_vcf_gz_tbi
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(pbsv_vcf_gz,"GB")) + ceil(size(sniffles_vcf_gz,"GB")) + ceil(size(pav_vcf_gz,"GB")) ) + 50
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 4 ))

        # Step 1 - Merging
        # Pastes the samples together in the order of the preferred genotypes.
        # That is to say, this creates a three sample VCF with sample columns
        # from pbsv, sniffles, pav_sv
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples -O z -o tmp.vcf.gz ~{pbsv_vcf_gz} ~{sniffles_vcf_gz} ~{pav_vcf_gz}
        tabix -f tmp.vcf.gz
        
        # Step 2 - Removing multiallelic records. We observed that they are
        # created in Step 1 sometimes, e.g.:
        #
        # 2024-03-19 20:52:28,548 [ERROR] Cannot compare multi-allelic records.
        # Please split
        # 2024-03-19 20:52:28,548 [ERROR] line
        # chr4	137168756	pbsv.INS.2751;chr4-137168757-DEL-52	A   ACGTATGTGTATACGTATACATATACGCGTATATACATACGTATACATATACG,A	4	PASS	SVTYPE=INS;SVLEN=52;SVANN=TANDEM;ID=chr4-137168757-DEL-52;TIG_REGION=h2tg007223l:91206-91206;QUERY_STRAND=+;HOM_REF=0,23;HOM_TIG=0,23;INVScore=0.981132;AC=2,1	GT:AD:DP:SAC	1/1:1,8,.:9:1,0,3,5	./.:.:.:.	0|2:.:.:.
        #
        ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z tmp.vcf.gz > ~{sample_id}.bcftools_merged.vcf.gz
        tabix -f ~{sample_id}.bcftools_merged.vcf.gz
        rm -f tmp.vcf.gz*

        # Step 3 - Collapsing
        ${TIME_COMMAND} truvari collapse -i ~{sample_id}.bcftools_merged.vcf.gz -c removed.vcf.gz --sizemin 0 --sizemax 1000000 --keep maxqual --gt het --intra --pctseq 0.90 --pctsize 0.90 --refdist 500 --output tmp.vcf
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z tmp.vcf > ~{sample_id}.truvari_collapsed.vcf.gz
        tabix -f ~{sample_id}.truvari_collapsed.vcf.gz
    >>>
    
    output {
    	File truvari_collapsed = work_dir + "/" + sample_id + ".truvari_collapsed.vcf.gz"
    	File truvari_collapsed_idx = work_dir + "/" + sample_id + ".truvari_collapsed.vcf.gz.tbi"
    	File bcftools_merged = work_dir + "/" + sample_id + ".bcftools_merged.vcf.gz"
    	File bcftools_merged_idx = work_dir + "/" + sample_id + ".bcftools_merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: 1
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
