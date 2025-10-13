version 1.0


# Given a VCF, counts all calls in TRs and of given lengths.
#
workflow GetNCalls2 {
    input {        
        String vcf_id
        File vcf_gz
        File tbi
        
        Array[Int] svlen_from
        Array[Int] svlen_to
        
        File tandem_bed
        File reference_fai
    }
    parameter_meta {
    }
    
    scatter (i in range(length(svlen_from))) {
        call GetNCalls {
            input:
                vcf_id = vcf_id,
                vcf_gz = vcf_gz,
                tbi = tbi,
                
                tandem_bed = tandem_bed,
                svlen_from = svlen_from[i],
                svlen_to = svlen_to[i]
        }
    }
    
    output {
        Array[Array[File]] out_txt = GetNCalls.out_txt
    }
}


#
task GetNCalls {
    input {
        String vcf_id
        File vcf_gz
        File tbi
        
        File tandem_bed
        Int svlen_from
        Int svlen_to
        
        Int n_cpu = 2
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 20*( ceil(size(vcf_gz,"GB")) ) + 100
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Making sure every call has the SVLEN and SVTYPE fields
        source activate truvari5
        date
        truvari anno svinfo -m 1 ~{vcf_gz} | bgzip > input.vcf.gz
        date
        tabix -f input.vcf.gz
        
        # Counting
        ${TIME_COMMAND} bcftools query --targets-file ~{tandem_bed} --targets-overlap pos --include 'SVTYPE=="INS" && ( (SVLEN>='~{svlen_from}' && SVLEN<'~{svlen_to}') || (SVLEN>-'~{svlen_to}' && SVLEN<=-'~{svlen_from}') )' --format '%CHROM' input.vcf.gz | wc -l > ~{vcf_id}_tr_ins_~{svlen_to}.txt &
        ${TIME_COMMAND} bcftools query --targets-file ~{tandem_bed} --targets-overlap pos --include 'SVTYPE=="DEL" && ( (SVLEN>='~{svlen_from}' && SVLEN<'~{svlen_to}') || (SVLEN>-'~{svlen_to}' && SVLEN<=-'~{svlen_from}') )' --format '%CHROM' input.vcf.gz | wc -l > ~{vcf_id}_tr_del_~{svlen_to}.txt &
        wait
    >>>
    
    output {
        Array[File] out_txt = glob("*.txt")
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 2
    }
}
