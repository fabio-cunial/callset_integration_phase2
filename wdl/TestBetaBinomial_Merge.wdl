version 1.0


# Just the merge part of `TestBetaBinomial.wdl`.
#
workflow TestBetaBinomial_Merge {
    input {
        String remote_output_dir
    }
    parameter_meta {
    }
    
    call Merge {
        input:
            remote_output_dir = remote_output_dir
    }
    
    output {
    }
}



task Merge {
    input {
        String remote_output_dir
        
        Int n_cpu = 4
        Int ram_size_gb = 32
        Int disk_size_gb = 200
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Merging
        gsutil -m cp ~{remote_output_dir}/'*_kanpig_betabinomial.vcf.gz*' .
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --output-type z *_kanpig_betabinomial.vcf.gz > merged.vcf.gz
        ${TIME_COMMAND} tabix -f merged.vcf.gz
        
        # Outputting
        while : ; do
            TEST=$(gsutil -m cp merged.vcf.'gz*' ~{remote_output_dir} && echo 0 || echo 1)
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
        docker: "fcunial/callset_integration_phase2_squish"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}