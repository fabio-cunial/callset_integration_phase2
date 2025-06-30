version 1.0


#
workflow Resolve {
    input {
        String sample_id
        File vcf_gz
        File vcf_gz_tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    call ResolveImpl {
        input:
            sample_id = sample_id,
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            reference_fa = reference_fa
    }
    
    output {
    	File resolved_vcf = ResolveImpl.resolved_vcf
    	File resolved_tbi = ResolveImpl.resolved_tbi
    }
}


#
task ResolveImpl {
    input {
        String sample_id
        File vcf_gz
        File vcf_gz_tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(vcf_gz,"GB")) ) + ceil(size(reference_fa,"GB")) + 50
    String docker_dir = "/callset_integration"
    String work_dir = "/mnt/disks/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Removing multiallelic records from the input
        bcftools norm --multiallelics - --output-type z ~{vcf_gz} > new.vcf.gz
        tabix -f new.vcf.gz
        
        # Step 1 - clean up the VCF
        # - Assigns quality scores to each SV caller's result
        #  - pav 4
        #  - pbsv 3
        #  - sniffles 2
        # - Resolves any symbolic alts (e.g. `<DEL>` with the sequence from the
        #   reference)
        #   - symbolic variants are given quality score of 1
        # - Filters out `<CNV>` from pbsv and `<INS>` from sniffles
        # - Filters out variants greater than 100Mbp
        # - Fills in blank genotypes with `0`
        # - Filters out BND variants
        # - Turn some deletion/insertion pairs into inversions
        # The quality scores are set based on which variant representations we
        # believe are generally more accurate with higher being better.
        python ~{docker_dir}/resolve.py new.vcf.gz ~{reference_fa} | bcftools norm --check-ref s --fasta-ref ~{reference_fa} -N -m-any | bcftools view -i "SVTYPE != 'BND'" -O z -o ~{sample_id}.resolved.vcf.gz
        tabix -f ~{sample_id}.resolved.vcf.gz
        
        # $inversion_guesser.py$ creates a DEL with wrong ALT when an INS
        # is followed by a DEL that cancels it out. I.e. this input:
        #
        # chr2    14066975        pbsv.INS.669    T       TATATATATGATATATATATCATATATATGATATATATGATATATATATCATATATATG
        # chr2    14066975        pbsv.DEL.670    TATATATATGATATATATATCATATATATGATATATATGATATATATATCATATATATG
        #
        # gives the following output:
        #
        # chr2    14066975        pbsv.DEL.670    TATATATATGATATATATATCATATATATGATATATATGATATATATATCATATATATG     TATATATATGATATATATATCATATATATGATATATATGATATATATATCATATATATG
        #python ~{docker_dir}/inversion_guesser.py -i $prename -o $outname
    >>>
    
    output {
    	File resolved_vcf = work_dir + "/" + sample_id + ".resolved.vcf.gz"
        File resolved_tbi = work_dir + "/" + sample_id + ".resolved.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_resolve"
        cpu: 1
        memory: "128GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
