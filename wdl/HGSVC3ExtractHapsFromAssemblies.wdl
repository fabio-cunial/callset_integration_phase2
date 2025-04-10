version 1.0


# 
#
workflow HGSVC3ExtractHapsFromAssemblies {
    input {
        String input_dir
        String output_dir
    }
    parameter_meta {
    }
    
    call ExtractHapsImpl {
        input:
            input_dir = input_dir,
            output_dir = output_dir
    }
    
    output {
    }
}


#
task ExtractHapsImpl {
    input {
        String input_dir
        String output_dir
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
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        SEQKIT_COMMAND="./seqkit"

        gsutil ls ~{input_dir} > assemblies.txt
        while read FILE; do
            LOCAL_FILE=$(basename ${FILE})
            gsutil -m cp ${FILE} .
            zcat ${LOCAL_FILE} | grep '^>' > headers.txt
            grep '.h1tg' headers.txt > tmp.txt
            cat tmp.txt
            cut -d ' ' -f 1 tmp.txt | tr -d '>' > contigs1.txt
            rm -f tmp.txt
            grep '.h2tg' headers.txt > tmp.txt
            cat tmp.txt
            cut -d ' ' -f 1 tmp.txt | tr -d '>' > contigs2.txt
            rm -f tmp.txt
            rm -f headers.txt
            ID=${LOCAL_FILE%.fna.gz}
            ID=${LOCAL_FILE#*_}
            ID=${ID#*_}
            ID=${ID%%.*}
            ~{docker_dir}/seqkit grep --pattern-file contigs1.txt ${LOCAL_FILE} | bgzip --compress-level 1 > ${ID}_hap1.fna.gz &
            ~{docker_dir}/seqkit grep --pattern-file contigs2.txt ${LOCAL_FILE} | bgzip --compress-level 2 > ${ID}_hap1.fna.gz &
            wait
            gsutil -m mv ${ID}_hap1.fna.gz ${OUTPUT_DIR}
            gsutil -m mv ${ID}_hap2.fna.gz ${OUTPUT_DIR}
            rm -f ${LOCAL_FILE}.fna.gz contigs*.txt
        done < assemblies.txt
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: 1
        memory: "16GB"
        disks: "local-disk 500 HDD"
        preemptible: 0
    }
}
