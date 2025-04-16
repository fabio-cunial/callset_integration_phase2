version 1.0


# Assume that we have built an inter-sample-merged VCF, that we have
# subsequently kept just one column of it (setting it to all 0/1), and that we
# have re-genotyped this single-column VCF using the BAM of every sample. The
# program pastes all the GT columns of such re-genotyped VCFs into an
# inter-sample VCF with updated genotypes.
#
workflow PasteGTs {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        Array[File] sample_gts
        File samples_file
    }
    parameter_meta {
        sample_gts: "The output of `KanpigMerged.wdl`. Remark: all such files have the same number of records; for every file, records appear in the order induced by a global sort of `intersample_vcf_gz`."
        intersample_vcf_gz: "The output of `TruvariIntersamplePhase2.wdl`. Remark: since this file may not have been globally sorted, the order of the records in this file might be different from the order of the records in each file of `sample_gts`."
        samples_file: "Order in which to store the samples in the output."
    }
    
    call Sort {
        input:
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi,
            samples_file = samples_file
    }
    call Paste {
        input:
            intersample_vcf_gz = Sort.output_vcf_gz,
            intersample_tbi = Sort.output_tbi,
            sample_gts = sample_gts,
            samples_file = samples_file
    }
    
    output {
        File output_vcf_gz = Paste.output_vcf_gz
        File output_tbi = Paste.output_tbi
    }
}


# Sorts `intersample_vcf_gz` and ensures that samples appear in the given order.
#
task Sort {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        File samples_file
        
        Int n_cpus = 4
        Int ram_size_gb = 64
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
        EFFECTIVE_MEM_GB=~{ram_size_gb}
        EFFECTIVE_MEM_GB=$(( ${EFFECTIVE_MEM_GB} - 4 ))
        
        ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z ~{intersample_vcf_gz} > intersample.vcf.gz
        tabix -f intersample.vcf.gz
        rm -f ~{intersample_vcf_gz}
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file ~{samples_file} --output-type z intersample.vcf.gz > output.vcf.gz
        tabix -f output.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/output.vcf.gz"
        File output_tbi = work_dir + "/output.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk 128 HDD"
        preemptible: 0
    }
}


task Paste {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        Array[File] sample_gts
        File samples_file
        
        Int n_cpus = 64
        Int ram_size_gb = 128
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
        EFFECTIVE_MEM_GB=~{ram_size_gb}
        EFFECTIVE_MEM_GB=$(( ${EFFECTIVE_MEM_GB} - 4 ))
        
        function pasteThread() {
            local THREAD_ID=$1
            
            OUTPUT_FILE="columns_${THREAD_ID}.txt"; touch ${OUTPUT_FILE}
            FIELDS_FILE="fields_${THREAD_ID}.txt"; touch ${FIELDS_FILE}
            LIST_FIELDS=""; LIST_BODIES=""; i="0";
            while read FILE; do
                i=$(( ${i} + 1 ))
                head -n 1 ${FILE} > ${THREAD_ID}_s_${i}.txt
                LIST_FIELDS="${LIST_FIELDS} ${THREAD_ID}_s_${i}.txt"
                tail -n +2 ${FILE} > ${THREAD_ID}_b_${i}.txt
                LIST_BODIES="${LIST_BODIES} ${THREAD_ID}_b_{i}.txt"
            done < list_${THREAD_ID}
            paste ${LIST_FIELDS} > ${FIELDS_FILE}
            ${TIME_COMMAND} paste ${LIST_BODIES} > ${OUTPUT_FILE}
            rm -f ${LIST_FIELDS}
            rm -f ${LIST_BODIES}
        }
        
        
        # Main program
        
        # - Initializing the output VCF with the input VCF.
        # Remark: this implies that call IDs in the output are identical to
        # truvari's, so they might not be all distinct.
        bcftools view --header-only ~{intersample_vcf_gz} > tmp.txt
        N_ROWS=$(wc -l < tmp.txt)
        head -n $(( ${N_ROWS} - 1 )) tmp.txt > header.txt
        tail -n 1 tmp.txt | cut -f 1,2,3,4,5,6,7,8,9 > fields.txt
        bcftools view --threads ${N_THREADS} --no-header ~{intersample_vcf_gz} | awk '{ printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tGT:FT:SQ:GQ:PS:NE:DP:AD:KS:SCORE:CALIBRATION_SENSITIVITY:SUPP_PBSV:SUPP_SNIFFLES:SUPP_PAV",$1,$2,$3,$4,$5,$6,$7,$8); }' > calls.txt
        
        # - Pasting the GT files in parallel (the order among the files/samples
        # is not important at this stage).
        INPUT_FILES=~{sep=',' sample_gts}
        echo ${INPUT_FILES} | tr ',' '\n' > files_list.txt
        N_ROWS=$(wc -l < files_list.txt)
        N_ROWS=$(( ${N_ROWS} / ${N_THREADS} ))
        if [ ${N_ROWS} -eq 0 ]; then
            N_ROWS=1
        fi 
        split -d -l ${N_ROWS} files_list.txt list_
        COLUMNS_FILES=""; FIELDS_FILES=""
        for LIST_FILE in $(find . -maxdepth 1 -name 'list_*'); do
            ID=${LIST_FILE#./list_}
            pasteThread ${ID} &
            COLUMNS_FILES="${COLUMNS_FILES} columns_${ID}.txt"
            FIELDS_FILES="${FIELDS_FILES} fields_${ID}.txt"
        done
        wait
        paste fields.txt ${FIELDS_FILES} > fields_all.txt
        ${TIME_COMMAND} paste calls.txt ${COLUMNS_FILES} > body.txt
        rm -f ${FIELDS_FILES} ${COLUMNS_FILES}
        
        cat fields_all.txt
        head -n 10 body.txt
        
        cat header.txt fields_all.txt body.txt > merged.vcf
        rm -f header.txt fields_all.txt body.txt
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} merged.vcf
        tabix -f merged.vcf.gz
        
        # - Ensuring that samples appear in the given order
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file ~{samples_file} --output-type z merged.vcf.gz > merged_prime.vcf.gz
        tabix -f merged_prime.vcf.gz
        rm -f merged.vcf.gz*
        
        # - Copying from `intersample_vcf_gz` the FORMAT fields that are
        # missing in `merged_prime.vcf.gz`.
        ${TIME_COMMAND} java -Xmx${EFFECTIVE_MEM_GB}G -cp ~{docker_dir} CopyFormat ~{intersample_vcf_gz} merged_prime.vcf.gz output.vcf.gz
        tabix -f output.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/output.vcf.gz"
        File output_tbi = work_dir + "/output.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk 128 HDD"
        preemptible: 0
    }
}