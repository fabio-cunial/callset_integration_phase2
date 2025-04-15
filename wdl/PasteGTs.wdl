version 1.0


# Assume that we have built an inter-sample-merged VCF, that we have
# subsequently kept just one column of it (setting it to all 0/1), and that we
# have re-genotyped this single-column VCF using the BAM of every sample. The
# program pastes all the GT columns of such re-genotyped VCFs into an
# inter-sample VCF with updated genotypes.
#
# Remark: every re-genotyped VCF file is assumed to have exactly the same set
# of calls in the same order. This is usually the case if the files are the
# result of re-genotyping the same inter-sample VCF with different BAMs.
# If this is not the case the program gives a wrong output, and it should be
# replaced with $bcftools merge --merge none$.
#
workflow PasteGTs {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        Array[File] sample_gts
        
        Int n_cpus
    }
    parameter_meta {
    }
    
    call Paste {
        input:
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi,
            sample_gts = sample_gts,
            n_cpus = n_cpus
    }
    
    output {
        File output_vcf_gz = Paste.output_vcf_gz
        File output_tbi = Paste.output_tbi
    }
}


task Paste {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        Array[File] sample_gts
        
        Int n_cpus
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/callset_integration"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        function pasteThread() {
            local THREAD_ID=$1
            
            OUTPUT_FILE="columns_${THREAD_ID}.txt"; touch ${OUTPUT_FILE}
            FIELDS_FILE="fields_${THREAD_ID}.txt"; touch ${FIELDS_FILE}
            TMP_PREFIX="tmp_${THREAD_ID}"
            i="0";
            while read FILE; do
                i=$(( $i + 1 ))
                # Adding the new sample to the set of columns                
                head -n 1 ${FILE} > sample_${THREAD_ID}.txt
                if [ $i = "1" ]; then
                    mv sample_${THREAD_ID}.txt ${FIELDS_FILE}
                else
                    paste ${FIELDS_FILE} sample_${THREAD_ID}.txt > ${FIELDS_FILE}.prime
                    mv ${FIELDS_FILE}.prime ${FIELDS_FILE}
                fi
                echo "Current fields of thread ${THREAD_ID}:"; cat ${FIELDS_FILE}
                # Adding the new column to the body
                tail -n +2 ${FILE} > ${TMP_PREFIX}.txt
                if [ $i = "1" ]; then
                    mv ${TMP_PREFIX}.txt ${OUTPUT_FILE}
                else
                    ${TIME_COMMAND} paste ${OUTPUT_FILE} ${TMP_PREFIX}.txt > ${OUTPUT_FILE}.prime
                    mv ${OUTPUT_FILE}.prime ${OUTPUT_FILE}
                fi
            done < list_${THREAD_ID}
        }
        
        
        # Main program
        
        # Initializing the output VCF with the input VCF
        bcftools view --header-only ~{intersample_vcf_gz} > tmp.txt
        N_ROWS=$(wc -l < tmp.txt)
        head -n $(( ${N_ROWS} - 1 )) ~{intersample_vcf_gz} > header.txt
        tail -n 1 tmp.txt | cut -f 1,2,3,4,5,6,7,8,9 > fields.txt
        bcftools view --no-header ~{intersample_vcf_gz} | cut -f 1,2,3,4,5,6,7,8,9 > calls.txt
        rm -f ~{intersample_vcf_gz}
        
        # Pasting GT files in parallel
        INPUT_FILES=~{sep=',' sample_gts}
        echo ${INPUT_FILES} | tr ',' '\n' > files_list.txt
        N_ROWS=$(wc -l < files_list.txt)
        N_ROWS=$(( ${N_ROWS} / ${N_THREADS} ))
        if [ ${N_ROWS} -eq 0 ]; then
            N_ROWS=1
        fi 
        split -d -l ${N_ROWS} files_list.txt list_
        COLUMNS_FILES=""; FIELDS_FILES=""
        for LIST_FILE in $(find . -maxdepth 1 -name 'list_*' | sort); do
            ID=${LIST_FILE#./list_}
            pasteThread ${ID} &
            COLUMNS_FILES="${COLUMNS_FILES} columns_${ID}.txt"
            FIELDS_FILES="${FIELDS_FILES} fields_${ID}.txt"
        done
        wait
        paste fields.txt ${FIELDS_FILES} > fields_all.txt
        rm -f fields_*.txt
        ${TIME_COMMAND} paste calls.txt ${COLUMNS_FILES} > body.txt
        rm -f columns_*.txt
        cat header.txt fields_all.txt body.txt > merged.vcf
        rm -f header.txt fields_all.txt body.txt
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} merged.vcf
        tabix -f merged.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpus
        disks: "local-disk 128 HDD"
        preemptible: 0
    }
}