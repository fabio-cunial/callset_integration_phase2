version 1.0


# Pastes all the GT column files coming from the inter-sample-re-genotyped VCF
# into a new inter-sample VCF that holds the updated genotypes.
#
# Remark: all GT column files have the same number of records, and records
# appear in the same order, which is identical to their order in
# `truvari_collapse_intersample_vcf_gz`.
#
workflow TestKanpigIntersample {
    input {
        File truvari_collapse_intersample_vcf_gz
        File truvari_collapse_intersample_tbi
        
        String remote_indir
        String remote_outdir
        
        File samples_file
        String format_field = "GT:FT:SQ:GQ:PS:NE:DP:AD:KS"
    }
    parameter_meta {
        truvari_collapse_intersample_vcf_gz: "The file given in input to inter-sample kanpig. Its FORMAT field contains just `GT`."
        remote_indir: "Contains GT column files, whose rows are in the same order as `truvari_collapse_intersample_vcf_gz`."
        samples_file: "Order in which to store the samples in the output VCF."
    }
    
    call PasteGTs {
        input:
            truvari_collapse_intersample_vcf_gz = truvari_collapse_intersample_vcf_gz,
            truvari_collapse_intersample_tbi = truvari_collapse_intersample_tbi,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            samples_file = samples_file,
            format_field = format_field
    }
    
    output {
    }
}


# Pastes the truvari intersample VCF and all the GT files.
#
# Remark: pasting is done sequentially, since it is fast enough in practice.
#
# Performance on 10'070 samples, 15x, GRCh38, stringent (_S) and lenient (_L):
#
# TOOL                      CPU_S   RAM_S   TIME_S  CPU_L   RAM_L   TIME_L
# bcftools view | awk                       3.5h                    6h
# paste                     100%    50M     2h      100%    50M     3h
# bgzip                     200%    20M     1h      200%    20M     2h
#
task PasteGTs {
    input {
        File truvari_collapse_intersample_vcf_gz
        File truvari_collapse_intersample_tbi
        
        String remote_indir
        String remote_outdir
        File samples_file
        String format_field
        
        Int n_cpu = 4
        Int ram_size_gb = 16
        Int disk_size_gb = 3000
    }
    parameter_meta {
        disk_size_gb: "3TB for stringent and 4TB for lenient suffice."
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
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        # Downloading and checking GT files
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir}/'*_sorted.txt' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        N_RECORDS=$(bcftools index --nrecords ~{truvari_collapse_intersample_tbi})
        i="0"; COLUMNS_FILES=""; FIELDS_FILES="";
        while read SAMPLE_ID; do
            N=$(wc -l < ${SAMPLE_ID}_sorted.txt)
            N=$(( ${N} - 1 ))
            if [ ${N} -ne ${N_RECORDS} ]; then
                echo "Error: file ${SAMPLE_ID}_sorted.txt has ${N} records, but the intersample VCF has ${N_RECORDS} records."
                exit 1
            fi
            head -n 1 ${SAMPLE_ID}_sorted.txt > s.${i}.txt
            tail -n +2 ${SAMPLE_ID}_sorted.txt > b.${i}.txt
            rm -f ${SAMPLE_ID}_sorted.txt
            i=$(( ${i} + 1 ))
        done < ~{samples_file}
        FIELDS_FILES=$( ls s.*.txt | sort --version-sort | tr '\n' ' ' )
        COLUMNS_FILES=$( ls b.*.txt | sort --version-sort | tr '\n' ' ' )
    
        # Pasting the GT files sequentially
        bcftools view --header-only ~{truvari_collapse_intersample_vcf_gz} > tmp.txt
        N_ROWS=$(wc -l < tmp.txt)
        head -n $(( ${N_ROWS} - 1 )) tmp.txt > header.txt
        tail -n 1 tmp.txt | cut -f 1-9 > fields.txt
        paste fields.txt ${FIELDS_FILES} > fields_all.txt
        cat fields_all.txt
        rm -f ${FIELDS_FILES}
        date
        bcftools view --threads ${N_THREADS} --no-header ~{truvari_collapse_intersample_vcf_gz} | awk -v format_field=~{format_field} -F '\t' 'BEGIN { } { printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,format_field); }' > calls.txt
        ls -lh calls.txt
        date
        cat header.txt fields_all.txt > out.vcf
        ${TIME_COMMAND} paste calls.txt ${COLUMNS_FILES} >> out.vcf
        ls -lh out.vcf
        rm -f calls.txt ${COLUMNS_FILES}
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} --compress-level 1 out.vcf
        ${TIME_COMMAND} tabix -f out.vcf.gz
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp out.vcf.'gz*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading output. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>

    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
