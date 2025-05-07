version 1.0


# Pastes all the GT column files coming from the inter-sample-re-genotyped VCF
# into a new inter-sample VCF that holds the updated genotypes.
#
# Remark: all GT column files have the same number of records, and records
# appear in the same order, which is identical to their order in
# `truvari_collapse_intersample_vcf_gz`.
#
workflow Workpackage11 {
    input {
        File truvari_collapse_intersample_vcf_gz
        File truvari_collapse_intersample_tbi
        
        String remote_indir
        String remote_outdir
        File samples_file
        
        Int n_cpu = 4
        Int ram_size_gb = 128
        Int disk_size_gb = 512
    }
    parameter_meta {
        truvari_collapse_intersample_vcf_gz: "The output of `Workpackage8.wdl`. The FORMAT field of this file is assumed to be already the correct one: `GT:FT:SQ:GQ:PS:NE:DP:AD:KS:SCORE:CALIBRATION_SENSITIVITY:SUPP_PBSV:SUPP_SNIFFLES:SUPP_PAV`."
        remote_indir: "Contains GT column files, whose rows are in the same order as `truvari_collapse_intersample_vcf_gz`."
        samples_file: "Order in which to store the samples in the output VCF."
    }
    
    call Workpackage11Impl {
        input:
            truvari_collapse_intersample_vcf_gz = truvari_collapse_intersample_vcf_gz,
            truvari_collapse_intersample_tbi = truvari_collapse_intersample_tbi,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            samples_file = samples_file,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}


# Performance on 10023 samples, 15x, GRCh38, <=0.7:
#
# TOOL                      CPU     RAM     TIME
# paste global              
# bgzip 1
# CopyFormat                
# bgzip 2
#
task Workpackage11Impl {
    input {
        File truvari_collapse_intersample_vcf_gz
        File truvari_collapse_intersample_tbi
        
        String remote_indir
        String remote_outdir
        File samples_file
        
        Int n_cpu
        Int ram_size_gb
        Int disk_size_gb
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
        while read SAMPLE_ID; do
            N=$(wc -l < ${SAMPLE_ID}_sorted.txt)
            N=$(( ${N} - 1 ))
            if [ ${N} -ne ${N_RECORDS} ]; then
                echo "Error: file ${SAMPLE_ID}_sorted.txt has ${N} records, but the intersample VCF has ${N_RECORDS} records."
                exit 1
            fi
        done < ~{samples_file}
        
        # Pasting the GT files sequentially
        bcftools view --header-only ~{truvari_collapse_intersample_vcf_gz} > tmp.txt
        N_ROWS=$(wc -l < tmp.txt)
        head -n $(( ${N_ROWS} - 1 )) tmp.txt > header.txt
        tail -n 1 tmp.txt | cut -f 1-9 > fields.txt
        COLUMNS_FILES=""; FIELDS_FILES=""
        while read SAMPLE_ID; do
            head -n 1 ${SAMPLE_ID}_sorted.txt > s_${i}.txt
            FIELDS_FILES="${FIELDS_FILES} s_${i}.txt"
            tail -n +2 ${SAMPLE_ID}_sorted.txt > b_${i}.txt
            COLUMNS_FILES="${COLUMNS_FILES} b_${i}.txt"
            rm -f ${SAMPLE_ID}_sorted.txt
        done < ~{samples_file}
        paste fields.txt ${FIELDS_FILES} > fields_all.txt
        rm -f ${FIELDS_FILES}
        cat fields_all.txt
        date
        bcftools view --threads ${N_THREADS} --no-header ~{truvari_collapse_intersample_vcf_gz} | cut -f 1-9 > calls.txt
        date
        ${TIME_COMMAND} paste calls.txt ${COLUMNS_FILES} > body.txt
        rm -f ${COLUMNS_FILES}
        head -n 10 body.txt
        cat header.txt fields_all.txt body.txt > merged.vcf
        rm -f header.txt fields_all.txt body.txt
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} merged.vcf
        tabix -f merged.vcf.gz
        
        # Copying from `truvari_collapse_intersample_vcf_gz` the FORMAT fields
        # that are missing in `merged.vcf.gz`.
        ${TIME_COMMAND} java -Xmx${EFFECTIVE_MEM_GB}G -cp ~{docker_dir} CopyFormat ~{truvari_collapse_intersample_vcf_gz} merged.vcf.gz > truvari_collapse_intersample_regenotyped.vcf
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} truvari_collapse_intersample_regenotyped.vcf
        tabix -f truvari_collapse_intersample_regenotyped.vcf.gz
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp truvari_collapse_intersample_regenotyped.vcf.'gz*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading the regenotyped file. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>

    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}