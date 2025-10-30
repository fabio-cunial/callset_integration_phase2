version 1.0


# Splits every XGBoost-scored single-sample VCF into ~100 pieces in order to run
# bcftools merge in parallel.
#
workflow SV_Integration_Workpackage4 {
    input {
        File sv_integration_chunk_tsv
        File split_for_bcftools_merge_csv
        
        String remote_indir
        String remote_outdir
    }
    parameter_meta {
        split_for_bcftools_merge_csv: "A partition that covers all chromosomes. Every line is a 0-based, half-open, consecutive chunk of a chromosome. Lines are assumed to be sorted."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


#
task Impl {
    input {
        File sv_integration_chunk_tsv
        File split_for_bcftools_merge_csv
        
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 2
        Int ram_size_gb = 4
        Int disk_size_gb = 20
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local REMOTE_DIR=$2
            
            while : ; do
                TEST=$(gsutil -m cp ${REMOTE_DIR}/${SAMPLE_ID}_scored.vcf.'gz*' . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${REMOTE_DIR}/${SAMPLE_ID}_scored.vcf.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        }
        
        
        # Deletes all files and directories related to the sample
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -rf ./${SAMPLE_ID}_*
        }
        
        
        # Transfers INFO annotations to FORMAT, so that they are preserved by
        # the inter-sample merge later.
        #
        # Remark: SCORE has already been transfered to FORMAT by
        # `Workpackage3.wdl`
        #
        function CopyInfoToFormat() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            echo '##FORMAT=<ID=CALIBRATION_SENSITIVITY,Number=1,Type=Float,Description="Calibration sensitivity according to the model applied by ScoreVariantAnnotations">' > ${SAMPLE_ID}_header.txt
            echo '##FORMAT=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> ${SAMPLE_ID}_header.txt
            echo '##FORMAT=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by sniffles">' >> ${SAMPLE_ID}_header.txt
            echo '##FORMAT=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by pav">' >> ${SAMPLE_ID}_header.txt
            bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CALIBRATION_SENSITIVITY\t%SUPP_PBSV\t%SUPP_SNIFFLES\t%SUPP_PAV\n' ${INPUT_VCF_GZ} | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
            tabix -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --header-lines ${SAMPLE_ID}_header.txt --annotations ${SAMPLE_ID}_annotations.tsv.gz --columns CHROM,POS,~ID,REF,ALT,FORMAT/CALIBRATION_SENSITIVITY,FORMAT/SUPP_PBSV,FORMAT/SUPP_SNIFFLES,FORMAT/SUPP_PAV ${INPUT_VCF_GZ} --output-type z > ${SAMPLE_ID}_annotated.vcf.gz
            tabix -f ${SAMPLE_ID}_annotated.vcf.gz
        }
        
        
        function ChunkVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            i="0"
            while read INTERVAL; do
                echo ${INTERVAL} | tr ',' '\t' > ${SAMPLE_ID}.bed
                ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --regions-file ${SAMPLE_ID}.bed --regions-overlap pos --output-type z ${INPUT_VCF_GZ} > ${SAMPLE_ID}_chunk_${i}.vcf.gz
                tabix -f ${SAMPLE_ID}_chunk_${i}.vcf.gz
                i=$(( ${i} + 1 ))
            done < ~{split_for_bcftools_merge_csv}
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${SAMPLE_ID}_chunk_'*'.vcf.'gz*' ~{remote_outdir}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading chunks for sample ${SAMPLE_ID}. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------

        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            LocalizeSample ${SAMPLE_ID} ~{remote_indir}
            
            CopyInfoToFormat ${SAMPLE_ID} ${SAMPLE_ID}_scored.vcf.gz
            ChunkVcf ${SAMPLE_ID} ${SAMPLE_ID}_annotated.vcf.gz
            
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done < chunk.csv
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
