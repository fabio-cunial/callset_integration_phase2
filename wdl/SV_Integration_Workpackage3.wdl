version 1.0


# XGBoost scoring
#
workflow SV_Integration_Workpackage3 {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir
        
        Array[String] annotations = ["KS_1","KS_2","SQ","GQ","DP","AD_NON_ALT","AD_ALL","GT_COUNT","SUPP_PAV","SUPP_SNIFFLES","SUPP_PBSV","SVLEN"]
        
        File training_resource_bed
        File training_python_script
        File scoring_python_script
        File hyperparameters_json
    }
    parameter_meta {
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
        training_resource_bed: "The same BED used in `SV_Integration_Workpackage2`."
        hyperparameters_json: "Parameters for `gatk TrainVariantAnnotationsModel`."
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            training_resource_bed = training_resource_bed,
            annotations = annotations,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json
    }
    
    output {
    }
}


#
task Impl {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir
        
        File training_resource_bed

        Array[String] annotations
        File training_python_script
        File scoring_python_script
        File hyperparameters_json
        
        Int n_cpu = 4
        Int ram_size_gb = 12
        Int disk_size_gb = 20
    }
    parameter_meta {
    }
    
    String docker_dir = "/root"
    
    command <<<
        set -euxo pipefail
        
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
                TEST=$(gcloud storage cp ${REMOTE_DIR}/${SAMPLE_ID}_preprocessed.vcf.'gz*' . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${REMOTE_DIR}/${SAMPLE_ID}_preprocessed.vcf.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gcloud storage cp ${REMOTE_DIR}/${SAMPLE_ID}_training.vcf.'gz*' . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${REMOTE_DIR}/${SAMPLE_ID}_training.vcf.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        }
        
        
        # Deletes all files and directories related to the sample.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -rf ./${SAMPLE_ID}_*
        }
        
        
        function JointVcfFiltering() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local RESOURCE_VCF_GZ=$3
            
            gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ExtractVariantAnnotations -V ${INPUT_VCF_GZ} -O ${SAMPLE_ID}_extract -A ~{sep=" -A " annotations} --resource:resource,training=true,calibration=true ${RESOURCE_VCF_GZ} --maximum-number-of-unlabeled-variants 10000000 --mode INDEL --mnp-type INDEL -L ~{training_resource_bed}
            ls -laht
            # Output:
            # ${SAMPLE_ID}_extract.annot.hdf5
            # ${SAMPLE_ID}_extract.unlabeled.annot.hdf5
            # ${SAMPLE_ID}_extract.vcf.gz
            # ${SAMPLE_ID}_extract.vcf.gz.tbi
            gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" TrainVariantAnnotationsModel --annotations-hdf5 ${SAMPLE_ID}_extract.annot.hdf5 --unlabeled-annotations-hdf5 ${SAMPLE_ID}_extract.unlabeled.annot.hdf5 --model-backend PYTHON_SCRIPT --python-script ~{training_python_script} --hyperparameters-json ~{hyperparameters_json} -O ${SAMPLE_ID}.train --mode INDEL --verbosity DEBUG
            ls -laht
            # Output: 
            # ${SAMPLE_ID}.train.*
            gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ScoreVariantAnnotations -V ${INPUT_VCF_GZ} -O ${SAMPLE_ID}_score -A ~{sep=" -A " annotations} --resource:resource,training=true,calibration=true ${RESOURCE_VCF_GZ} --resource:extracted,extracted=true ${SAMPLE_ID}_extract.vcf.gz --model-prefix ${SAMPLE_ID}.train --model-backend PYTHON_SCRIPT --python-script ~{scoring_python_script} --mode INDEL --mnp-type INDEL --ignore-all-filters --verbosity DEBUG
            ls -laht
            # Output:
            # ${SAMPLE_ID}_score.vcf.gz
            # ${SAMPLE_ID}_score.vcf.gz.tbi
            # ${SAMPLE_ID}_score.annot.hdf5
            # ${SAMPLE_ID}_score.scores.hdf5
            
            # Removing temporary files
            rm -f ${SAMPLE_ID}_extract.annot.hdf5 ${SAMPLE_ID}_extract.unlabeled.annot.hdf5 ${SAMPLE_ID}_extract.vcf.gz* ${SAMPLE_ID}.train.* ${SAMPLE_ID}_score.annot.hdf5 ${SAMPLE_ID}_score.scores.hdf5
        }
        
        
        # Copies SCORE from INFO to FORMAT
        #
        function CopyInfoToFormat() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            

            echo '##FORMAT=<ID=SCORE,Number=1,Type=Float,Description="Score according to the XGBoost model">' > ${SAMPLE_ID}_header.txt
            bcftools query --format '%CHROM\t%POS\t%ID\t%SCORE\n' ${INPUT_VCF_GZ} | bgzip -c > ${SAMPLE_ID}_format.tsv.gz
            tabix -s1 -b2 -e2 ${SAMPLE_ID}_format.tsv.gz
            bcftools annotate --threads ${N_THREADS} --header-lines ${SAMPLE_ID}_header.txt --annotations ${SAMPLE_ID}_format.tsv.gz --columns CHROM,POS,~ID,FORMAT/SCORE --output-type z ${INPUT_VCF_GZ} > ${SAMPLE_ID}_scored.vcf.gz
            bcftools index --threads ${N_THREADS} --tbi ${SAMPLE_ID}_scored.vcf.gz
            bcftools view --no-header ${SAMPLE_ID}_scored.vcf.gz | head -n 5 || echo "0"
            
            # Removing temporary files
            rm -f ${SAMPLE_ID}_header.txt ${SAMPLE_ID}_format.tsv.gz*
        }

        

        
        # ---------------------------- Main program ----------------------------
        
        export GATK_LOCAL_JAR="/root/gatk.jar"
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            LocalizeSample ${SAMPLE_ID} ~{remote_indir}
            
            JointVcfFiltering ${SAMPLE_ID} ${SAMPLE_ID}_preprocessed.vcf.gz ${SAMPLE_ID}_training.vcf.gz
            CopyInfoToFormat ${SAMPLE_ID} ${SAMPLE_ID}_score.vcf.gz
            
            gcloud storage mv ${SAMPLE_ID}_scored.vcf.'gz*' ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done < chunk.csv
    >>>
    
    output {
    }
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
