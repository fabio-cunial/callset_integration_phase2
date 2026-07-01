version 1.0


# ---------->
#
workflow SV_Integration_BndScore {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File resource_vcf_gz
        File resource_vcf_gz_tbi
        String remote_outdir
                
        String exclude_chromosomes_string = " "
        File reference_fa
        File reference_fai

        File training_python_script
        File scoring_python_script
        File hyperparameters_json
        
        String docker_image = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb"
    }
    parameter_meta {
        input_vcf_gz: "Assumed to contain all the annotations used in this workflow."
        remote_outdir: "Without final slash"
        exclude_chromosomes_string: "Example: -XL chr1 -XL chr2 -XL chr3 -XL chr4 -XL chr5"
        hyperparameters_json: "Parameters for `gatk TrainVariantAnnotationsModel`."
    }
    
    
    call Score as sniffles {
        input:
            caller_id = "sniffles",
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            resource_vcf_gz = resource_vcf_gz,
            resource_vcf_gz_tbi = resource_vcf_gz_tbi,
            remote_outdir = remote_outdir,
            exclude_chromosomes_string = exclude_chromosomes_string,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
    }
    call Score as pbsv {
        input:
            caller_id = "pbsv",
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            resource_vcf_gz = resource_vcf_gz,
            resource_vcf_gz_tbi = resource_vcf_gz_tbi,
            remote_outdir = remote_outdir,
            exclude_chromosomes_string = exclude_chromosomes_string,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
    }
    
    output {
    }
}


#
task Score {
    input {
        String caller_id
        
        File input_vcf_gz
        File input_vcf_gz_tbi
        File resource_vcf_gz
        File resource_vcf_gz_tbi
        String remote_outdir
        
        String exclude_chromosomes_string
        File reference_fa
        File reference_fai

        File training_python_script
        File scoring_python_script
        File hyperparameters_json
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 10
        Int disk_size_gb = 10
    }
    parameter_meta {
    }
    
    String docker_dir = "/root"
    
    command <<<
        set -euxo pipefail
        
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        export GATK_LOCAL_JAR="/root/gatk.jar"
        



        # ----------------------- Steps of the pipeline ------------------------

        # Extracts the BNDs created by a given caller and annotates them with
        # corresponding fields for training.
        #
        function ExtractBndsOfCaller() {
            local INPUT_VCF_GZ=$1
            local CALLER_ID=$2
            local OUTPUT_VCF_GZ=$3

            bcftools filter --include 'ID~"'${CALLER_ID}'"' --output-type v ${INPUT_VCF_GZ} --output ${CALLER_ID}.vcf
            if [ ${CALLER_ID} = "sniffles" ]; then
                bcftools query --format '%CHROM\t%POS\t%ID\t%INFO/SUPPORT\t%INFO/COVERAGE\t%INFO/STRAND\t%INFO/STDEV_POS\t[%GT]\t[%GQ]\t[%DR]\t[%DV]\n' ${CALLER_ID}.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    STRAND=-1; \
                    if ($6=="--") STRAND=0; \
                    else if ($6=="-+") STRAND=1; \
                    else if ($6=="+-") STRAND=2; \
                    else if ($6=="++") STRAND=3; \
                    GT_COUNT=-1; \
                    if ($8=="0/0" || $8=="0|0" || $8=="./."  || $8==".|." || $8=="./0" || $8==".|0" || $8=="0/." || $8=="0|." || $8=="0" || $8==".") GT_COUNT=0; \
                    else if ($8=="0/1" || $8=="0|1" || $8=="1/0" || $8=="1|0" || $8=="./1" || $8==".|1" || $8=="1/." || $8=="1|." || $8=="1") GT_COUNT=1; \
                    else if ($8=="1/1" || $8=="1|1") GT_COUNT=2; \
                    printf("%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,STRAND,$7,GT_COUNT,$9,$10,$11); \
                }' | tr ',' '\t' | bgzip -c > ${CALLER_ID}_annotations.tsv.gz
                echo '##INFO=<ID=SNIFFLES_SUPPORT,Number=1,Type=Integer,Description="...">' > ${CALLER_ID}_header.txt
                echo '##INFO=<ID=SNIFFLES_COVERAGE_1,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=SNIFFLES_COVERAGE_2,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=SNIFFLES_COVERAGE_3,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=SNIFFLES_COVERAGE_4,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=SNIFFLES_COVERAGE_5,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=SNIFFLES_STRAND,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=SNIFFLES_STDEV_POS,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=SNIFFLES_GT_COUNT,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=SNIFFLES_GQ,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=SNIFFLES_DR,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=SNIFFLES_DV,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                local COLUMNS='CHROM,POS,~ID,INFO/SNIFFLES_SUPPORT,INFO/SNIFFLES_COVERAGE_1,INFO/SNIFFLES_COVERAGE_2,INFO/SNIFFLES_COVERAGE_3,INFO/SNIFFLES_COVERAGE_4,INFO/SNIFFLES_COVERAGE_5,INFO/SNIFFLES_STRAND,INFO/SNIFFLES_STDEV_POS,INFO/SNIFFLES_GT_COUNT,INFO/SNIFFLES_GQ,INFO/SNIFFLES_DR,INFO/SNIFFLES_DV'
            else
                bcftools query --format '%CHROM\t%POS\t%ID\t%INFO/CIPOS\t[%GT]\t[%AD]\t[%DP]\n' ${CALLER_ID}.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    GT_COUNT=-1; \
                    if ($5=="0/0" || $5=="0|0" || $5=="./."  || $5==".|." || $5=="./0" || $5==".|0" || $5=="0/." || $5=="0|." || $5=="0" || $5==".") GT_COUNT=0; \
                    else if ($5=="0/1" || $5=="0|1" || $5=="1/0" || $5=="1|0" || $5=="./1" || $5==".|1" || $5=="1/." || $5=="1|." || $5=="1") GT_COUNT=1; \
                    else if ($5=="1/1" || $5=="1|1") GT_COUNT=2; \
                    printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,GT_COUNT,$6,$7); \
                }' | tr ',' '\t' | bgzip -c > ${CALLER_ID}_annotations.tsv.gz
                echo '##INFO=<ID=PBSV_CIPOS_1,Number=1,Type=Integer,Description="...">' > ${CALLER_ID}_header.txt
                echo '##INFO=<ID=PBSV_CIPOS_2,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=PBSV_GT_COUNT,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=PBSV_AD_1,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=PBSV_AD_2,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                echo '##INFO=<ID=PBSV_DP,Number=1,Type=Integer,Description="...">' >> ${CALLER_ID}_header.txt
                local COLUMNS='CHROM,POS,~ID,INFO/PBSV_CIPOS_1,INFO/PBSV_CIPOS_2,INFO/PBSV_GT_COUNT,INFO/PBSV_AD_1,INFO/PBSV_AD_2,INFO/PBSV_DP'
            fi
            tabix -f -s1 -b2 -e2 ${CALLER_ID}_annotations.tsv.gz
            bcftools annotate --annotations ${CALLER_ID}_annotations.tsv.gz --header-lines ${CALLER_ID}_header.txt --columns ${COLUMNS} --output-type z ${CALLER_ID}.vcf --output ${OUTPUT_VCF_GZ}
            rm -f ${CALLER_ID}.vcf ${CALLER_ID}_annotations.tsv.gz ${CALLER_ID}_header.txt
            bcftools index -f -t ${OUTPUT_VCF_GZ}
        }




        # ---------------------------- Main program ----------------------------

        ExtractBndsOfCaller ~{input_vcf_gz} ~{caller_id} input_cleaned.vcf.gz
        ExtractBndsOfCaller ~{resource_vcf_gz} ~{caller_id} resource_cleaned.vcf.gz
        if [ ~{caller_id} = "sniffles" ]; then
            ANNOTATIONS_STRING="-A SNIFFLES_SUPPORT -A SNIFFLES_COVERAGE_1 -A SNIFFLES_COVERAGE_2 -A SNIFFLES_COVERAGE_3 -A SNIFFLES_COVERAGE_4 -A SNIFFLES_COVERAGE_5 -A SNIFFLES_STRAND -A SNIFFLES_STDEV_POS -A SNIFFLES_GT_COUNT -A SNIFFLES_GQ -A SNIFFLES_DR -A SNIFFLES_DV"
        else
            ANNOTATIONS_STRING="-A PBSV_CIPOS_1 -A PBSV_CIPOS_2 -A PBSV_GT_COUNT -A PBSV_AD_1 -A PBSV_AD_2 -A PBSV_DP"
        fi

        gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ExtractVariantAnnotations -V input_cleaned.vcf.gz ~{exclude_chromosomes_string} -O extract ${ANNOTATIONS_STRING} --resource:resource,training=true,calibration=true resource_cleaned.vcf.gz --maximum-number-of-unlabeled-variants 1000000000 --mode INDEL --mnp-type INDEL
        ls -laht 1>&2
        # Output:
        # extract.annot.hdf5
        # extract.unlabeled.annot.hdf5
        # extract.vcf.gz
        # extract.vcf.gz.tbi
        gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" TrainVariantAnnotationsModel --annotations-hdf5 extract.annot.hdf5 --unlabeled-annotations-hdf5 extract.unlabeled.annot.hdf5 --model-backend PYTHON_SCRIPT --python-script ~{training_python_script} --hyperparameters-json ~{hyperparameters_json} -O train.train --mode INDEL --verbosity DEBUG
        ls -laht 1>&2
        # Output: 
        # train.train.indel.unlabeledScores.hdf5
        # train.train.indel.calibrationScores.hdf5
        # train.train.indel.trainingScores.hdf5
        # train.train.indel.scorer.pkl
        gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ScoreVariantAnnotations -V input_cleaned.vcf.gz -O score ${ANNOTATIONS_STRING} --resource:resource,training=true,calibration=true resource_cleaned.vcf.gz --resource:extracted,extracted=true extract.vcf.gz --model-prefix train.train --model-backend PYTHON_SCRIPT --python-script ~{scoring_python_script} --mode INDEL --mnp-type INDEL --ignore-all-filters --verbosity DEBUG
        ls -laht 1>&2
        # Output:
        # score.vcf.gz
        # score.vcf.gz.tbi
        # score.annot.hdf5
        # score.scores.hdf5
        gsutil -m mv score.vcf.gz ~{remote_outdir}/~{caller_id}_score.vcf.gz
        gsutil -m mv score.vcf.gz.tbi ~{remote_outdir}/~{caller_id}_score.vcf.gz.tbi
    >>>
    
    output {
        File extract_annot_hdf5 = "extract.annot.hdf5"
        File extract_unlabeled_annot_hdf5 = "extract.unlabeled.annot.hdf5"

        File train_indel_unlabeled_scores_hdf5 = "train.train.indel.unlabeledScores.hdf5"
        File train_indel_calibration_scores_hdf5 = "train.train.indel.calibrationScores.hdf5"
        File train_indel_training_scores_hdf5 = "train.train.indel.trainingScores.hdf5"
        File train_indel_scorer_pkl = "train.train.indel.scorer.pkl"

        File score_annot_hdf5 = "score.annot.hdf5"
        File score_scores_hdf5 = "score.scores.hdf5"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
