version 1.0


# A version of `SV_Integration_Workpackage3.wdl` that takes in input the 
# ultralong calls annotated by `SV_Integration_UltralongAnnotate.wdl`, filters
# them using global models, and produces output compatible with 
# `SV_Integration_Workpackage12.wdl`.
#
# Remark: three VCFs are produced for each sample: a lenient one, a stringent 
# one, and an unfiltered one that just contains model scores.
#
# Remark: INSDUPs are converted back to INS, to preserve as much information as 
# possible downstream. The breakpoints of the original INSDUPs are saved in
# INFO.
#
workflow SV_Integration_Workpackage3_Ultralong {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir_lenient
        String remote_outdir_stringent
        String remote_outdir_all

        Int remove_symbolic_ins = 1

        File del_indel_scorer_15x_pkl
        File ins_indel_scorer_15x_pkl
        File dup_indel_scorer_15x_pkl
        File insdup_indel_scorer_15x_pkl
        File inv_indel_scorer_15x_pkl

        File del_indel_calibrationScores_15x_hdf5
        File ins_indel_calibrationScores_15x_hdf5
        File dup_indel_calibrationScores_15x_hdf5
        File insdup_indel_calibrationScores_15x_hdf5
        File inv_indel_calibrationScores_15x_hdf5

        File del_indel_scorer_30x_pkl
        File ins_indel_scorer_30x_pkl
        File dup_indel_scorer_30x_pkl
        File insdup_indel_scorer_30x_pkl
        File inv_indel_scorer_30x_pkl

        File del_indel_calibrationScores_30x_hdf5
        File ins_indel_calibrationScores_30x_hdf5
        File dup_indel_calibrationScores_30x_hdf5
        File insdup_indel_calibrationScores_30x_hdf5
        File inv_indel_calibrationScores_30x_hdf5

        File sample_coverages_csv

        Array[String] annotations_interval = [ "GT_COUNT","SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
                                               "BIN_BEFORE_COVERAGE","BIN_LEFT_COVERAGE","BIN_1_COVERAGE","BIN_2_COVERAGE","BIN_3_COVERAGE","BIN_4_COVERAGE","BIN_5_COVERAGE","BIN_6_COVERAGE","BIN_7_COVERAGE","BIN_8_COVERAGE","BIN_9_COVERAGE","BIN_10_COVERAGE","BIN_RIGHT_COVERAGE","BIN_AFTER_COVERAGE","BIN_LEFT_MAPQ","BIN_RIGHT_MAPQ","BIN_LEFT_SECONDARY","BIN_RIGHT_SECONDARY","LL","LR","RL","RR","LL_RL_1","LL_RL_2","LL_RL_3","LL_RL_4","LL_RR_1","LL_RR_2","LL_RR_3","LL_RR_4","LR_RL_1","LR_RL_2","LR_RL_3","LR_RL_4","LR_RR_1","LR_RR_2","LR_RR_3","LR_RR_4",
                                               "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK" 
                                             ]
        Array[String] annotations_point = [ "GT_COUNT","SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV",
                                            "BIN_POS","BIN_POINT_MAPQ","BIN_POINT_SECONDARY","PL","PR","PL_PL_1","PL_PL_2","PL_PL_3","PL_PL_4","PL_PR_1","PL_PR_2","PL_PR_3","PL_PR_4","PR_PR_1","PR_PR_2","PR_PR_3","PR_PR_4",
                                            "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK"
                                        ]
        Int annotations_have_gt_count = 1
        File scoring_python_script

        String filter_string_lenient
        String filter_string_stringent
        
        String docker_image = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb"
        File UltralongInsdups2Ins_java
        File AddSvlenToSymbolicAlt_java
    }
    parameter_meta {
        remote_indir: "Without final slash"
        remote_outdir_lenient: "Without final slash"
        filter_string_lenient: "Example: FORMAT/CALIBRATION_SENSITIVITY<=0.9"
        sample_coverages_csv: "One line per sample, with columns: `SAMPLE_ID,COVERAGE`. Used to select the appropriate model for each sample."
        remove_symbolic_ins: "It might happen that some records with ALT=<INS> were still present in the VCFs in input to `SV_Integration_UltralongAnnotate.wdl`. The current version of the pipeline supports them but an older version didn't: for the latter, we can discard them here before further processing."
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            remote_indir = remote_indir,
            remote_outdir_lenient = remote_outdir_lenient,
            remote_outdir_stringent = remote_outdir_stringent,
            remote_outdir_all = remote_outdir_all,

            remove_symbolic_ins = remove_symbolic_ins,

            del_indel_scorer_15x_pkl = del_indel_scorer_15x_pkl,
            ins_indel_scorer_15x_pkl = ins_indel_scorer_15x_pkl,
            dup_indel_scorer_15x_pkl = dup_indel_scorer_15x_pkl,
            insdup_indel_scorer_15x_pkl = insdup_indel_scorer_15x_pkl,
            inv_indel_scorer_15x_pkl = inv_indel_scorer_15x_pkl,

            del_indel_scorer_30x_pkl = del_indel_scorer_30x_pkl,
            ins_indel_scorer_30x_pkl = ins_indel_scorer_30x_pkl,
            dup_indel_scorer_30x_pkl = dup_indel_scorer_30x_pkl,
            insdup_indel_scorer_30x_pkl = insdup_indel_scorer_30x_pkl,
            inv_indel_scorer_30x_pkl = inv_indel_scorer_30x_pkl,

            del_indel_calibrationScores_15x_hdf5 = del_indel_calibrationScores_15x_hdf5,
            ins_indel_calibrationScores_15x_hdf5 = ins_indel_calibrationScores_15x_hdf5,
            dup_indel_calibrationScores_15x_hdf5 = dup_indel_calibrationScores_15x_hdf5,
            insdup_indel_calibrationScores_15x_hdf5 = insdup_indel_calibrationScores_15x_hdf5,
            inv_indel_calibrationScores_15x_hdf5 = inv_indel_calibrationScores_15x_hdf5,

            del_indel_calibrationScores_30x_hdf5 = del_indel_calibrationScores_30x_hdf5,
            ins_indel_calibrationScores_30x_hdf5 = ins_indel_calibrationScores_30x_hdf5,
            dup_indel_calibrationScores_30x_hdf5 = dup_indel_calibrationScores_30x_hdf5,
            insdup_indel_calibrationScores_30x_hdf5 = insdup_indel_calibrationScores_30x_hdf5,
            inv_indel_calibrationScores_30x_hdf5 = inv_indel_calibrationScores_30x_hdf5,

            sample_coverages_csv = sample_coverages_csv,

            annotations_interval = annotations_interval,
            annotations_point = annotations_point,
            annotations_have_gt_count = annotations_have_gt_count,
            scoring_python_script = scoring_python_script,

            filter_string_lenient = filter_string_lenient,
            filter_string_stringent = filter_string_stringent,
            
            docker_image = docker_image,
            UltralongInsdups2Ins_java = UltralongInsdups2Ins_java,
            AddSvlenToSymbolicAlt_java = AddSvlenToSymbolicAlt_java
    }
    
    output {
    }
}


# Remark: we use gsutil instead of gcloud since we found the latter to have
# issues in practice (maybe the gcloud version in the docker is not up to
# date?). 
#
# Memory bottlenecks (measured on a 4GB VM):
#
# ScoreVariantAnnotations               200 MB
#
task Impl {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir_lenient
        String remote_outdir_stringent
        String remote_outdir_all

        Int remove_symbolic_ins

        File del_indel_scorer_15x_pkl
        File ins_indel_scorer_15x_pkl
        File dup_indel_scorer_15x_pkl
        File insdup_indel_scorer_15x_pkl
        File inv_indel_scorer_15x_pkl

        File del_indel_calibrationScores_15x_hdf5
        File ins_indel_calibrationScores_15x_hdf5
        File dup_indel_calibrationScores_15x_hdf5
        File insdup_indel_calibrationScores_15x_hdf5
        File inv_indel_calibrationScores_15x_hdf5

        File del_indel_scorer_30x_pkl
        File ins_indel_scorer_30x_pkl
        File dup_indel_scorer_30x_pkl
        File insdup_indel_scorer_30x_pkl
        File inv_indel_scorer_30x_pkl

        File del_indel_calibrationScores_30x_hdf5
        File ins_indel_calibrationScores_30x_hdf5
        File dup_indel_calibrationScores_30x_hdf5
        File insdup_indel_calibrationScores_30x_hdf5
        File inv_indel_calibrationScores_30x_hdf5

        File sample_coverages_csv

        Array[String] annotations_interval
        Array[String] annotations_point
        Int annotations_have_gt_count
        File scoring_python_script

        String filter_string_lenient
        String filter_string_stringent
        
        String docker_image
        File UltralongInsdups2Ins_java
        File AddSvlenToSymbolicAlt_java
        Int n_cpu = 2
        Int ram_size_gb = 4
        Int disk_size_gb = 20
        Int preemptible_number = 4
    }
    parameter_meta {
    }
    
    String docker_dir = "/root"
    
    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export GATK_LOCAL_JAR="/root/gatk.jar"
        RAM_PER_THREAD_MB=$(( ~{ram_size_gb} * 1024 - 500 ))
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local REMOTE_DIR=$2
            
            gsutil cp ${REMOTE_DIR}/${SAMPLE_ID}_del.vcf.'gz*' ${REMOTE_DIR}/${SAMPLE_ID}_ins.vcf.'gz*' ${REMOTE_DIR}/${SAMPLE_ID}_dup.vcf.'gz*' ${REMOTE_DIR}/${SAMPLE_ID}_insdup.vcf.'gz*' ${REMOTE_DIR}/${SAMPLE_ID}_inv.vcf.'gz*' .
        }
        

        # Removes records that had ALT=<INS> in the VCFs in input to 
        # `SV_Integration_UltralongAnnotate.wdl`.
        #
        function RemoveSymbolicIns() {
            local SAMPLE_ID=$1

            N_RECORDS_BEFORE=$(bcftools index --nrecords ${SAMPLE_ID}_ins.vcf.gz.tbi)
            bcftools filter --exclude 'ALT="<INS>"' --output-type z ${SAMPLE_ID}_ins.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_ins.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_ins.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_ins.vcf.gz
            N_RECORDS_AFTER=$(bcftools index --nrecords ${SAMPLE_ID}_ins.vcf.gz.tbi)

            N_RECORDS_BEFORE=$(bcftools index --nrecords ${SAMPLE_ID}_insdup.vcf.gz.tbi)
            bcftools filter --exclude 'INS_ALT="<INS>"' --output-type z ${SAMPLE_ID}_insdup.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_insdup.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_insdup.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_insdup.vcf.gz
            N_RECORDS_AFTER=$(bcftools index --nrecords ${SAMPLE_ID}_insdup.vcf.gz.tbi)
        }

        
        # Deletes all files and directories related to the sample.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -rf ./${SAMPLE_ID}_*
        }


        # Adds field `INFO/GT_COUNT` to the input VCF, which is overwritten.
        #
        function AddGtCount() {
            local SAMPLE_ID=$1
            local SVTYPE=$2

            bcftools query --format '%CHROM\t%POS\t%ID\t[%GT]\n' ${SAMPLE_ID}_${SVTYPE}.vcf.gz | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                GT_COUNT=-1; \
                if ($4=="0/0" || $4=="0|0" || $4=="./."  || $4==".|." || $4=="./0" || $4==".|0" || $4=="0/." || $4=="0|." || $4=="0" || $4==".") GT_COUNT=0; \
                else if ($4=="0/1" || $4=="0|1" || $4=="1/0" || $4=="1|0" || $4=="./1" || $4==".|1" || $4=="1/." || $4=="1|." || $4=="1") GT_COUNT=1; \
                else if ($4=="1/1" || $4=="1|1") GT_COUNT=2; \
                printf("%s\t%d\t%s\t%d\n",$1,$2,$3,GT_COUNT); \
            }' | bgzip -c > ${SAMPLE_ID}_${SVTYPE}_annotations.tsv.gz
            tabix -f -s1 -b2 -e2 ${SAMPLE_ID}_${SVTYPE}_annotations.tsv.gz
            echo '##INFO=<ID=GT_COUNT,Number=1,Type=Integer,Description="Original GT converted to an integer in {0,1,2}.">' > ${SAMPLE_ID}_${SVTYPE}_header.txt
            local COLUMNS='CHROM,POS,~ID,INFO/GT_COUNT'
            bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_${SVTYPE}_annotations.tsv.gz --header-lines ${SAMPLE_ID}_${SVTYPE}_header.txt --columns ${COLUMNS} --output-type z ${SAMPLE_ID}_${SVTYPE}.vcf.gz --output ${SAMPLE_ID}_${SVTYPE}_annotated.vcf.gz
            rm -f ${SAMPLE_ID}_${SVTYPE}_annotations.tsv.gz ${SAMPLE_ID}_${SVTYPE}_header.txt ${SAMPLE_ID}_${SVTYPE}.vcf.gz
            mv ${SAMPLE_ID}_${SVTYPE}_annotated.vcf.gz ${SAMPLE_ID}_${SVTYPE}.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${SVTYPE}.vcf.gz
        }


        # Copies the following fields from INFO to FORMAT, so that they are
        # preserved by the inter-sample merge downstream:
        #
        # SUPP_*, SCORE, CALIBRATION_SENSITIVITY
        #
        # @param 2 A VCF where all IDs are distinct. This is guaranteed by the
        # workpackages upstream.
        #
        function CopyInfoToFormat() {
            local SAMPLE_ID=$1
            local SVTYPE=$2
            
            echo '##FORMAT=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> ${SAMPLE_ID}_${SVTYPE}_header.txt
            echo '##FORMAT=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by sniffles">' >> ${SAMPLE_ID}_${SVTYPE}_header.txt
            echo '##FORMAT=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by pav">' >> ${SAMPLE_ID}_${SVTYPE}_header.txt
            echo '##FORMAT=<ID=SCORE,Number=1,Type=Float,Description="Score according to the XGBoost model">' >> ${SAMPLE_ID}_${SVTYPE}_header.txt
            echo '##FORMAT=<ID=CALIBRATION_SENSITIVITY,Number=1,Type=Float,Description="Calibration sensitivity according to the model applied by ScoreVariantAnnotations">' >> ${SAMPLE_ID}_${SVTYPE}_header.txt
            bcftools query --format '%CHROM\t%POS\t%ID\t%SUPP_PBSV\t%SUPP_SNIFFLES\t%SUPP_PAV\t%SCORE\t%CALIBRATION_SENSITIVITY\n' ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz | bgzip -c > ${SAMPLE_ID}_${SVTYPE}_format.tsv.gz
            tabix -f -s1 -b2 -e2 ${SAMPLE_ID}_${SVTYPE}_format.tsv.gz
            bcftools annotate --threads ${N_THREADS} --header-lines ${SAMPLE_ID}_${SVTYPE}_header.txt --annotations ${SAMPLE_ID}_${SVTYPE}_format.tsv.gz --columns CHROM,POS,~ID,FORMAT/SUPP_PBSV,FORMAT/SUPP_SNIFFLES,FORMAT/SUPP_PAV,FORMAT/SCORE,FORMAT/CALIBRATION_SENSITIVITY --output-type z ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz --output ${SAMPLE_ID}_${SVTYPE}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz* ; mv ${SAMPLE_ID}_${SVTYPE}_out.vcf.gz ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz
            (bcftools view --no-header ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz | head -n 1 || echo "0") 1>&2
            
            # Removing temporary files
            rm -f ${SAMPLE_ID}_${SVTYPE}_header.txt ${SAMPLE_ID}_${SVTYPE}_format.tsv.gz*
        }


        # Assumes that `CopyInfoToFormat()` has already been executed.
        #
        function PrintDebugInformation() {
            local SAMPLE_ID=$1
            local SVTYPE=$2
            
            rm -f ${SAMPLE_ID}_${SVTYPE}_xgboost.csv
            local N_RECORDS_BEFORE_FILTERING=$(bcftools index --nrecords ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz)
            for THRESHOLD in 0.7 0.8 0.9 0.95 ; do
                local N_RECORDS_AFTER_FILTERING=$( bcftools query --format '%ID\n' --include "FORMAT/CALIBRATION_SENSITIVITY<=${THRESHOLD}" ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz | wc -l )
                if [ ${N_RECORDS_BEFORE_FILTERING} -eq 0 ]; then
                    local PERCENT=0
                else
                    local PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER_FILTERING} / ${N_RECORDS_BEFORE_FILTERING}" | bc )
                fi
                echo "${N_RECORDS_AFTER_FILTERING},${N_RECORDS_BEFORE_FILTERING},${PERCENT},Number of records with CALIBRATION_SENSITIVITY<=${THRESHOLD}" >> ${SAMPLE_ID}_${SVTYPE}_xgboost.csv
            done
        }


        function ScoreAndFilter() {
            local SAMPLE_ID=$1
            local COVERAGE=$2
            local SVTYPE=$3
            local RAM_PER_THREAD_MB=$4

            # Scoring
            if [ ${SVTYPE} = "ins" ]; then
                local ANNOTATIONS="~{sep=" -A " annotations_point}"
            else
                local ANNOTATIONS="~{sep=" -A " annotations_interval}"
            fi
            gatk --java-options "-Xmx${RAM_PER_THREAD_MB}m" ScoreVariantAnnotations -V ${SAMPLE_ID}_${SVTYPE}.vcf.gz -O ${SAMPLE_ID}_${SVTYPE}_score -A ${ANNOTATIONS} --model-prefix ${COVERAGE}.${SVTYPE} --model-backend PYTHON_SCRIPT --python-script ~{scoring_python_script} --mode INDEL --mnp-type INDEL --ignore-all-filters --verbosity DEBUG
            CopyInfoToFormat ${SAMPLE_ID} ${SVTYPE}
            PrintDebugInformation ${SAMPLE_ID} ${SVTYPE}

            # Converting INSDUP to INS
            if [ ${SVTYPE} = "insdup" ]; then
                java UltralongInsdups2Ins ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz | bcftools sort --max-mem ${RAM_PER_THREAD_MB}M --output-type z --output ${SAMPLE_ID}_${SVTYPE}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz*
                mv ${SAMPLE_ID}_${SVTYPE}_out.vcf.gz ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz
            fi

            # Filtering
            bcftools view --threads ${N_THREADS} --include "~{filter_string_lenient}"   --output-type b ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz --output ${SAMPLE_ID}_${SVTYPE}_lenient.bcf
            bcftools view --threads ${N_THREADS} --include "~{filter_string_stringent}" --output-type b ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz --output ${SAMPLE_ID}_${SVTYPE}_stringent.bcf
            bcftools view --threads ${N_THREADS}                                        --output-type b ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz --output ${SAMPLE_ID}_${SVTYPE}_all.bcf
            bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_${SVTYPE}_lenient.bcf
            bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_${SVTYPE}_stringent.bcf
            bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_${SVTYPE}_all.bcf
            rm -f ${SAMPLE_ID}_${SVTYPE}_score.vcf.gz*

            # Adding debug information
            local N_RECORDS_BEFORE_FILTERING=$(bcftools index --nrecords ${SAMPLE_ID}_${SVTYPE}_all.bcf)
            local N_RECORDS_AFTER_FILTERING=$(bcftools index --nrecords ${SAMPLE_ID}_${SVTYPE}_lenient.bcf)
            if [ ${N_RECORDS_BEFORE_FILTERING} -eq 0 ]; then
                local PERCENT=0
            else
                local PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER_FILTERING} / ${N_RECORDS_BEFORE_FILTERING}" | bc )
            fi
            echo "${N_RECORDS_AFTER_FILTERING},${N_RECORDS_BEFORE_FILTERING},${PERCENT},Number of records in the lenient VCF" >> ${SAMPLE_ID}_${SVTYPE}_xgboost.csv
            local N_RECORDS_AFTER_FILTERING=$(bcftools index --nrecords ${SAMPLE_ID}_${SVTYPE}_stringent.bcf)
            if [ ${N_RECORDS_BEFORE_FILTERING} -eq 0 ]; then
                local PERCENT=0
            else
                local PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER_FILTERING} / ${N_RECORDS_BEFORE_FILTERING}" | bc )
            fi
            echo "${N_RECORDS_AFTER_FILTERING},${N_RECORDS_BEFORE_FILTERING},${PERCENT},Number of records in the stringent VCF" >> ${SAMPLE_ID}_${SVTYPE}_xgboost.csv
            cat ${SAMPLE_ID}_${SVTYPE}_xgboost.csv 1>&2
        }




        # ---------------------------- Main program ----------------------------

        # Compiling input scripts
        mv ~{UltralongInsdups2Ins_java} UltralongInsdups2Ins.java
        mv ~{AddSvlenToSymbolicAlt_java} AddSvlenToSymbolicAlt.java
        javac *.java

        # Enforcing a consistent naming scheme on all model files
        mv ~{del_indel_scorer_15x_pkl} 15x.del.indel.scorer.pkl
        mv ~{ins_indel_scorer_15x_pkl} 15x.ins.indel.scorer.pkl
        mv ~{dup_indel_scorer_15x_pkl} 15x.dup.indel.scorer.pkl
        mv ~{insdup_indel_scorer_15x_pkl} 15x.insdup.indel.scorer.pkl
        mv ~{inv_indel_scorer_15x_pkl} 15x.inv.indel.scorer.pkl
        mv ~{del_indel_calibrationScores_15x_hdf5} 15x.del.indel.calibrationScores.hdf5
        mv ~{ins_indel_calibrationScores_15x_hdf5} 15x.ins.indel.calibrationScores.hdf5
        mv ~{dup_indel_calibrationScores_15x_hdf5} 15x.dup.indel.calibrationScores.hdf5
        mv ~{insdup_indel_calibrationScores_15x_hdf5} 15x.insdup.indel.calibrationScores.hdf5
        mv ~{inv_indel_calibrationScores_15x_hdf5} 15x.inv.indel.calibrationScores.hdf5

        mv ~{del_indel_scorer_30x_pkl} 30x.del.indel.scorer.pkl
        mv ~{ins_indel_scorer_30x_pkl} 30x.ins.indel.scorer.pkl
        mv ~{dup_indel_scorer_30x_pkl} 30x.dup.indel.scorer.pkl
        mv ~{insdup_indel_scorer_30x_pkl} 30x.insdup.indel.scorer.pkl
        mv ~{inv_indel_scorer_30x_pkl} 30x.inv.indel.scorer.pkl
        mv ~{del_indel_calibrationScores_30x_hdf5} 30x.del.indel.calibrationScores.hdf5
        mv ~{ins_indel_calibrationScores_30x_hdf5} 30x.ins.indel.calibrationScores.hdf5
        mv ~{dup_indel_calibrationScores_30x_hdf5} 30x.dup.indel.calibrationScores.hdf5
        mv ~{insdup_indel_calibrationScores_30x_hdf5} 30x.insdup.indel.calibrationScores.hdf5
        mv ~{inv_indel_calibrationScores_30x_hdf5} 30x.inv.indel.calibrationScores.hdf5

        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            # Skipping the sample if it has already been processed
            TEST=$( gsutil ls ~{remote_outdir_all}/${SAMPLE_ID}.done || echo "0" )
            if [ "${TEST}" != "0" ]; then
                continue
            fi
            LocalizeSample ${SAMPLE_ID} ~{remote_indir}

            # Removing symbolic <INS>
            if [ ~{remove_symbolic_ins} -eq 1 ]; then
                RemoveSymbolicIns ${SAMPLE_ID}
            fi

            # Adding GT_COUNT
            if [ ~{annotations_have_gt_count} -eq 1 ]; then
                AddGtCount ${SAMPLE_ID} del
                AddGtCount ${SAMPLE_ID} ins
                AddGtCount ${SAMPLE_ID} dup
                AddGtCount ${SAMPLE_ID} insdup
                AddGtCount ${SAMPLE_ID} inv
            fi
            
            # Filtering
            awk -F ',' -v sample="${SAMPLE_ID}" '$1 == sample { print $2 }' ~{sample_coverages_csv} > ${SAMPLE_ID}_coverage.txt
            N_ROWS=$(wc -l < ${SAMPLE_ID}_coverage.txt)
            if [ ${N_ROWS} -ne 1 ]; then
                echo "Expected exactly one row for sample ${SAMPLE_ID}, found ${N_ROWS}" 1>&2
                exit 1
            fi
            COVERAGE=$(head -n 1 ${SAMPLE_ID}_coverage.txt)
            if [ $(echo "${COVERAGE} < 22.5" | bc) -eq 1 ]; then
                PREFIX="15x"
            else
                PREFIX="30x"
            fi
            ScoreAndFilter ${SAMPLE_ID} ${PREFIX} del ${RAM_PER_THREAD_MB}
            ScoreAndFilter ${SAMPLE_ID} ${PREFIX} ins ${RAM_PER_THREAD_MB}
            ScoreAndFilter ${SAMPLE_ID} ${PREFIX} dup ${RAM_PER_THREAD_MB}
            ScoreAndFilter ${SAMPLE_ID} ${PREFIX} insdup ${RAM_PER_THREAD_MB}
            ScoreAndFilter ${SAMPLE_ID} ${PREFIX} inv ${RAM_PER_THREAD_MB}

            # Assembling a single VCF.
            # Adding SVLEN to symbolic ALTs, to avoid overcollapse in `bcftools 
            # merge` downstream.
            for SUFFIX in lenient stringent all ; do
                bcftools concat --threads ${N_THREADS} --allow-overlaps --remove-duplicates --output-type v ${SAMPLE_ID}_del_${SUFFIX}.bcf ${SAMPLE_ID}_ins_${SUFFIX}.bcf ${SAMPLE_ID}_insdup_${SUFFIX}.bcf ${SAMPLE_ID}_dup_${SUFFIX}.bcf ${SAMPLE_ID}_inv_${SUFFIX}.bcf --output ${SAMPLE_ID}_${SUFFIX}.vcf
                rm -f ${SAMPLE_ID}_*_${SUFFIX}.bcf*
                java AddSvlenToSymbolicAlt ${SAMPLE_ID}_${SUFFIX}.vcf > ${SAMPLE_ID}_${SUFFIX}_out.vcf
                rm -f ${SAMPLE_ID}_${SUFFIX}.vcf ; mv ${SAMPLE_ID}_${SUFFIX}_out.vcf ${SAMPLE_ID}_${SUFFIX}_in.vcf
                bcftools sort --max-mem ${RAM_PER_THREAD_MB}M --output-type b ${SAMPLE_ID}_${SUFFIX}_in.vcf --output ${SAMPLE_ID}_${SUFFIX}.bcf
                rm -f ${SAMPLE_ID}_${SUFFIX}_in.vcf ; bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_${SUFFIX}.bcf
            done

            # Uploading
            gsutil mv ./${SAMPLE_ID}_lenient.'bcf*' ~{remote_outdir_lenient}/
            gsutil mv ./${SAMPLE_ID}_stringent.'bcf*' ~{remote_outdir_stringent}/
            gsutil mv ./${SAMPLE_ID}_all.'bcf*' ./${SAMPLE_ID}_'*_xgboost.csv' ~{remote_outdir_all}/
            touch ${SAMPLE_ID}.done
            gsutil mv ${SAMPLE_ID}.done ~{remote_outdir_all}/ && echo 0 || echo 1
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done 3< chunk.csv
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
