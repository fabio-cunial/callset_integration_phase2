version 1.0


# Merges a cohort VCF and a single-sample VCF before re-genotyping using the
# reads of every sample.
#
# We work on 5 HPRC samples: HG03579 NA18906 NA19240 NA20129 NA21309
#
workflow PersonalizedCohortVcf {
    input {
        File sv_integration_chunk_tsv
        Int n_rows
        String filter_string = "FORMAT/CALIBRATION_SENSITIVITY<=0.7"
        
        File cohort_vcf_gz
        String remote_outdir
        Int only_50_bp
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        
        Int n_cpu = 6
        Int ram_size_gb = 20
        Int disk_size_gb = 100
        String kanpig_params_multisample = "--sizemin 20 --sizemax 10000 --neighdist 500 --gpenalty 0.04 --hapsim 0.97"
    }
    parameter_meta {
        sv_integration_chunk_tsv: "A subset of the rows of table `sv_integration_hg38`, without the header."
    }
    
    scatter (i in range(n_rows)) {
        call Impl {
            input:
                sv_integration_chunk_tsv = sv_integration_chunk_tsv,
                selected_row = i+1,
                filter_string = filter_string,
                
                cohort_vcf_gz = cohort_vcf_gz,
                remote_outdir = remote_outdir,
                only_50_bp = only_50_bp,
                
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                ploidy_bed_female = ploidy_bed_female,
                ploidy_bed_male = ploidy_bed_male,
                
                n_cpu = n_cpu,
                ram_size_gb = ram_size_gb,
                disk_size_gb = disk_size_gb,
                kanpig_params_multisample = kanpig_params_multisample
        }
    }
    
    output {
    }
}


# Performance on 10'070 samples, 15x, GRCh38:
#
# CAL_SENS  CPU     RAM     TIME
# <=0.7     500%    16G     30m
#
task Impl {
    input {
        File sv_integration_chunk_tsv
        Int selected_row
        String filter_string
        
        File cohort_vcf_gz
        String remote_outdir
        Int only_50_bp
        
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
        
        Int n_cpu
        Int ram_size_gb
        Int disk_size_gb
        String kanpig_params_multisample
    }
    parameter_meta {
        sv_integration_chunk_tsv: "The single-sample VCFs in this file are assumed to have been scored by xgboost but not to have been filtered based on such scores."
        filter_string: "The same single-sample filters used to build `cohort_vcf_gz`."
    }
    
    String docker_dir = "/callset_integration"
    Int compression_level = 1
    
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
            local LINE=$2
            
            ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 3)
            ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 4)
            SAMPLE_VCF_GZ=$(echo ${LINE} | cut -d , -f 5)
            SAMPLE_TBI=$(echo ${LINE} | cut -d , -f 6)
            while : ; do
                TEST=$(gsutil -m cp ${ALIGNED_BAM} ./${SAMPLE_ID}_aligned.bam && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${ALIGNED_BAM}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gsutil -m cp ${ALIGNED_BAI} ./${SAMPLE_ID}_aligned.bam.bai && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${ALIGNED_BAI}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gsutil -m cp ${SAMPLE_VCF_GZ} ./${SAMPLE_ID}.vcf.gz && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${SAMPLE_VCF_GZ}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gsutil -m cp ${SAMPLE_TBI} ./${SAMPLE_ID}.vcf.gz.tbi && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${SAMPLE_TBI}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        }
        
        
        # Deletes all the files related to a sample
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}_* ${SAMPLE_ID}.*
        }
        
        
        # Merges a sample VCF with the cohort VCF.
        #
        # Remark: `bcftools concat` and `bcftools merge` take approx. the same
        # amount of time.
        #
        function MergeVcfs() {
            local SAMPLE_ID=$1
            local SAMPLE_VCF_GZ=$2
            local COHORT_VCF_GZ=$3
            
            
            if [ ~{only_50_bp} -ne 0 ]; then
                ${TIME_COMMAND} truvari anno svinfo -m 1 ${SAMPLE_VCF_GZ} | bgzip > ${SAMPLE_ID}_annotated.vcf.gz
                tabix -f ${SAMPLE_ID}_annotated.vcf.gz
                rm -f ${SAMPLE_VCF_GZ}*
                bcftools filter --include 'SVLEN>=50 || SVLEN<=-50' --output-type z ${SAMPLE_ID}_annotated.vcf.gz > ${SAMPLE_ID}_filtered.vcf.gz
                rm -f ${SAMPLE_ID}_annotated.vcf.gz*
                mv ${SAMPLE_ID}_filtered.vcf.gz ${SAMPLE_VCF_GZ}
                tabix -f ${SAMPLE_VCF_GZ}
            fi
            
            FILTER_STRING="~{filter_string}"
            INCLUDE_STR="--include ${FILTER_STRING}"
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} ${INCLUDE_STR} --output-type z ${SAMPLE_VCF_GZ} > ${SAMPLE_ID}_filtered.vcf.gz
            tabix -f ${SAMPLE_ID}_filtered.vcf.gz
            rm -f ${SAMPLE_VCF_GZ}*
            
            ${TIME_COMMAND} bcftools reheader --threads ${N_THREADS} --samples samples.txt ${SAMPLE_ID}_filtered.vcf.gz > ${SAMPLE_ID}_reheader.vcf.gz
            tabix -f ${SAMPLE_ID}_reheader.vcf.gz
            rm -f ${SAMPLE_ID}_filtered.vcf.gz*
            
            ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --rm-dups exact --output-type z ${COHORT_VCF_GZ} ${SAMPLE_ID}_reheader.vcf.gz > ${SAMPLE_ID}_merged.vcf.gz
            tabix -f ${SAMPLE_ID}_merged.vcf.gz
            rm -f ${SAMPLE_ID}_reheader.vcf.gz*
        }
        
        
        function Kanpig() {
            local SAMPLE_ID=$1
            local SEX=$2
            local ALIGNMENTS_BAM=$3
            local INPUT_VCF_GZ=$4

            touch ${ALIGNMENTS_BAM}.bai
            if [ ${SEX} == "M" ]; then
                PLOIDY_BED=$(echo ~{ploidy_bed_male})
            else
                PLOIDY_BED=$(echo ~{ploidy_bed_female})
            fi
            export RUST_BACKTRACE="full"
            ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_multisample} --reference ~{reference_fa} --input ${INPUT_VCF_GZ} --reads ${ALIGNMENTS_BAM} --out ${SAMPLE_ID}_kanpig.vcf
            ${TIME_COMMAND} bgzip ${SAMPLE_ID}_kanpig.vcf
            tabix -f ${SAMPLE_ID}_kanpig.vcf.gz
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        # Preprocessing the cohort VCF
        echo "SAMPLE" > samples.txt
        ${TIME_COMMAND} bcftools reheader --threads ${N_THREADS} --samples samples.txt ~{cohort_vcf_gz} > cohort.vcf.gz
        tabix -f cohort.vcf.gz
        rm -f ~{cohort_vcf_gz}
        if [ ~{only_50_bp} -ne 0 ]; then
            ${TIME_COMMAND} truvari anno svinfo -m 1 cohort.vcf.gz | bgzip > cohort_annotated.vcf.gz
            rm -f cohort.vcf.gz*
            mv cohort_annotated.vcf.gz cohort.vcf.gz
            tabix -f cohort.vcf.gz
            bcftools filter --include 'SVLEN>=50 || SVLEN<=-50' --output-type z cohort.vcf.gz > tmp.vcf.gz
            rm -f cohort.vcf.gz*
            mv tmp.vcf.gz cohort.vcf.gz
            tabix -f cohort.vcf.gz
        fi
        
        # Re-genotyping
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' | head -n ~{selected_row} | tail -n 1 > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            SEX=$(echo ${LINE} | cut -d , -f 2)

            LocalizeSample ${SAMPLE_ID} ${LINE}
            MergeVcfs ${SAMPLE_ID} ${SAMPLE_ID}.vcf.gz cohort.vcf.gz
            Kanpig ${SAMPLE_ID} ${SEX} ${SAMPLE_ID}_aligned.bam ${SAMPLE_ID}_merged.vcf.gz
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${SAMPLE_ID}_kanpig.vcf.gz ~{remote_outdir}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading the GT file. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done < chunk.csv
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_squish"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
