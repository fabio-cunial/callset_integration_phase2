version 1.0


# Re-genotypes an inter-sample VCF using the reads of every sample.
#
workflow Workpackage9 {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir
        
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
    
    call Workpackage9Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb,
            kanpig_params_multisample = kanpig_params_multisample
    }
    
    output {
    }
}


# Performance on 10'070 samples, 15x, GRCh38:
#
# CAL_SENS  CPU     RAM     TIME
# <=0.7     500%    16G     30m
#
task Workpackage9Impl {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir
        
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
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int compression_level = 1
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
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
        }
        
        
        # Deletes all the files related to a sample
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}_*
        }
        
        
        function Kanpig() {
            local SAMPLE_ID=$1
            local SEX=$2
            local ALIGNMENTS_BAM=$3

            # Re-genotyping
            touch ${ALIGNMENTS_BAM}.bai
            if [ ${SEX} == "M" ]; then
                PLOIDY_BED=$(echo ~{ploidy_bed_male})
            else
                PLOIDY_BED=$(echo ~{ploidy_bed_female})
            fi
            export RUST_BACKTRACE="full"
            ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_multisample} --reference ~{reference_fa} --input truvari_collapsed_for_kanpig.vcf.gz --reads ${ALIGNMENTS_BAM} --out ${SAMPLE_ID}_kanpig.vcf

            # Building the GT-only output file
            echo "ID\t${SAMPLE_ID}" > ${SAMPLE_ID}_gts.txt
            bcftools view --no-header ${SAMPLE_ID}_kanpig.vcf | cut -f 3,10 >> ${SAMPLE_ID}_gts.txt
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        while : ; do
            TEST=$(gsutil -m cp ~{remote_indir}/truvari_collapsed_for_kanpig.vcf.'gz*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading the VCF to genotype. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            SEX=$(echo ${LINE} | cut -d , -f 2)
            
            TEST=$(gsutil ls ~{remote_outdir}/${SAMPLE_ID}_gts.txt && echo "" || echo "ERROR")
            if [ ${TEST} = "ERROR" ]; then
                # Proceeding only if genotypes have not already been computed
                LocalizeSample ${SAMPLE_ID} ${LINE}
                Kanpig ${SAMPLE_ID} ${SEX} ${SAMPLE_ID}_aligned.bam
                while : ; do
                    TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${SAMPLE_ID}_gts.txt ~{remote_outdir}/ && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error uploading the GT file. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                DelocalizeSample ${SAMPLE_ID}
                ls -laht
            fi
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
