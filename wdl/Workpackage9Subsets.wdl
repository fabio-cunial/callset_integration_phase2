version 1.0


# Re-genotypes an inter-sample VCF using the reads of every sample.
#
workflow Workpackage9Subsets {
    input {
        File sv_integration_chunk_tsv
        File truvari_intersample_vcf_gz
        File truvari_intersample_tbi
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
    
    
    call Subset as present {
        input:
            truvari_intersample_vcf_gz = truvari_intersample_vcf_gz,
            truvari_intersample_tbi = truvari_intersample_tbi,
            only_present = 1,
            only_50bp = 0
    }
    call Subset as present_and_50 {
        input:
            truvari_intersample_vcf_gz = truvari_intersample_vcf_gz,
            truvari_intersample_tbi = truvari_intersample_tbi,
            only_present = 1,
            only_50bp = 1
    }
    call Workpackage9Impl as present_2 {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            intersample_vcf_gz = present.out_vcf_gz,
            intersample_tbi = present.out_tbi,
            intersample_type = "present",
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
    call Workpackage9Impl as present_and_50_2 {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            intersample_vcf_gz = present_and_50.out_vcf_gz,
            intersample_tbi = present_and_50.out_tbi,
            intersample_type = "present_and_50",
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


#
task Subset {
    input {
        File truvari_intersample_vcf_gz
        File truvari_intersample_tbi
        
        Int only_50bp
        Int only_present
        
        Int n_cpu = 8
        Int ram_size_gb = 16
        Int disk_size_gb = 500
    }
    parameter_meta {
    }
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Filtering
        if [ ~{only_50bp} -eq 1 ]; then
            if [ ~{only_present} -eq 1 ]; then
                ${TIME_COMMAND} bcftools filter --include '(SVLEN>=50 || SVLEN<=-50) && COUNT(GT="alt")>0' --output-type z ~{truvari_intersample_vcf_gz} > tmp.vcf.gz
            else
                ${TIME_COMMAND} bcftools filter --include 'SVLEN>=50 || SVLEN<=-50' --output-type z ~{truvari_intersample_vcf_gz} > tmp.vcf.gz
            fi
        else
            if [ ~{only_present} -eq 1 ]; then
                ${TIME_COMMAND} bcftools filter --include 'COUNT(GT="alt")>0' --output-type z ~{truvari_intersample_vcf_gz} > tmp.vcf.gz
            fi
        fi
        ${TIME_COMMAND} tabix -f tmp.vcf.gz
        
        # Preparing for kanpig
        bcftools view --header-only tmp.vcf.gz > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > out.vcf
        echo '##INFO=<ID=ORIGINAL_ID,Number=1,Type=String,Description="Original ID from truvari collapse">' >> out.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> out.vcf
        bcftools view --no-header tmp.vcf.gz | awk 'BEGIN { i=0; } { gsub(/;/,"_",$3); printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s;ORIGINAL_ID=%s\tGT\t0/1\n",$1,$2,++i,$4,$5,$6,$7,$8,$3); }' >> out.vcf
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} out.vcf
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration_phase2_squish"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
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
        File intersample_vcf_gz
        File intersample_tbi
        String intersample_type
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
            ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_multisample} --reference ~{reference_fa} --input ~{intersample_vcf_gz} --reads ${ALIGNMENTS_BAM} --out ${SAMPLE_ID}_kanpig.vcf

            # Building the GT-only output file
            echo "ID\t${SAMPLE_ID}" > ${SAMPLE_ID}_gts_~{intersample_type}.txt
            bcftools view --no-header ${SAMPLE_ID}_kanpig.vcf | cut -f 3,10 >> ${SAMPLE_ID}_gts_~{intersample_type}.txt
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            SEX=$(echo ${LINE} | cut -d , -f 2)
            
            TEST=$(gsutil ls ~{remote_outdir}/${SAMPLE_ID}_gts_~{intersample_type}.txt && echo "" || echo "ERROR")
            if [ ${TEST} = "ERROR" ]; then
                # Proceeding only if genotypes have not already been computed
                LocalizeSample ${SAMPLE_ID} ${LINE}
                Kanpig ${SAMPLE_ID} ${SEX} ${SAMPLE_ID}_aligned.bam
                while : ; do
                    TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${SAMPLE_ID}_gts_~{intersample_type}.txt ~{remote_outdir}/ && echo 0 || echo 1)
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
        docker: "fcunial/callset_integration_phase2_squish"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
