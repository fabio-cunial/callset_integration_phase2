version 1.0


# Runs `FilterIntrasampleDevPhase2.wdl` up to task `IdentifyTrainingSites`
# (included) in the same VM for multiple samples.
#
workflow Workpackage2 {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir
        
        File training_resource_vcf_gz
        File training_resource_tbi
        File training_resource_bed

        String identify_training_sites_extra_args = "--sizemin 20 --sizemax 1000000 --sizefilt 20 --pctsize 0.9 --pctseq 0.9 --pick multi"
        
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 50
    }
    parameter_meta {
    }
    
    call Workpackage2Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            training_resource_vcf_gz = training_resource_vcf_gz,
            training_resource_tbi = training_resource_tbi,
            training_resource_bed = training_resource_bed,
            identify_training_sites_extra_args = identify_training_sites_extra_args,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}


#
task Workpackage2Impl {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir
        
        File training_resource_vcf_gz
        File training_resource_tbi
        File training_resource_bed

        String identify_training_sites_extra_args
        
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local REMOTE_DIR=$2
            
            while : ; do
                TEST=$(gsutil -m cp ${REMOTE_DIR}/${SAMPLE_ID}_kanpig.vcf.'gz*' . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${REMOTE_DIR}/${SAMPLE_ID}_kanpig.vcf.gz>. Trying again..."
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
        
        
        # Adds the following INFO fields to the kanpig VCF (using corresponding
        # FORMAT fields):
        #
        # KS_1, KS_2, SQ, GQ, DP, AD_NON_ALT, AD_ALL
        #
        # Remark: the kanpig VCF is assumed to already contain the following
        # INFO fields:
        #
        # SUPP_PAV, SUPP_SNIFFLES, SUPP_PBSV
        #
        function PreprocessVCF() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            # Creating new header lines
            touch ${SAMPLE_ID}_format.hdr.txt
            for format_annot in SQ GQ DP
            do
                bcftools view -h ${INPUT_VCF_GZ} | grep ID=$format_annot, | sed -e 's/FORMAT/INFO/g' >> ${SAMPLE_ID}_format.hdr.txt
            done
            echo '##INFO=<ID=AD_NON_ALT,Number=1,Type=Integer,Description="Coverage for non-alternate alleles">' >> ${SAMPLE_ID}_format.hdr.txt
            echo '##INFO=<ID=AD_ALL,Number=1,Type=Integer,Description="Coverage for all alleles">' >> ${SAMPLE_ID}_format.hdr.txt
            echo '##INFO=<ID=KS_1,Number=1,Type=Integer,Description="Kanpig score 1">' >> ${SAMPLE_ID}_format.hdr.txt
            echo '##INFO=<ID=KS_2,Number=1,Type=Integer,Description="Kanpig score 2">' >> ${SAMPLE_ID}_format.hdr.txt
            echo '##INFO=<ID=GT_COUNT,Number=1,Type=Integer,Description="GT converted into an int in {0,1,2}.">' >> ${SAMPLE_ID}_format.hdr.txt

            # Ensuring that every record has a unique ID. We need to join by
            # CHROM,POS,ID in what follows, since using CHROM,POS,REF,ALT makes
            # `bcftools annotate` segfault.
            bcftools view --header-only ${INPUT_VCF_GZ} > ${SAMPLE_ID}_tmp.vcf
            bcftools view --no-header ${INPUT_VCF_GZ} | awk 'BEGIN { i=0; } { printf("%s\t%s\t%d-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,++i,$3,$4,$5,$6,$7,$8,$9,$10); }' >> ${SAMPLE_ID}_tmp.vcf
            bgzip --compress-level 1 ${SAMPLE_ID}_tmp.vcf
            tabix -f ${SAMPLE_ID}_tmp.vcf.gz
            bcftools view --no-header ${SAMPLE_ID}_tmp.vcf.gz | head -n 10 || echo "0"
            bcftools query -f '%CHROM\t%POS\t%ID\t[%KS]\t[%SQ]\t[%GQ]\t[%DP]\t[%AD]\t[%GT]\t%INFO/SUPP_PBSV\t%INFO/SUPP_SNIFFLES\t%INFO/SUPP_PAV\n' ${SAMPLE_ID}_tmp.vcf.gz | awk '{ \
                KS_1=-1; KS_2=-1; \
                p=0; \
                for (i=1; i<=length($4); i++) { \
                    if (substr($4,i,1)==",") { p=i; break; } \
                } \
                if (p==0) { KS_1=$4; KS_2=$4; } \
                else { KS_1=substr($4,1,p-1); KS_2=substr($4,p+1); } \
                if (KS_1==".") KS_1=-1; \
                if (KS_2==".") KS_2=-1; \
                \
                SQ=$5; \
                if (SQ==".") SQ=-1; \
                \
                GQ=$6; \
                if (GQ==".") GQ=-1; \
                \
                DP=$7; \
                if (DP==".") DP=-1; \
                \
                AD_NON_ALT=-1; AD_ALL=1; \
                p=0; \
                for (i=1; i<=length($8); i++) { \
                    if (substr($8,i,1)==",") { p=i; break; } \
                } \
                if (p==0) { AD_NON_ALT=$8; AD_ALL=$8; } \
                else { AD_NON_ALT=substr($8,1,p-1); AD_ALL=substr($8,p+1); } \
                if (AD_NON_ALT==".") AD_NON_ALT=-1; \
                if (AD_ALL==".") AD_ALL=-1; \
                \
                GT_COUNT=-1; \
                if ($9=="0/0" || $9=="0|0" || $9=="./."  || $9==".|." || $9=="./0" || $9==".|0" || $9=="0/." || $9=="0|.") GT_COUNT=0; \
                else if ($9=="0/1" || $9=="0|1" || $9=="1/0" || $9=="1|0" || $9=="./1" || $9==".|1" || $9=="1/." || $9=="1|.") GT_COUNT=1; \
                else if ($9=="1/1" || $9=="1|1") GT_COUNT=2; \
                \
                printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$1,$2,$3,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,$10,$11,$12); \
            }' | bgzip -c > ${SAMPLE_ID}_format.tsv.gz
            tabix -s1 -b2 -e2 ${SAMPLE_ID}_format.tsv.gz
            bcftools annotate --threads ${N_THREADS} -a ${SAMPLE_ID}_format.tsv.gz -h ${SAMPLE_ID}_format.hdr.txt -c CHROM,POS,ID,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,SUPP_PBSV,SUPP_SNIFFLES,SUPP_PAV ${SAMPLE_ID}_tmp.vcf.gz -Oz -o ${SAMPLE_ID}_preprocessed.vcf.gz
            bcftools view --no-header ${SAMPLE_ID}_preprocessed.vcf.gz | head -n 10 || echo "0"
            bcftools index -t ${SAMPLE_ID}_preprocessed.vcf.gz
            # TODO do we still need to hard-filter SVLEN >= 50 for extract, and
            # should we clear existing filters?
            
            # Removing temporary files
            rm -f ${SAMPLE_ID}_format.hdr.txt ${SAMPLE_ID}_tmp.vcf* ${SAMPLE_ID}_format.tsv.gz 
        }
        
        
        function IdentifyTrainingSites() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            ${TIME_COMMAND} truvari bench -b ~{training_resource_vcf_gz} -c ${INPUT_VCF_GZ} --includebed ~{training_resource_bed} ~{identify_training_sites_extra_args} -o ${SAMPLE_ID}_truvari/
            mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_training.vcf.gz
            mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz.tbi ${SAMPLE_ID}_training.vcf.gz.tbi
            
            # Removing temporary files
            rm -rf ./${SAMPLE_ID}_truvari/
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        source activate truvari5
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            LocalizeSample ${SAMPLE_ID} ~{remote_indir}
            
            PreprocessVCF ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz
            IdentifyTrainingSites ${SAMPLE_ID} ${SAMPLE_ID}_preprocessed.vcf.gz
            
            gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ${SAMPLE_ID}_preprocessed.vcf.gz ~{remote_outdir}
            gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ${SAMPLE_ID}_preprocessed.vcf.gz.tbi ~{remote_outdir}
            gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ${SAMPLE_ID}_training.vcf.gz ~{remote_outdir}
            gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv ${SAMPLE_ID}_training.vcf.gz.tbi ~{remote_outdir}
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
