version 1.0


# Annotates the INFO field of every intra-sample VCF with FORMAT fields
# annotated by kanpig, and marks true-positive records using a training
# resource.
#
workflow SV_Integration_Workpackage2 {
    input {
        File sv_integration_chunk_tsv
        String remote_indir
        String remote_outdir
        
        File training_resource_vcf_gz
        File training_resource_tbi
        File training_resource_bed
    }
    parameter_meta {
        sv_integration_chunk_tsv: "We assume that every intra-sample-merged VCF has already been subset to the correct length range upstream."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
        training_resource_vcf_gz: "We assume that the training resource VCF has already been subset to the correct length range upstream."
        training_resource_bed: "Training resource calls can belong only to these regions. Typically a high-confidence dipcall BED, or a BED derived from dipcall's."
    }
    
    call Impl {
        input:
            sv_integration_chunk_tsv = sv_integration_chunk_tsv,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir,
            training_resource_vcf_gz = training_resource_vcf_gz,
            training_resource_tbi = training_resource_tbi,
            training_resource_bed = training_resource_bed
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
        
        File training_resource_vcf_gz
        File training_resource_tbi
        File training_resource_bed
        
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 20
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
        
        
        # Adds the following INFO fields to the kanpig VCF, copying their values
        # from corresponding FORMAT fields:
        #
        # KS_1, KS_2, SQ, GQ, DP, AD_NON_ALT, AD_ALL
        #
        # Remark: the kanpig VCF is assumed to already contain the following
        # INFO fields:
        #
        # SUPP_PAV, SUPP_SNIFFLES, SUPP_PBSV
        #
        # and to be such that every record has a distinct ID.
        #
        function CopyFormatToInfo() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            # Creating new header lines
            touch ${SAMPLE_ID}_header.txt
            for FIELD in SQ GQ DP
            do
                bcftools view --header-only ${INPUT_VCF_GZ} | grep ID="${FIELD}," | sed -e 's/FORMAT/INFO/g' >> ${SAMPLE_ID}_header.txt
            done
            echo '##INFO=<ID=AD_NON_ALT,Number=1,Type=Integer,Description="Coverage for non-alternate alleles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=AD_ALL,Number=1,Type=Integer,Description="Coverage for all alleles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=KS_1,Number=1,Type=Integer,Description="Kanpig score 1">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=KS_2,Number=1,Type=Integer,Description="Kanpig score 2">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=GT_COUNT,Number=1,Type=Integer,Description="GT converted to an integer in {0,1,2}.">' >> ${SAMPLE_ID}_header.txt
            
            # Copying fields from FORMAT to INFO. Every record is assumed to
            # have a distinct ID, which is enforced by the steps of the
            # pipeline upstream.
            bcftools query -f '%CHROM\t%POS\t%ID\t[%KS]\t[%SQ]\t[%GQ]\t[%DP]\t[%AD]\t[%GT]\t%INFO/SUPP_PBSV\t%INFO/SUPP_SNIFFLES\t%INFO/SUPP_PAV\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
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
                if ($9=="0/0" || $9=="0|0" || $9=="./."  || $9==".|." || $9=="./0" || $9==".|0" || $9=="0/." || $9=="0|." || $9=="0" || $9==".") GT_COUNT=0; \
                else if ($9=="0/1" || $9=="0|1" || $9=="1/0" || $9=="1|0" || $9=="./1" || $9==".|1" || $9=="1/." || $9=="1|." || $9=="1") GT_COUNT=1; \
                else if ($9=="1/1" || $9=="1|1") GT_COUNT=2; \
                \
                printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$1,$2,$3,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,$10,$11,$12); \
            }' | bgzip -c > ${SAMPLE_ID}_format.tsv.gz
            tabix -s1 -b2 -e2 ${SAMPLE_ID}_format.tsv.gz
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ${SAMPLE_ID}_format.tsv.gz --header-lines ${SAMPLE_ID}_header.txt --columns CHROM,POS,~ID,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,SUPP_PBSV,SUPP_SNIFFLES,SUPP_PAV --output-type z ${INPUT_VCF_GZ} > ${SAMPLE_ID}_out.vcf.gz
            mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -f ${SAMPLE_ID}_in.vcf.gz
            bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 5 || echo "0"
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_preprocessed.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_preprocessed.vcf.gz.tbi
            
            # Removing temporary files
            rm -f ${SAMPLE_ID}_header.txt ${SAMPLE_ID}_format.tsv.gz 
        }
        
        
        # Every record that has a stringent `truvari bench` match with some
        # records in the resource.
        #
        function GetTrainingRecords() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            
            ${TIME_COMMAND} truvari bench -b ~{training_resource_vcf_gz} -c ${INPUT_VCF_GZ} --includebed ~{training_resource_bed} --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick multi -o ${SAMPLE_ID}_truvari/
            mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz ${SAMPLE_ID}_training.vcf.gz
            mv ${SAMPLE_ID}_truvari/tp-comp.vcf.gz.tbi ${SAMPLE_ID}_training.vcf.gz.tbi
            
            # Removing temporary files
            rm -rf ./${SAMPLE_ID}_truvari/
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            
            # Skipping the sample if it has already been processed
            TEST1=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}_preprocessed.vcf.gz || echo "0" )
            TEST2=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}_training.vcf.gz || echo "0" )
            if [ ${TEST1} != "0" -a ${TEST2} != "0" ]; then
                continue
            fi
            
            # Annotating and marking training records
            LocalizeSample ${SAMPLE_ID} ~{remote_indir}
            CopyFormatToInfo ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz
            GetTrainingRecords ${SAMPLE_ID} ${SAMPLE_ID}_preprocessed.vcf.gz
            
            # Uploading
            rm -f ${SAMPLE_ID}_list.txt
            echo "${SAMPLE_ID}_preprocessed.vcf.gz" >> ${SAMPLE_ID}_list.txt
            echo "${SAMPLE_ID}_preprocessed.vcf.gz.tbi" >> ${SAMPLE_ID}_list.txt 
            echo "${SAMPLE_ID}_training.vcf.gz" >> ${SAMPLE_ID}_list.txt
            echo "${SAMPLE_ID}_training.vcf.gz.tbi" >> ${SAMPLE_ID}_list.txt
            cat ${SAMPLE_ID}_list.txt | gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} mv -I ~{remote_outdir}/
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
