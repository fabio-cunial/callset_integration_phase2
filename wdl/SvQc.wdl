version 1.0


# 
#
workflow SvQc {
    input {
        File chunk_csv
        File tr_bed
        File reference_fai
        String remote_outdir

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        chunk_csv: "Format: ID,COVERAGE,PAV_TBI,PAV_VCF_GZ,PBSV_TBI,PBSV_VCF_GZ,SNIFFLES_TBI,SNIFFLES_VCF_GZ"
        tr_bed: "Assumed to contain at least one record for each standard chromosome"
    }

    call Impl {
        input:
            chunk_csv = chunk_csv,
            tr_bed = tr_bed,
            reference_fai = reference_fai,
            remote_outdir = remote_outdir,
            docker_image = docker_image
    }
    
    output {
    }
}


# ~3m per sample, including localization.
#
task Impl {
    input {
        File chunk_csv
        File tr_bed
        File reference_fai
        String remote_outdir

        String docker_image
        Int n_cpu = 1
        Int ram_size_gb = 4
        Int disk_size_gb = 10
        Int preemptible_number = 4
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




        # ----------------------- Steps of the pipeline ------------------------

        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2
            
            local COVERAGE=$(echo ${LINE} | cut -d , -f 2)
            local PAV_TBI=$(echo ${LINE} | cut -d , -f 3)
            local PAV_VCF_GZ=$(echo ${LINE} | cut -d , -f 4)
            local PBSV_TBI=$(echo ${LINE} | cut -d , -f 5)
            local PBSV_VCF_GZ=$(echo ${LINE} | cut -d , -f 6)
            local SNIFFLES_TBI=$(echo ${LINE} | cut -d , -f 7)
            local SNIFFLES_VCF_GZ=$(echo ${LINE} | cut -d , -f 8)
            
            gcloud storage cp ${PAV_VCF_GZ} ./${SAMPLE_ID}_tmp.vcf.gz
            ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="SNV"' --output-type v ${SAMPLE_ID}_tmp.vcf.gz --output ${SAMPLE_ID}_pav.vcf
            rm -f ${SAMPLE_ID}_tmp.vcf.gz
            gcloud storage cp ${PBSV_VCF_GZ} ./${SAMPLE_ID}_pbsv.vcf.gz
            gunzip ${SAMPLE_ID}_pbsv.vcf.gz
            gcloud storage cp ${SNIFFLES_VCF_GZ} ./${SAMPLE_ID}_sniffles.vcf.gz
            gunzip ${SAMPLE_ID}_sniffles.vcf.gz
        }


        # Puts in canonical form a raw VCF from an SV caller.
        #
        function CanonizeVcf() {
            local SAMPLE_ID=$1
            local CALLER_ID=$2

            mv ${SAMPLE_ID}_${CALLER_ID}.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # Ensuring that SVLEN has the correct type for bcftools norm
            bcftools view --header-only ${SAMPLE_ID}_${CALLER_ID}_in.vcf | sed 's/ID=SVLEN,Number=.,/ID=SVLEN,Number=A,/g' > ${SAMPLE_ID}_${CALLER_ID}_header.txt
            ${TIME_COMMAND} bcftools reheader --header ${SAMPLE_ID}_${CALLER_ID}_header.txt --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics -any --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # Removing SNVs and small INDELs to speed up downstream steps
            if [ ${CALLER_ID} = 'pav' ]; then
                ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="SNV" || ABS(SVLEN)<20' --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            fi

            # Making sure SVLEN and SVTYPE are consistently annotated
            ${TIME_COMMAND} java -cp ~{docker_dir} AddSvtypeSvlen ${SAMPLE_ID}_${CALLER_ID}_in.vcf > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${SAMPLE_ID}_${CALLER_ID}.vcf
        }
        
        
        # Deletes all and only the files downloaded by `LocalizeSample()`.
        #
        function DelocalizeSample() {
            local SAMPLE_ID=$1
            
            rm -f ${SAMPLE_ID}_pav.vcf* ${SAMPLE_ID}_pbsv.vcf* ${SAMPLE_ID}_sniffles.vcf*
        }


        function Count() {
            local SAMPLE_ID=$1
            local MIN_SV_LENGTH=$2
            local SVTYPE=$3
            local REGIONS_MODE=$4
            
            if [ ${REGIONS_MODE} = "0" ]; then
                local REGIONS_ARG=" "
            elif [ ${REGIONS_MODE} = "1" ]; then
                local REGIONS_ARG="--targets-file tr.bed --targets-overlap pos"
            elif [ ${REGIONS_MODE} = "2" ]; then
                local REGIONS_ARG="--targets-file not_tr.bed --targets-overlap pos"
            fi
            if [ ${SVTYPE} = "BND" -o ${SVTYPE} = "UNK" -o ${SVTYPE} = "SUB" ]; then
                local N_RECORDS_PAV=$(bcftools query --format '%ID\n' --include "SVTYPE='${SVTYPE}'" ${REGIONS_ARG} ${SAMPLE_ID}_pav.vcf | wc -l)
                local N_RECORDS_PBSV=$(bcftools query --format '%ID\n' --include "SVTYPE='${SVTYPE}'" ${REGIONS_ARG} ${SAMPLE_ID}_pbsv.vcf | wc -l)
                local N_RECORDS_SNIFFLES=$(bcftools query --format '%ID\n' --include "SVTYPE='${SVTYPE}'" ${REGIONS_ARG} ${SAMPLE_ID}_sniffles.vcf | wc -l)
            else
                local N_RECORDS_PAV=$(bcftools query --format '%ID\n' --include "SVTYPE='${SVTYPE}' && ABS(SVLEN)>=${MIN_SV_LENGTH}" ${REGIONS_ARG} ${SAMPLE_ID}_pav.vcf | wc -l)
                local N_RECORDS_PBSV=$(bcftools query --format '%ID\n' --include "SVTYPE='${SVTYPE}' && ABS(SVLEN)>=${MIN_SV_LENGTH}" ${REGIONS_ARG} ${SAMPLE_ID}_pbsv.vcf | wc -l)
                local N_RECORDS_SNIFFLES=$(bcftools query --format '%ID\n' --include "SVTYPE='${SVTYPE}' && ABS(SVLEN)>=${MIN_SV_LENGTH}" ${REGIONS_ARG} ${SAMPLE_ID}_sniffles.vcf | wc -l)
            fi
            COUNTS="${COUNTS},${N_RECORDS_PAV},${N_RECORDS_PBSV},${N_RECORDS_SNIFFLES}"
        }




        # ---------------------------- Main program ----------------------------

        # Initializing BED files
        mv ~{tr_bed} ./tr.bed
        bedtools complement -L -i tr.bed -g ~{reference_fai} > not_tr.bed

        # Counting
        SV_TYPES_PRIMARY="DEL INS DUP INV"
        SV_TYPES_SECONDARY="BND SUB UNK"
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            COVERAGE=$(echo ${LINE} | cut -d , -f 2)

            # Skipping the sample if it has already been processed
            TEST=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}.counts || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi

            # Putting the raw VCFs in canonical form
            LocalizeSample ${SAMPLE_ID} ${LINE}
            CanonizeVcf ${SAMPLE_ID} "pbsv"
            CanonizeVcf ${SAMPLE_ID} "sniffles"
            CanonizeVcf ${SAMPLE_ID} "pav"

            # Counting
            COUNTS="${SAMPLE_ID},${COVERAGE}"
            for SVLEN in 20 50; do
                for REGIONS in 0 1 2; do
                    for SVTYPE in ${SV_TYPES_PRIMARY}; do
                        Count ${SAMPLE_ID} ${SVLEN} ${SVTYPE} ${REGIONS}
                    done
                done
            done
            for REGIONS in 0 1 2; do
                for SVTYPE in ${SV_TYPES_SECONDARY}; do
                    Count ${SAMPLE_ID} 0 ${SVTYPE} ${REGIONS}
                done
            done
            echo "${COUNTS}" > ${SAMPLE_ID}.counts

            # Uploading
            gcloud storage cp ${SAMPLE_ID}.counts ~{remote_outdir}/${SAMPLE_ID}.counts
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done 3< ~{chunk_csv}
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
