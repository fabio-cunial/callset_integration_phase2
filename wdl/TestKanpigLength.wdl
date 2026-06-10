version 1.0


# Measures the effect of re-genotyping on true SVs of different lengths.
#
workflow TestKanpigLength {
    input {
        String sample_id
        String sex
        
        File dipcall_vcf_gz
        File dipcall_bed

        File alignments_bam
        File alignments_bai

        String kanpig_params_singlesample = "--neighdist 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"

        File ploidy_bed_male
        File ploidy_bed_female
        File reference_fa
        File reference_fai
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        kanpig_params_singlesample: "The same used in the main AoU pipeline."
    }
    
    call Impl {
        input:
            sample_id = sample_id,
            sex = sex,

            dipcall_vcf_gz = dipcall_vcf_gz,
            dipcall_bed = dipcall_bed,

            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,

            kanpig_params_singlesample = kanpig_params_singlesample,

            ploidy_bed_male = ploidy_bed_male,
            ploidy_bed_female = ploidy_bed_female,
            reference_fa = reference_fa,
            reference_fai = reference_fai,

            docker_image = docker_image
    }
    
    output {
        File counts_csv = Impl.counts_csv
    }
}


#
task Impl {
    input {
        String sample_id
        String sex
        
        File dipcall_vcf_gz
        File dipcall_bed

        File alignments_bam
        File alignments_bai

        String kanpig_params_singlesample = "--neighdist 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"

        File ploidy_bed_male
        File ploidy_bed_female
        File reference_fa
        File reference_fai
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 4
        Int disk_size_gb = 100
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
        



        # ---------------------- Steps of the pipeline -------------------------

        # Puts in canonical form a raw VCF from dipcall.
        #
        function CanonizeDipcallVcf() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local MIN_SV_LENGTH=$3
            local MAX_SV_LENGTH=$4
            
            
            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_in.vcf.gz
            tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only records in the dipcall BED
            ${TIME_COMMAND} bcftools filter --regions-file ${SAMPLE_ID}.bed --regions-overlap pos --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics -any --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Removing SNVs, replacement records, records that are not marked
            # as ALT, records with a FILTER, and records with unresolved
            # REF/ALT.
            ${TIME_COMMAND} bcftools filter --exclude '(STRLEN(REF)=1 && STRLEN(ALT)=1) || (STRLEN(REF)>1 && STRLEN(ALT)>1) || GT!="alt" || (FILTER!="PASS" && FILTER!=".") || REF="*" || ALT="*"' --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Making sure SVLEN and SVTYPE are consistently annotated        
            ${TIME_COMMAND} java -cp ~{docker_dir} AddSvtypeSvlen ${SAMPLE_ID}_in.vcf.gz | bgzip > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            # Keeping only INS and DEL in the given length range
            ${TIME_COMMAND} bcftools filter --include '(SVTYPE="INS" || SVTYPE="DEL") && ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type z ${SAMPLE_ID}_in.vcf.gz > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; tabix -@ ${N_THREADS} -f ${SAMPLE_ID}_in.vcf.gz
            
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_canonized.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_canonized.vcf.gz.tbi
        }




        # ------------------------- Main program -------------------------------

        INFINITY="1000000000"
        ~{docker_dir}/kanpig --version 1>&2

        mv ~{dipcall_bed} ~{sample_id}.bed
        CanonizeDipcallVcf ~{sample_id} ~{dipcall_vcf_gz} 50 ${INFINITY}

        # Backing up the original GT in INFO
        bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT]\n' ~{sample_id}_canonized.vcf.gz | bgzip -c > ~{sample_id}_annotations.tsv.gz
        tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ~{sample_id}_annotations.tsv.gz
        echo '##INFO=<ID=ORIGINAL_GT,Number=1,Type=String,Description="Original GT.">' > ~{sample_id}_header.txt
        COLUMNS='CHROM,POS,~ID,REF,ALT,INFO/ORIGINAL_GT'
        ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ~{sample_id}_annotations.tsv.gz --header-lines ~{sample_id}_header.txt --columns ${COLUMNS} --output-type z ~{sample_id}_canonized.vcf.gz --output ~{sample_id}_annotated.vcf.gz
        tabix -@ ${N_THREADS} -f ~{sample_id}_annotated.vcf.gz
        rm -f ~{sample_id}_annotations.tsv.gz ~{sample_id}_header.txt

        # Re-genotyping
        if [ ~{sex} == "M" ]; then
            PLOIDY_BED=$(echo ~{ploidy_bed_male})
        else
            PLOIDY_BED=$(echo ~{ploidy_bed_female})
        fi
        ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_singlesample} --sizemin 10 --sizemax ${INFINITY} --reference ~{reference_fa} --input ~{sample_id}_annotated.vcf.gz --reads ~{alignments_bam} --out ~{sample_id}_kanpig.vcf
        bcftools query --format '%INFO/SVTYPE,%INFO/SVLEN,%INFO/ORIGINAL_GT,[%GT]\n' ~{sample_id}_kanpig.vcf > ~{sample_id}_counts.csv
    >>>
    
    output {
        File counts_csv = sample_id + "_counts.csv"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
