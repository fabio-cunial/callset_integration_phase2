version 1.0


# Concatenates the `truvari collapse` chunks. Ensures that every record in
# output has a globally unique ID (this is necessary for kanpig downstream;
# duplicated IDs may arise naturally from the previous steps of the pipeline),
# and an INFO field that counts the number of samples it was discovered in. The
# original ID is copied to an INFO field.
#
# Partitions the result of the concatenation into the following files:
# - Frequent: the subset of all records that were discovered in at least the
#   specified fraction of all samples. This VCF has no FORMAT and SAMPLE
#   columns.
# - Infrequent: the subset of all records that were discovered in less than the
#   specified fraction of all samples. This VCF has all the original SAMPLE
#   columns.
#
workflow SV_Integration_Workpackage8 {
    input {
        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"]
        String remote_indir
        String remote_workpackages1_dir
        String remote_outdir
        Float n_samples_fraction_frequent = 0.1
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        chromosomes: "The order of the chromosomes becomes their order in the output VCF."
        remote_indir: "Without final slash"
        remote_workpackages1_dir: "Contains the TSV workpackage files used by Workpackage1.wdl, which contain the SAB of each sample in the second column. Without final slash."
        remote_outdir: "Without final slash"
        n_samples_fraction_frequent: "A record is considered frequent iff it was discovered in at least this fraction of the total number of samples."
    }
    
    scatter (chr in chromosomes) {
        call SingleChromosome {
            input:
                chromosome = chr,
                remote_indir = remote_indir,
                remote_workpackages1_dir = remote_workpackages1_dir,
                remote_outdir = remote_outdir,
                n_samples_fraction_frequent = n_samples_fraction_frequent,
                docker_image = docker_image
        }
    }
    call AllChromosomes {
        input:
            chromosomes = chromosomes,
            out_txt = SingleChromosome.out_txt,
            remote_outdir = remote_outdir,
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, HDD:
#
# TOOL                           CPU     RAM   TIME
# download                      250%    100M    15s
# bcftools concat                30%     20M     1m
# bcftools query                100%    150M     4m
# bcftools annotate             100%    150M    10m
# bcftools view frequent        150%    160M     2m
# bcftools view infrequent      300%    140M     5m
#
task SingleChromosome {
    input {
        String chromosome
        String remote_indir
        String remote_workpackages1_dir
        String remote_outdir
        Float n_samples_fraction_frequent
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 4
        Int disk_size_gb = 200
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
        
        #
        function GetNMaleSamples() {
            gcloud storage cp ~{remote_workpackages1_dir}/workpackage1'_*' .
            N_MALE_SAMPLES="0"
            for FILE in $(ls workpackage1_* | sort -V); do
                N=$(cut -f 2 ${FILE} | grep M | wc -l)
                N_MALE_SAMPLES=$(( ${N_MALE_SAMPLES} + ${N} ))
            done
            rm -f workpackage1_*
            echo ${N_MALE_SAMPLES}
        }
        
        
        # ---------------------------- Main program ----------------------------
        
        TEST=$( gcloud storage ls ~{remote_outdir}/~{chromosome}/~{chromosome}.done || echo "0" )
        if [ ${TEST} != "0" ]; then
            # Skipping the chromosome if it has already been processed
            :
        else
            # Localizing all chunks
            gcloud storage ls ~{remote_indir}/~{chromosome}/chunk_'*.bcf' > test.txt
            if grep -q '.bcf' test.txt ; then
                :
            else
                echo "ERROR: ~{chromosome} has no truvari collapse chunks."
                exit
            fi
            ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/~{chromosome}/chunk_'*.bcf*' .
            ls chunk_*.bcf | sort -V > chunk_list.txt
            cat chunk_list.txt
            df -h 1>&2
        
            # Concatenating all chunks
            ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list chunk_list.txt --output-type b --output out.bcf
            df -h 1>&2
            rm -rf chunk_* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
            
            # Enforcing a distinct ID in every record, and annotating every
            # record with the number of samples it occurs in. Note that the
            # latter is not equal to the QUAL field in input to truvari collapse
            # upstream, so we have to recompute this number.
            CHR=~{chromosome}
            CHR=${CHR#chr}
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%COUNT(GT="alt")\n' in.bcf | awk -v id=${CHR} 'BEGIN { FS="\t"; OFS="\t"; i=0; } { $3=sprintf("%s_%d",id,i++); print $0 }' | bgzip -c > annotations.tsv.gz
            tabix -@ ${N_THREADS} -s1 -b2 -e2 annotations.tsv.gz
            echo '##INFO=<ID=N_DISCOVERY_SAMPLES,Number=1,Type=Integer,Description="Number of samples where the record was discovered">' > header.txt
            ${TIME_COMMAND} bcftools annotate --header-lines header.txt --annotations annotations.tsv.gz --columns CHROM,POS,ID,REF,ALT,N_DISCOVERY_SAMPLES --output-type b in.bcf --output out.bcf
            df -h 1>&2
            rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
            gcloud storage cp in.bcf ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf
            gcloud storage cp in.bcf.csi ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf.csi
        
            # Separating frequent and infrequent records
            if [ ~{chromosome} = "chrY" ]; then
                N_SAMPLES=$(GetNMaleSamples)
            else
                N_SAMPLES=$(bcftools view --header-only in.bcf | tail -n 1 | tr '\t' '\n' | tail -n +10 | wc -l)
            fi
            MIN_N_SAMPLES=$(echo "scale=2; ~{n_samples_fraction_frequent} * ${N_SAMPLES}" | bc)
            MIN_N_SAMPLES=$(echo ${MIN_N_SAMPLES} | cut -d . -f 1)
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --drop-genotypes --include 'N_DISCOVERY_SAMPLES>='${MIN_N_SAMPLES} --output-type b in.bcf --output frequent.bcf &
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS}                  --include 'N_DISCOVERY_SAMPLES<'${MIN_N_SAMPLES}  --output-type b in.bcf --output infrequent.bcf &
            wait
            df -h 1>&2
            ${TIME_COMMAND} bcftools index --threads $(( ${N_THREADS} / 2 )) -f frequent.bcf &
            ${TIME_COMMAND} bcftools index --threads $(( ${N_THREADS} / 2 )) -f infrequent.bcf &
            wait
            gcloud storage mv frequent.'bcf*' infrequent.'bcf*' ~{remote_outdir}/~{chromosome}/
            touch ~{chromosome}.done
            gcloud storage mv ~{chromosome}.done ~{remote_outdir}/~{chromosome}/
        fi
        echo "~{chromosome}" > out.txt
    >>>
    
    output {
        File out_txt = "out.txt"
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


# Performance on 12'680 samples, 15x, GRCh38, CAL_SENS<=0.999, HDD:
#
# TOOL                           CPU     RAM        TIME
# concat truvari                  1%     9KB        50m 
# concat frequent                 3%     8KB        16s
# concat infrequent               1%     9KB        40m
#
task AllChromosomes {
    input {
        Array[String] chromosomes
        Array[File] out_txt
        String remote_outdir
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 4
        Int disk_size_gb = 200
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
        
        
        # Localizing
        CHROMOSOMES=~{sep=',' chromosomes}
        echo ${CHROMOSOMES} | tr ',' '\n' > chr_list.txt
        rm -f file_list1.txt file_list2.txt file_list3.txt
        while read CHROMOSOME; do
            TEST=$( gcloud storage ls ~{remote_outdir}/${CHROMOSOME}/truvari_collapsed.bcf || echo 1 )
            if [ ${TEST} -eq 1 ]; then
                echo "ERROR: ${CHROMOSOME} has not been truvari collapsed."
                exit
            fi
            gcloud storage cp ~{remote_outdir}/${CHROMOSOME}/'*.bcf*' .
            mv truvari_collapsed.bcf ${CHROMOSOME}_truvari_collapsed.bcf
            mv truvari_collapsed.bcf.csi ${CHROMOSOME}_truvari_collapsed.bcf.csi
            mv frequent.bcf ${CHROMOSOME}_frequent.bcf
            mv frequent.bcf.csi ${CHROMOSOME}_frequent.bcf.csi
            mv infrequent.bcf ${CHROMOSOME}_infrequent.bcf
            mv infrequent.bcf.csi ${CHROMOSOME}_infrequent.bcf.csi
            echo ${CHROMOSOME}_truvari_collapsed.bcf >> file_list1.txt
            echo ${CHROMOSOME}_frequent.bcf >> file_list2.txt
            echo ${CHROMOSOME}_infrequent.bcf >> file_list3.txt
        done < chr_list.txt
        
        # Ensuring that all files have exactly the same header
        FIRST_CHROMOSOME=$(head -n chr_list.txt)
        bcftools view --header-only ${FIRST_CHROMOSOME}_truvari_collapsed.bcf > header.txt
        while read CHROMOSOME; do
            ${TIME_COMMAND} bcftools reheader --header header.txt ${CHROMOSOME}_truvari_collapsed.bcf --output tmp1.bcf &
            ${TIME_COMMAND} bcftools reheader --header header.txt ${CHROMOSOME}_frequent.bcf --output tmp2.bcf &
            ${TIME_COMMAND} bcftools reheader --header header.txt ${CHROMOSOME}_infrequent.bcf --output tmp3.bcf &
            wait
            rm -f ${CHROMOSOME}_truvari_collapsed.bcf* ; mv tmp1.bcf ${CHROMOSOME}_truvari_collapsed.bcf
            rm -f ${CHROMOSOME}_frequent.bcf* ; mv tmp2.bcf ${CHROMOSOME}_frequent.bcf
            rm -f ${CHROMOSOME}_infrequent.bcf* ; mv tmp3.bcf ${CHROMOSOME}_infrequent.bcf
            bcftools index --threads ${N_THREADS} -f ${CHROMOSOME}_truvari_collapsed.bcf &
            bcftools index --threads ${N_THREADS} -f ${CHROMOSOME}_frequent.bcf &
            bcftools index --threads ${N_THREADS} -f ${CHROMOSOME}_infrequent.bcf &
            wait
        done < chr_list.txt
        
        # Concatenating
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list file_list1.txt --output-type b --output truvari_collapsed.bcf &
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list file_list2.txt --output-type b --output frequent.bcf &
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list file_list3.txt --output-type b --output infrequent.bcf &
        wait
        bcftools index --threads ${N_THREADS} -f truvari_collapsed.bcf &
        bcftools index --threads ${N_THREADS} -f frequent.bcf &
        bcftools index --threads ${N_THREADS} -f infrequent.bcf &
        wait
        
        # Uploading
        gcloud storage mv truvari_collapsed.'bcf*' frequent.'bcf*' infrequent.'bcf*' ~{remote_outdir}/
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
