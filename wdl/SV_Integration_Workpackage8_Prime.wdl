version 1.0

# Transforms `cohort_truvari_vcf_gz` into the following BCFs:
# - Frequent: the subset of all records that occur in at least the specified 
#   fraction of all samples.
# - Infrequent: a single BCF per sample, containing all and only the records
#   that occur in the sample and that occur in less than the specified 
#   fraction of all samples.
# Every record in every BCF above has a globally unique integer ID (this is
# necessary for kanpig downstream), and an INFO field that counts the number of
# samples it occurs in, but it does not have FORMAT and SAMPLE columns.
#
# Remark: this is a temporary step to go from `SV_Integration_Workpackage8.wdl`
# to personalized-VCF re-genotyping. It should eventually be merged with
# `SV_Integration_Workpackage8.wdl`.
#
workflow SV_Integration_Workpackage8_Prime {
    input {
        File cohort_truvari_vcf_gz
        File cohort_truvari_tbi
        
        Int min_sv_length
        Int max_sv_length
        Float n_samples_fraction_frequent
        
        String remote_outdir
    }
    parameter_meta {
        n_samples_fraction_frequent: "A record is considered frequent iff it occurs in at least this fraction of the total number of samples."
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            cohort_truvari_vcf_gz = cohort_truvari_vcf_gz,
            cohort_truvari_tbi = cohort_truvari_tbi,
            min_sv_length = min_sv_length,
            max_sv_length = max_sv_length,
            n_samples_fraction_frequent = n_samples_fraction_frequent,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


# Remark: this should probably be parallelized by chromosome in the final
# version.
#
# Performance on 12'680 samples, 15x, GRCh38, chr6, CAL_SENS<=0.999, SSD:
#
# TOOL                              CPU     RAM     TIME
# bcftools view b                   200%    600M    15m
# bcftools query                    100%    600M    3m
# bcftools annotate                 100%    1G      15m
# bcftools view --drop-genotypes    100%    500M    3m
# bcftools +split                   100%    500M    3m
# bcftools view --drop-genotypes    100%    15M     10s
#
task Impl {
    input {
        File cohort_truvari_vcf_gz
        File cohort_truvari_tbi
        
        Int min_sv_length
        Int max_sv_length
        Float n_samples_fraction_frequent
        
        String remote_outdir
        
        Int n_cpu = 4
        Int ram_size_gb = 4
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 10*ceil(size(cohort_truvari_vcf_gz,"GB"))
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"

        
        # Formats a chunk of infrequent VCFs
        function FormatFiles() {
            local CHUNK_FILE=$1
            
            while read FILE; do
                SAMPLE_ID=$(basename ${FILE} .bcf)
                ${TIME_COMMAND} bcftools view --drop-genotypes --include 'COUNT(GT="alt")>0' --output-type b ${FILE} > ${SAMPLE_ID}_infrequent.bcf
                bcftools index ${SAMPLE_ID}_infrequent.bcf
                rm -f ${FILE}
            done < ${CHUNK_FILE}
        }
        
        
        # -------------------------- Main program ------------------------------
        
        mv ~{cohort_truvari_vcf_gz} in.vcf.gz
        mv ~{cohort_truvari_tbi} in.vcf.gz.tbi
        
        # Converting to BCF, to speed up all steps downstream.
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'ABS(SVLEN)>='~{min_sv_length}' && ABS(SVLEN)<='~{max_sv_length} --output-type b in.vcf.gz > out.bcf
        rm -f in.vcf.gz* ; mv out.bcf in.bcf ; bcftools index in.bcf
        
        # Enforcing a distinct integer ID in every record, and annotating every
        # record with the number of samples it occurs in.
        ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%ID\t%COUNT(GT="alt")\n' in.bcf | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { $3=++i; gsub(/;/,"_",$6); print $0 }' | bgzip -c > annotations.tsv.gz
        tabix -s1 -b2 -e2 annotations.tsv.gz
        echo '##INFO=<ID=N_DISCOVERY_SAMPLES,Number=1,Type=Integer,Description="Number of samples where the record was discovered">' > header.txt
        echo '##INFO=<ID=ORIGINAL_ID,Number=1,Type=String,Description="Original ID from truvari collapse">' >> header.txt
        ${TIME_COMMAND} bcftools annotate --header-lines header.txt --annotations annotations.tsv.gz --columns CHROM,POS,ID,REF,ALT,ORIGINAL_ID,N_DISCOVERY_SAMPLES --output-type z in.bcf > out.bcf
        rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index in.bcf
        
        # Separating frequent and infrequent records, and projecting to a
        # distinct VCF every column of the infrequent VCF.
        bcftools view --header-only in.bcf | tail -n 1 | tr '\t' '\n' | tail -n +10 > samples.txt
        N_SAMPLES=$(wc -l < samples.txt)
        MIN_N_SAMPLES=$(echo "scale=2; ~{n_samples_fraction_frequent} * ${N_SAMPLES}" | bc)
        MIN_N_SAMPLES=$(echo ${MIN_N_SAMPLES} | cut -d . -f 1)
        mkdir ./infrequent/
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --drop-genotypes --include 'N_DISCOVERY_SAMPLES>='${MIN_N_SAMPLES} --output-type b in.bcf > frequent.bcf &
        ${TIME_COMMAND} bcftools +split --include 'N_DISCOVERY_SAMPLES<'${MIN_N_SAMPLES} --samples-file samples.txt --output-type b --output ./infrequent/ in.bcf &
        wait
        ${TIME_COMMAND} bcftools index frequent.bcf
        rm -f in.bcf
        
        # Formatting the infrequent BCFs
        ls ./infrequent/*.bcf > infrequent_files.txt
        N_FILES_PER_THREAD=$(( ${N_SAMPLES} / ${N_THREADS} ))
        split -l ${N_FILES_PER_THREAD} -d -a 4 infrequent_files.txt chunk_
        for CHUNK_FILE in $(ls chunk_*); do
            FormatFiles ${CHUNK_FILE} &
        done
        wait
        
        # Uploading
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ./frequent.'bcf*' './infrequent/*_infrequent.bcf*' ~{remote_outdir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading frequent and infrequent VCFs. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
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
