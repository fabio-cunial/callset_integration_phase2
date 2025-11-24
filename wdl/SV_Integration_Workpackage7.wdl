version 1.0


# Runs truvari collapse on a small chunk of a bcftools merged cohort VCF.
# Default arguments are optimal for 10k samples.
#
workflow SV_Integration_Workpackage7 {
    input {
        String chromosome_id
        File chunks_ids
        Boolean use_bed
        
        String remote_indir
        String remote_outdir
    }
    parameter_meta {
        chunks_ids: "One chunk ID per line"
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
    }
    
    call Impl {
        input:
            chromosome_id = chromosome_id,
            chunks_ids = chunks_ids,
            use_bed = use_bed,
            remote_indir = remote_indir,
            remote_outdir = remote_outdir
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, one 5 MB chunk of chr6,
# CAL_SENS<=0.999:
#
# TOOL               CPU     RAM     TIME   OUTPUT VCF
# truvari collapse   100%    2G      3m     500M
# bcftools sort      100%    300M    10s     
#
task Impl {
    input {
        String chromosome_id
        File chunks_ids
        Boolean use_bed
        
        String remote_indir
        String remote_outdir
        
        Int n_cpu = 2
        Int ram_size_gb = 8
        Int disk_size_gb = 50
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        function LocalizeBed() {
            while : ; do
                TEST=$(gsutil -m cp ~{remote_indir}/~{chromosome_id}_included.bed . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading ~{chromosome_id} included BED. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        }
        
        
        function LocalizeChunk() {
            local CHUNK_ID=$1
            
            while : ; do
                TEST=$(gsutil -m cp ~{remote_indir}/~{chromosome_id}_chunk_${CHUNK_ID}.vcf.'gz*' . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading ~{chromosome_id} chunk ${CHUNK_ID}. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        }
        
        
        # Sets QUAL to the number of samples where a record was discovered, to
        # simulate `--keep common` (which is slow on 10k samples) with `--keep
        # maxqual`. See e.g.:
        #
        # https://github.com/ACEnglish/truvari/issues/220#issuecomment-
        # 2830920205
        #
        # Remark: we use the number of samples rather than AC, since we don't
        # trust genotypes at this stage, and since a record being discovered
        # independently in more samples is more informative than it being
        # genotyped more times in fewer samples.
        #
        function CopyNSamplesToQual() {
            local CHUNK_ID=$1
        
            mv ~{chromosome_id}_chunk_${CHUNK_ID}.vcf.gz chunk_${CHUNK_ID}_in.vcf.gz
            mv ~{chromosome_id}_chunk_${CHUNK_ID}.vcf.gz.tbi chunk_${CHUNK_ID}_in.vcf.gz.tbi
        
            # Ensuring that all records are consistently annotated
            ${TIME_COMMAND} bcftools +fill-tags chunk_${CHUNK_ID}_in.vcf.gz -Oz -o chunk_${CHUNK_ID}_out.vcf.gz -- --tags AC_Hom,AC_Het
            rm -f chunk_${CHUNK_ID}_in.vcf.gz* ; mv chunk_${CHUNK_ID}_out.vcf.gz chunk_${CHUNK_ID}_in.vcf.gz ; tabix -f chunk_${CHUNK_ID}_in.vcf.gz
        
            # Copying the number of samples to QUAL.
            # Remark: joining annotations on REF and ALT in addition to ID is
            # important, since IDs are not necessarily all distinct at this
            # stage.
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC_Hom\t%AC_Het\n' chunk_${CHUNK_ID}_in.vcf.gz | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                printf("%s\t%s\t%s\t%s\t%s\t%d\n",$1,$2,$3,$4,$5,$6 / 2 + $7); \
            }' | bgzip -c > chunk_${CHUNK_ID}_annotations.tsv.gz
            tabix -s1 -b2 -e2 chunk_${CHUNK_ID}_annotations.tsv.gz
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations chunk_${CHUNK_ID}_annotations.tsv.gz --columns CHROM,POS,~ID,REF,ALT,QUAL --output-type z chunk_${CHUNK_ID}_in.vcf.gz > chunk_${CHUNK_ID}_out.vcf.gz
            rm -f chunk_${CHUNK_ID}_in.vcf.gz* ; mv chunk_${CHUNK_ID}_out.vcf.gz chunk_${CHUNK_ID}_in.vcf.gz ; tabix -f chunk_${CHUNK_ID}_in.vcf.gz
            rm -f chunk_${CHUNK_ID}_annotations.tsv.gz
            
            mv chunk_${CHUNK_ID}_in.vcf.gz chunk_${CHUNK_ID}_annotated.vcf.gz
            mv chunk_${CHUNK_ID}_in.vcf.gz.tbi chunk_${CHUNK_ID}_annotated.vcf.gz.tbi
        }
        
        
        # Remark: we would like to set `--gt all` to avoid collapsing records
        # that are present in the same sample, since we assume that
        # intra-sample merging has already been upstream. However, `--gt all`
        # is slow on 10k samples.
        #
        # Remark: to further improve speed we could think of dropping genotypes
        # before running truvari collapse, since AC is the key information we
        # need during re-genotyping downstream. See e.g.:
        #
        # https://github.com/ACEnglish/truvari/issues/220#issuecomment-
        # 2830920205
        #
        # However, this would also discard e.g. SUPP fields that were copied to
        # FORMAT upstream, so it is not correct for our setup. It would also
        # make it impossible e.g. to compare precision/recall after collapse to
        # precision/recall after cohort re-genotyping.
        #
        function Collapse() {
            local CHUNK_ID=$1
            
            mv chunk_${CHUNK_ID}_annotated.vcf.gz chunk_${CHUNK_ID}_in.vcf.gz
            mv chunk_${CHUNK_ID}_annotated.vcf.gz.tbi chunk_${CHUNK_ID}_in.vcf.gz.tbi
            
            ${TIME_COMMAND} truvari collapse --sizemin 0 --sizemax ${INFINITY} --keep maxqual --gt off ${BED_FLAGS} --input chunk_${CHUNK_ID}_in.vcf.gz --output chunk_${CHUNK_ID}_out.vcf
            ls -laht
            rm -f chunk_${CHUNK_ID}_in.vcf.gz* ; mv chunk_${CHUNK_ID}_out.vcf chunk_${CHUNK_ID}_in.vcf
        
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z chunk_${CHUNK_ID}_in.vcf > chunk_${CHUNK_ID}_out.vcf.gz
            rm -f chunk_${CHUNK_ID}_in.vcf ; mv chunk_${CHUNK_ID}_out.vcf.gz chunk_${CHUNK_ID}_in.vcf.gz ; tabix -f chunk_${CHUNK_ID}_in.vcf.gz
        
            mv chunk_${CHUNK_ID}_in.vcf.gz ~{chromosome_id}_chunk_${CHUNK_ID}_truvari.vcf.gz
            mv chunk_${CHUNK_ID}_in.vcf.gz.tbi ~{chromosome_id}_chunk_${CHUNK_ID}_truvari.vcf.gz.tbi
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        INFINITY="1000000000"
        if ~{use_bed} ; then
            LocalizeBed
            BED_FLAGS="--bed ~{chromosome_id}_included.bed"
        else 
            BED_FLAGS=" "
        fi
        while read CHUNK_ID; do
            LocalizeChunk ${CHUNK_ID}
            CopyNSamplesToQual ${CHUNK_ID}
            Collapse ${CHUNK_ID}
            
            # Uploading
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{chromosome_id}_chunk_${CHUNK_ID}_truvari.vcf.'gz*' ~{remote_outdir}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading a truvari collapse chunk. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            rm -rf *chunk_${CHUNK_ID}_*
            ls -laht
        done < ~{chunks_ids}
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
