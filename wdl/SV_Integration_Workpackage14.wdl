version 1.0


# Runs `truvari collapse` on a small chunk of an ultralong bcftools merged
# cohort VCF. Default truvari parameters have been tested on 10k samples.
#
# Remark: this is almost identical to `SV_Integration_Workpackage7.wdl`: use the
# latter for BNDs.
#
workflow SV_Integration_Workpackage14 {
    input {
        String remote_indir
        String chromosome_id
        File chunks_ids
        String remote_outdir
        
        String truvari_matching_parameters = "--refdist 500 --pctseq 0 --pctsize 0.95 --pctovl 0.0"
        Boolean use_bed
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
        truvari_matching_parameters: "Truvari's definition of a match. The default string contains truvari's defaults. Sequence similarity is turned off for speed reasons."
    }
    
    call Impl {
        input:
            remote_indir = remote_indir,
            chromosome_id = chromosome_id,
            chunks_ids = chunks_ids,
            remote_outdir = remote_outdir,
            
            truvari_matching_parameters = truvari_matching_parameters,
            use_bed = use_bed,
            
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, all calls, one chunk with the 
# whole chr1:
#
# TOOL                      CPU     RAM     TIME
# bcftools view | awk                        20m
# bcftools query            30%      2G       3m
# bcftools annotate        200%      5G       5m
# truvari collapse         100%     27G       2h
# bcftools sort            100%      3G       2m
#
# Peak disk: 10G
#
task Impl {
    input {
        String remote_indir
        String chromosome_id
        File chunks_ids
        String remote_outdir
        
        String truvari_matching_parameters
        Boolean use_bed
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 32
        Int disk_size_gb = 20
        Int preemptible_number = 3
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
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        
        
        
        
        # ----------------------- Steps of the pipeline ------------------------
        
        # Sets QUAL to the number of samples where a record was discovered, to
        # simulate `--keep common` (which is slow on 10k samples) with `--keep
        # maxqual` in `truvari collapse`. See e.g.:
        #
        # https://github.com/ACEnglish/truvari/issues/220#issuecomment-
        # 2830920205
        #
        # Remark: we use the number of samples rather than AC, since we don't
        # trust genotypes at this stage, and since a record being discovered
        # independently in more samples is more informative than it being
        # genotyped more times in fewer samples.
        #
        # Remark: we cannot join annotations just by ID at this stage, since the
        # IDs in the output of bcftools merge are not necessarily all distinct.
        # So it is crucial for this step to be executed on a VCF whose ALTs
        # still contain SVLEN.
        #
        function CopyNSamplesToQual() {
            local CHUNK_ID=$1
        
            mv chunk_${CHUNK_ID}.bcf chunk_${CHUNK_ID}_in.bcf
            mv chunk_${CHUNK_ID}.bcf.csi chunk_${CHUNK_ID}_in.bcf.csi
        
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%COUNT(GT="alt")\n' chunk_${CHUNK_ID}_in.bcf | bgzip -c > chunk_${CHUNK_ID}_annotations.tsv.gz
            tabix -@ ${N_THREADS} -s1 -b2 -e2 chunk_${CHUNK_ID}_annotations.tsv.gz
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations chunk_${CHUNK_ID}_annotations.tsv.gz --columns CHROM,POS,~ID,REF,ALT,QUAL --output-type b chunk_${CHUNK_ID}_in.bcf --output chunk_${CHUNK_ID}_out.bcf
            rm -f chunk_${CHUNK_ID}_in.bcf* ; mv chunk_${CHUNK_ID}_out.bcf chunk_${CHUNK_ID}_in.bcf
            rm -f chunk_${CHUNK_ID}_annotations.tsv.gz
            
            mv chunk_${CHUNK_ID}_in.bcf chunk_${CHUNK_ID}.bcf
        }


        # Removes SVLEN from symbolic ALTs, in order not to interfere with
        # `truvari collapse`.
        #
        function ResetAlts() {
            local CHUNK_ID=$1
            
            date 1>&2
            ( bcftools view --header-only chunk_${CHUNK_ID}.bcf ; bcftools view --no-header chunk_${CHUNK_ID}.bcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                if (substr($5,1,1)=="<") $5 = substr($5,1,4) ">"; \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' ) | bcftools view --output-type z --output chunk_${CHUNK_ID}_out.vcf.gz
            date 1>&2
            rm -f chunk_${CHUNK_ID}.bcf* ; mv chunk_${CHUNK_ID}_out.vcf.gz chunk_${CHUNK_ID}.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t chunk_${CHUNK_ID}.vcf.gz
        }
        
        
        # Remark: in theory we should set `--gt all` to avoid collapsing records
        # that are present in the same sample, since we assume that intra-
        # sample merging has already done that upstream. In practice `--gt all`
        # is too slow on 10k samples.
        #
        # Remark: to further improve speed we could think of dropping genotypes
        # before running truvari collapse. See e.g.:
        #
        # https://github.com/ACEnglish/truvari/issues/220#issuecomment-
        # 2830920205
        #
        # However, this would also discard e.g. SUPP fields that were copied to
        # FORMAT upstream, so it is not correct for our setup.
        #
        function Collapse() {
            local CHUNK_ID=$1
            
            mv chunk_${CHUNK_ID}.vcf.gz chunk_${CHUNK_ID}_in.vcf.gz
            mv chunk_${CHUNK_ID}.vcf.gz.tbi chunk_${CHUNK_ID}_in.vcf.gz.tbi
            
            # Remark: we do not store `removed.vcf` since it's not needed and
            # it can be much bigger than the collapsed output.
            ${TIME_COMMAND} truvari collapse --sizemin 0 --sizemax ${INFINITY} --keep maxqual --gt off ~{truvari_matching_parameters} ${BED_FLAGS} --input chunk_${CHUNK_ID}_in.vcf.gz --output chunk_${CHUNK_ID}_out.vcf --removed-output /dev/null
            df -h 1>&2
            ls -laht 1>&2
            rm -f chunk_${CHUNK_ID}_in.vcf.gz* ; mv chunk_${CHUNK_ID}_out.vcf chunk_${CHUNK_ID}_in.vcf
        
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type b chunk_${CHUNK_ID}_in.vcf --output chunk_${CHUNK_ID}_out.bcf
            df -h 1>&2
            ls -laht 1>&2
            rm -f chunk_${CHUNK_ID}_in.vcf ; mv chunk_${CHUNK_ID}_out.bcf chunk_${CHUNK_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f chunk_${CHUNK_ID}_in.bcf
                
            # Setting to `.` every ID written by truvari collapse, since these
            # can be very long in a large cohort and needlessly inflate output
            # size. Unique IDs genome-wide will be assigned downstream.
            ${TIME_COMMAND} bcftools annotate --remove ID --output-type b chunk_${CHUNK_ID}_in.bcf --output chunk_${CHUNK_ID}_out.bcf
            rm -f chunk_${CHUNK_ID}_in.bcf* ; mv chunk_${CHUNK_ID}_out.bcf chunk_${CHUNK_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f chunk_${CHUNK_ID}_in.bcf
                
            mv chunk_${CHUNK_ID}_in.bcf chunk_${CHUNK_ID}_truvari.bcf
            mv chunk_${CHUNK_ID}_in.bcf.csi chunk_${CHUNK_ID}_truvari.bcf.csi
        }
        
        
        
        
        # ---------------------------- Main program ----------------------------
        
        INFINITY="1000000000"
        truvari --help 1>&2
        df -h 1>&2
        
        if ~{use_bed} ; then
            gcloud storage cp ~{remote_indir}/~{chromosome_id}/included.bed .
            BED_FLAGS="--bed included.bed"
        else 
            BED_FLAGS=" "
        fi
        while read -u 3 CHUNK_ID || [ -n "${CHUNK_ID}" ]; do
            # Skipping the chunk if it has already been processed
            TEST=$( gcloud storage ls ~{remote_outdir}/~{chromosome_id}/chunk_${CHUNK_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi
            
            # Collapsing
            gcloud storage cp ~{remote_indir}/~{chromosome_id}/chunk_${CHUNK_ID}.'bcf*' .
            CopyNSamplesToQual ${CHUNK_ID}
            ResetAlts ${CHUNK_ID}
            Collapse ${CHUNK_ID}
            
            # Uploading
            gcloud storage mv chunk_${CHUNK_ID}_truvari.bcf'*' ~{remote_outdir}/~{chromosome_id}/
            touch chunk_${CHUNK_ID}.done
            gcloud storage mv chunk_${CHUNK_ID}.done ~{remote_outdir}/~{chromosome_id}/
            rm -rf chunk_${CHUNK_ID}*
            ls -laht 1>&2
        done 3< ~{chunks_ids}
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_number
    }
}
