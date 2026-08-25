version 1.0


# Essentially identical to `Workpackage8.wdl`, but simplified for ultralong and
# BND VCFs. Concatenates the `truvari collapse` chunks. Ensures that every
# record in output has a globally unique ID (duplicated IDs may arise naturally
# from the previous steps of the pipeline), and an INFO field that counts the
# number of samples it was discovered in.
#
# For ultralong: chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
# For BND:       chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM","chr1_KI270706v1_random","chr1_KI270707v1_random","chr1_KI270708v1_random","chr1_KI270709v1_random","chr1_KI270710v1_random","chr1_KI270711v1_random","chr1_KI270712v1_random","chr1_KI270713v1_random","chr1_KI270714v1_random","chr2_KI270715v1_random","chr2_KI270716v1_random","chr3_GL000221v1_random","chr4_GL000008v2_random","chr5_GL000208v1_random","chr9_KI270717v1_random","chr9_KI270718v1_random","chr9_KI270719v1_random","chr9_KI270720v1_random","chr11_KI270721v1_random","chr14_GL000009v2_random","chr14_GL000225v1_random","chr14_KI270722v1_random","chr14_GL000194v1_random","chr14_KI270723v1_random","chr14_KI270724v1_random","chr14_KI270725v1_random","chr14_KI270726v1_random","chr15_KI270727v1_random","chr16_KI270728v1_random","chr17_GL000205v2_random","chr17_KI270729v1_random","chr17_KI270730v1_random","chr22_KI270731v1_random","chr22_KI270732v1_random","chr22_KI270733v1_random","chr22_KI270734v1_random","chr22_KI270735v1_random","chr22_KI270736v1_random","chr22_KI270737v1_random","chr22_KI270738v1_random","chr22_KI270739v1_random","chrY_KI270740v1_random","chrUn_KI270302v1","chrUn_KI270304v1","chrUn_KI270303v1","chrUn_KI270305v1","chrUn_KI270322v1","chrUn_KI270320v1","chrUn_KI270310v1","chrUn_KI270316v1","chrUn_KI270315v1","chrUn_KI270312v1","chrUn_KI270311v1","chrUn_KI270317v1","chrUn_KI270412v1","chrUn_KI270411v1","chrUn_KI270414v1","chrUn_KI270419v1","chrUn_KI270418v1","chrUn_KI270420v1","chrUn_KI270424v1","chrUn_KI270417v1","chrUn_KI270422v1","chrUn_KI270423v1","chrUn_KI270425v1","chrUn_KI270429v1","chrUn_KI270442v1","chrUn_KI270466v1","chrUn_KI270465v1","chrUn_KI270467v1","chrUn_KI270435v1","chrUn_KI270438v1","chrUn_KI270468v1","chrUn_KI270510v1","chrUn_KI270509v1","chrUn_KI270518v1","chrUn_KI270508v1","chrUn_KI270516v1","chrUn_KI270512v1","chrUn_KI270519v1","chrUn_KI270522v1","chrUn_KI270511v1","chrUn_KI270515v1","chrUn_KI270507v1","chrUn_KI270517v1","chrUn_KI270529v1","chrUn_KI270528v1","chrUn_KI270530v1","chrUn_KI270539v1","chrUn_KI270538v1","chrUn_KI270544v1","chrUn_KI270548v1","chrUn_KI270583v1","chrUn_KI270587v1","chrUn_KI270580v1","chrUn_KI270581v1","chrUn_KI270579v1","chrUn_KI270589v1","chrUn_KI270590v1","chrUn_KI270584v1","chrUn_KI270582v1","chrUn_KI270588v1","chrUn_KI270593v1","chrUn_KI270591v1","chrUn_KI270330v1","chrUn_KI270329v1","chrUn_KI270334v1","chrUn_KI270333v1","chrUn_KI270335v1","chrUn_KI270338v1","chrUn_KI270340v1","chrUn_KI270336v1","chrUn_KI270337v1","chrUn_KI270363v1","chrUn_KI270364v1","chrUn_KI270362v1","chrUn_KI270366v1","chrUn_KI270378v1","chrUn_KI270379v1","chrUn_KI270389v1","chrUn_KI270390v1","chrUn_KI270387v1","chrUn_KI270395v1","chrUn_KI270396v1","chrUn_KI270388v1","chrUn_KI270394v1","chrUn_KI270386v1","chrUn_KI270391v1","chrUn_KI270383v1","chrUn_KI270393v1","chrUn_KI270384v1","chrUn_KI270392v1","chrUn_KI270381v1","chrUn_KI270385v1","chrUn_KI270382v1","chrUn_KI270376v1","chrUn_KI270374v1","chrUn_KI270372v1","chrUn_KI270373v1","chrUn_KI270375v1","chrUn_KI270371v1","chrUn_KI270448v1","chrUn_KI270521v1","chrUn_GL000195v1","chrUn_GL000219v1","chrUn_GL000220v1","chrUn_GL000224v1","chrUn_KI270741v1","chrUn_GL000226v1","chrUn_GL000213v1","chrUn_KI270743v1","chrUn_KI270744v1","chrUn_KI270745v1","chrUn_KI270746v1","chrUn_KI270747v1","chrUn_KI270748v1","chrUn_KI270749v1","chrUn_KI270750v1","chrUn_KI270751v1","chrUn_KI270752v1","chrUn_KI270753v1","chrUn_KI270754v1","chrUn_KI270755v1","chrUn_KI270756v1","chrUn_KI270757v1","chrUn_GL000214v1","chrUn_KI270742v1","chrUn_GL000216v2","chrUn_GL000218v1","chrEBV"]
#
workflow SV_Integration_Workpackage15 {
    input {
        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        String remote_indir
        String remote_outdir
        Int concat_all_naive = 1
        
        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
        chromosomes: "The order of the chromosomes becomes their order in the output VCF."
        remote_indir: "Without final slash"
        remote_outdir: "Without final slash"
        concat_all_naive: "Concatenate chromosomes in a naive (1, default) or non-naive (0) way. Non-naive is necessary when different chromosomes were built by different versions of the pipeline, with slightly different code, and their headers are not exactly identical. Intra-chromosome concatenation is always performed in a naive way."
    }
    
    scatter (chr in chromosomes) {
        call SingleChromosome {
            input:
                chromosome = chr,
                remote_indir = remote_indir,
                remote_outdir = remote_outdir,
                docker_image = docker_image
        }
    }
    call AllChromosomes {
        input:
            chromosomes = chromosomes,
            out_txt = SingleChromosome.out_txt,
            remote_outdir = remote_outdir,
            naive = concat_all_naive,
            docker_image = docker_image
    }
    
    output {
    }
}


# Performance on 12'680 samples, 15x, GRCh38, chr1, HDD:
#
# TOOL                           CPU     RAM   TIME
# bcftools concat                
# bcftools query                
# bcftools annotate             
#
task SingleChromosome {
    input {
        String chromosome
        String remote_indir
        String remote_outdir
        
        String docker_image
        Int n_cpu = 1
        Int ram_size_gb = 2
        Int disk_size_gb = 10
        Int preemptible_number = 0
    }
    parameter_meta {
    }
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        

        # Skipping the chromosome if it has already been concatenated
        TEST=$( gcloud storage ls ~{remote_outdir}/~{chromosome}/~{chromosome}.done || echo "0" )
        if [ "${TEST}" != "0" ]; then
            echo "Warning: ~{chromosome} has already been concatenated."
            echo "~{chromosome}" > out.txt
            exit 0
        fi

        # Localizing all chunks
        TEST=$( gcloud storage ls ~{remote_indir}/~{chromosome}/chunk_'*.bcf' || echo "0" )
        if [ "${TEST}" = "0" ]; then
            echo "Warning: ~{chromosome} has no truvari collapse chunks."
            echo "~{chromosome}" > out.txt
            exit 0
        fi
        ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/~{chromosome}/chunk_'*.bcf*' .
        ls chunk_*.bcf | sort -V > chunk_list.txt
        cat chunk_list.txt
        df -h 1>&2
    
        # Concatenating all chunks
        N_CHUNKS=$(wc -l < chunk_list.txt)
        if [ ${N_CHUNKS} -gt 1 ]; then
            ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list chunk_list.txt --output-type b --output out.bcf
            df -h 1>&2
            rm -rf chunk_* ; mv out.bcf in.bcf
        else
            CHUNK_FILE=$(head -n 1 chunk_list.txt)
            mv ${CHUNK_FILE} in.bcf
            mv ${CHUNK_FILE}.csi in.bcf.csi
        fi
        
        # Enforcing a distinct ID in every record
        CHR=~{chromosome}
        CHR=${CHR#chr}
        ( bcftools view --header-only in.bcf ; bcftools view --no-header in.bcf | awk -v id=${CHR} 'BEGIN { FS="\t"; OFS="\t"; i=0; } { $3=sprintf("%s_%d",id,i++); print $0 }' ) | bgzip -c > out.vcf.gz
        rm in.bcf* ; mv out.vcf.gz in.vcf.gz

        # Annotating every record with the number of samples it occurs in.
        # Note that this is not equal to the QUAL field in input to `truvari
        # collapse` upstream, so it has to be recomputed here.
        ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%COUNT(GT="alt")\n' in.vcf.gz | bgzip -c > annotations.tsv.gz
        tabix -@ ${N_THREADS} -s1 -b2 -e2 annotations.tsv.gz
        echo '##INFO=<ID=N_DISCOVERY_SAMPLES,Number=1,Type=Integer,Description="Number of samples where the record was discovered">' > header.txt
        COLUMNS='CHROM,POS,REF,ALT,~ID,N_DISCOVERY_SAMPLES'
        ${TIME_COMMAND} bcftools annotate --header-lines header.txt --annotations annotations.tsv.gz --columns ${COLUMNS} --output-type b in.vcf.gz --output out.bcf
        df -h 1>&2
        rm -f in.vcf.gz* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf

        # Uploading
        gcloud storage mv in.bcf ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf
        gcloud storage mv in.bcf.csi ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf.csi
        touch ~{chromosome}.done
        gcloud storage mv ~{chromosome}.done ~{remote_outdir}/~{chromosome}/~{chromosome}.done

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
    }
}


# Performance on 12'680 samples, 15x, GRCh38, HDD:
#
# TOOL                           CPU     RAM        TIME
# concat --naive truvari         
#
task AllChromosomes {
    input {
        Array[String] chromosomes
        Array[File] out_txt
        String remote_outdir
        Int naive
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 4
        Int disk_size_gb = 100
        Int preemptible_number = 0
    }
    
    parameter_meta {
    }
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        # Localizing
        CHROMOSOMES=~{sep=',' chromosomes}
        echo ${CHROMOSOMES} | tr ',' '\n' > chr_list.txt
        rm -f file_list.txt
        while read -u 3 CHROMOSOME || [ -n "${CHROMOSOME}" ]; do
            TEST=$( gcloud storage ls ~{remote_outdir}/${CHROMOSOME}/truvari_collapsed.bcf || echo "0" )
            if [ "${TEST}" = "0" ]; then
                echo "Warning: ${CHROMOSOME} has no truvari collapse chunks."
            else
                gcloud storage cp ~{remote_outdir}/${CHROMOSOME}/truvari_collapsed.'bcf*' .
                mv truvari_collapsed.bcf ${CHROMOSOME}_truvari_collapsed.bcf
                mv truvari_collapsed.bcf.csi ${CHROMOSOME}_truvari_collapsed.bcf.csi
                echo ${CHROMOSOME}_truvari_collapsed.bcf >> file_list.txt
            fi
        done 3< chr_list.txt
        
        # Concatenating
        if [ ~{naive} -eq 1 ]; then
            CONCAT_FLAGS="--naive"
        else
            CONCAT_FLAGS=" "
        fi
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} ${CONCAT_FLAGS} --file-list file_list.txt --output-type b --output truvari_collapsed.bcf
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f truvari_collapsed.bcf
        
        # Uploading
        gcloud storage mv truvari_collapsed.'bcf*' ~{remote_outdir}/
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
