version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/slee:FilterIntrasampleJointVcfFilteringDev/versions/2/plain-WDL/descriptor" as JointVcfFiltering

struct RuntimeAttributes {
    Int? cpu
    Int? command_mem_gb
    Int? additional_mem_gb
    Int? disk_size_gb
    Int? boot_disk_size_gb
    Boolean? use_ssd
    Int? preemptible
    Int? max_retries
}

workflow FilterIntrasampleDevPhase2 {
    input {
        String sample_name
        
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_vcf_gz_tbi
        File regenotyped_vcf_gz
        File regenotyped_vcf_gz_tbi
        
        File training_resource_vcf_gz
        File training_resource_vcf_gz_tbi
        File training_resource_bed

        Array[String] annotations = ["KS_1","KS_2","SQ","GQ","DP","AD_NON_ALT","AD_ALL","GT_COUNT","SUPP_PAV","SUPP_SNIFFLES","SUPP_PBSV","SVLEN"]
        String? model_backend = "PYTHON_SCRIPT"
        File? training_python_script
        File? scoring_python_script
        File? hyperparameters_json
        String? extract_extra_args = "--mode INDEL"
        String? train_extra_args = "--mode INDEL --verbosity DEBUG"
        String? score_extra_args = "--mode INDEL --mnp-type INDEL --verbosity DEBUG"
        String gatk_docker
        File? gatk_override
        File? monitoring_script
        RuntimeAttributes? extract_runtime_attributes
        RuntimeAttributes? train_runtime_attributes
        RuntimeAttributes? score_runtime_attributes
    }

    call PreprocessVCF {
        input:
            sample_name = sample_name,
            truvari_collapsed_vcf_gz = truvari_collapsed_vcf_gz,
            truvari_collapsed_vcf_gz_tbi = truvari_collapsed_vcf_gz_tbi,
            regenotyped_vcf_gz = regenotyped_vcf_gz,
            regenotyped_vcf_gz_tbi = regenotyped_vcf_gz_tbi,
            training_resource_vcf_gz = training_resource_vcf_gz,
            training_resource_vcf_gz_tbi = training_resource_vcf_gz_tbi,
            output_prefix = sample_name,
            monitoring_script = monitoring_script
    }

    call IdentifyTrainingSites {
        input:
            preprocessed_vcf_gz = PreprocessVCF.preprocessed_vcf_gz,
            preprocessed_vcf_gz_tbi = PreprocessVCF.preprocessed_vcf_gz_tbi,
            preprocessed_training_resource_vcf_gz = PreprocessVCF.preprocessed_training_resource_vcf_gz,
            preprocessed_training_resource_vcf_gz_tbi = PreprocessVCF.preprocessed_training_resource_vcf_gz_tbi,
            training_resource_bed = training_resource_bed,
            monitoring_script = monitoring_script
    }

    call JointVcfFiltering.JointVcfFiltering {
        input:
            input_vcfs = [PreprocessVCF.preprocessed_vcf_gz],
            input_vcf_idxs = [PreprocessVCF.preprocessed_vcf_gz_tbi],
            sites_only_vcf = PreprocessVCF.preprocessed_vcf_gz,
            sites_only_vcf_idx = PreprocessVCF.preprocessed_vcf_gz_tbi,
            output_prefix = sample_name,
            annotations = annotations,
            resource_vcf = IdentifyTrainingSites.training_vcf_gz,
            resource_vcf_idx = IdentifyTrainingSites.training_vcf_gz_tbi,
            model_backend = model_backend,
            training_python_script = training_python_script,
            hyperparameters_json = hyperparameters_json,
            scoring_python_script = scoring_python_script,
            extract_extra_args = extract_extra_args,
            train_extra_args = train_extra_args,
            score_extra_args = score_extra_args,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            extract_runtime_attributes = extract_runtime_attributes,
            train_runtime_attributes = train_runtime_attributes,
            score_runtime_attributes = score_runtime_attributes,
            monitoring_script = monitoring_script
    }

    call PostprocessVCF {
        input:
            scored_vcf_gz = JointVcfFiltering.scored_vcfs[0],
            scored_vcf_gz_tbi = JointVcfFiltering.scored_vcf_idxs[0],
            output_prefix = sample_name,
            monitoring_script = monitoring_script
    }

    output {
        File extracted_annotations_hdf5 = JointVcfFiltering.extracted_annotations_hdf5
        File? extracted_unlabeled_annotations_hdf5 = JointVcfFiltering.extracted_unlabeled_annotations_hdf5
        File extracted_vcf_gz = JointVcfFiltering.extracted_vcf
        File extracted_vcf_gz_tbi = JointVcfFiltering.extracted_vcf_idx

        Array[File] model_files = JointVcfFiltering.model_files

        File? annotations_hdf5 = JointVcfFiltering.annotations_hdf5s[0]
        File? scores_hdf5 = JointVcfFiltering.scores_hdf5s[0]

        File scored_vcf_gz = PostprocessVCF.postprocessed_vcf_gz
        File scored_vcf_gz_tbi = PostprocessVCF.postprocessed_vcf_gz_tbi

        Array[File?] monitoring_logs = flatten(
          [
            [PreprocessVCF.monitoring_log],
            [IdentifyTrainingSites.monitoring_log],
            JointVcfFiltering.monitoring_logs,
            [PostprocessVCF.monitoring_log],
          ])
    }
}


# Adds to the truvari-collapsed VCF a number of INFO fields whose values are
# taken from the following FORMAT fields in the kanpig-regenotyped VCF:
#
# KS_1, KS_2, SQ, GQ, DP, AD_NON_ALT, AD_ALL
#
# and from the following INFO fields in the kanpig-regenotyped VCF:
#
# SUPP_PAV, SUPP_SNIFFLES, SUPP_PBSV
#
task PreprocessVCF {
    input {
        String sample_name
        
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_vcf_gz_tbi
        File regenotyped_vcf_gz
        File regenotyped_vcf_gz_tbi
        
        File training_resource_vcf_gz
        File training_resource_vcf_gz_tbi
        
        String output_prefix

        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    String truvari_collapsed_prefix = basename(truvari_collapsed_vcf_gz, ".vcf.gz")
    String training_resource_prefix = basename(training_resource_vcf_gz, ".vcf.gz")

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        if [ -s ~{monitoring_script} ]; then
          bash ~{monitoring_script} > monitoring.log &
        fi

        # The training resource is already correct
        mv ~{training_resource_vcf_gz} ~{training_resource_prefix}.preprocessed.vcf.gz
        mv ~{training_resource_vcf_gz_tbi} ~{training_resource_prefix}.preprocessed.vcf.gz.tbi

        # Creating supplemental header
        touch format.hdr.txt
        for format_annot in SQ GQ DP
        do
            bcftools view -h ~{regenotyped_vcf_gz} | grep ID=$format_annot, | sed -e 's/FORMAT/INFO/g' >> format.hdr.txt
        done
        echo '##INFO=<ID=AD_NON_ALT,Number=1,Type=Integer,Description="Coverage for non-alternate alleles">' >> format.hdr.txt
        echo '##INFO=<ID=AD_ALL,Number=1,Type=Integer,Description="Coverage for all alleles">' >> format.hdr.txt
        echo '##INFO=<ID=KS_1,Number=1,Type=Integer,Description="Kanpig score 1">' >> format.hdr.txt
        echo '##INFO=<ID=KS_2,Number=1,Type=Integer,Description="Kanpig score 2">' >> format.hdr.txt
        echo '##INFO=<ID=GT_COUNT,Number=1,Type=Integer,Description="GT converted into an int in {0,1,2}.">' >> format.hdr.txt
        echo '##INFO=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by PAV">' >> format.hdr.txt
        echo '##INFO=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by Sniffles2">' >> format.hdr.txt
        echo '##INFO=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> format.hdr.txt

        # Ensuring that every record has a unique ID. We need to join by
        # CHROM,POS,ID in what follows, since using CHROM,POS,REF,ALT makes
        # bcftools annotate segfault.
        bcftools view --header-only ~{truvari_collapsed_vcf_gz} > tmp.vcf
        bcftools view --no-header ~{truvari_collapsed_vcf_gz} | awk 'BEGIN { i=0; } { printf("%s\t%s\t%d-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,++i,$3,$4,$5,$6,$7,$8,$9,$10); }' >> tmp.vcf
        bgzip --compress-level 1 tmp.vcf
        tabix -f tmp.vcf.gz
        bcftools view --no-header tmp.vcf.gz | head -n 10 || echo "0"
        echo 'Annotating ...'
        bcftools query -f '%CHROM\t%POS\t%ID\t[%KS]\t[%SQ]\t[%GQ]\t[%DP]\t[%AD]\t[%GT]\t%INFO/SUPP_PBSV\t%INFO/SUPP_SNIFFLES\t%INFO/SUPP_PAV\n' tmp.vcf.gz | awk '{ \
            KS_1=0; KS_2=0; \
            p=0; \
            for (i=1; i<=length($4); i++) { \
                if (substr($4,i,1)==",") { p=i; break; } \
            } \
            if (p==0) { KS_1=$4; KS_2=$4; } \
            else { KS_1=substr($4,1,p-1); KS_2=substr($4,p+1); } \
            if (KS_1==".") KS_1=0; \
            if (KS_2==".") KS_2=0; \
            \
            SQ=$5; \
            if (SQ==".") SQ=0; \
            \
            GQ=$6; \
            if (GQ==".") GQ=0; \
            \
            DP=$7; \
            if (DP==".") DP=0; \
            \
            AD_NON_ALT=0; AD_ALL=0; \
            p=0; \
            for (i=1; i<=length($8); i++) { \
                if (substr($8,i,1)==",") { p=i; break; } \
            } \
            if (p==0) { AD_NON_ALT=$8; AD_ALL=$8; } \
            else { AD_NON_ALT=substr($8,1,p-1); AD_ALL=substr($8,p+1); } \
            if (AD_NON_ALT==".") AD_NON_ALT=0; \
            if (AD_ALL==".") AD_ALL=0; \
            \
            GT_COUNT=0; \
            if ($9=="0/0" || $9=="0|0" || $9=="./."  || $9==".|." || $9=="./0" || $9==".|0" || $9=="0/." || $9=="0|.") GT_COUNT=0; \
            else if ($9=="0/1" || $9=="0|1" || $9=="1/0" || $9=="1|0" || $9=="./1" || $9==".|1" || $9=="1/." || $9=="1|.") GT_COUNT=1; \
            else if ($9=="1/1" || $9=="1|1") GT_COUNT=2; \
            \
            printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$1,$2,$3,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,$10,$11,$12); \
        }' | bgzip -c > format.tsv.gz
        tabix -s1 -b2 -e2 format.tsv.gz
        bcftools annotate --threads ${N_THREADS} -a format.tsv.gz -h format.hdr.txt -c CHROM,POS,ID,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,SUPP_PBSV,SUPP_SNIFFLES,SUPP_PAV tmp.vcf.gz -Oz -o ~{output_prefix}.vcf.gz
        bcftools view --no-header ~{output_prefix}.vcf.gz | head -n 10 || echo "0"
        bcftools index -t ~{output_prefix}.vcf.gz

        # TODO do we still need to hard-filter SVLEN >= 50 for extract, and should we clear existing filters?
    >>>

    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File preprocessed_vcf_gz = output_prefix+".vcf.gz"
        File preprocessed_vcf_gz_tbi = output_prefix+".vcf.gz.tbi"
        File preprocessed_training_resource_vcf_gz = training_resource_prefix+".preprocessed.vcf.gz"
        File preprocessed_training_resource_vcf_gz_tbi = training_resource_prefix+".preprocessed.vcf.gz.tbi"
        File? monitoring_log = "monitoring.log"
    }
}


task IdentifyTrainingSites {
    input {
        File preprocessed_vcf_gz
        File preprocessed_vcf_gz_tbi
        
        File preprocessed_training_resource_vcf_gz
        File preprocessed_training_resource_vcf_gz_tbi
        File training_resource_bed
        
        String truvari_extra_args = "--sizemin 20 --sizemax 1000000 --sizefilt 20 --pctsize 0.9 --pctseq 0.9 --pick multi"

        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
        set -eou pipefail

        if [ -s ~{monitoring_script} ]; then
          bash ~{monitoring_script} > monitoring.log &
        fi

        truvari bench \
            -b ~{preprocessed_training_resource_vcf_gz} \
            -c ~{preprocessed_vcf_gz} \
            --includebed ~{training_resource_bed} \
            ~{truvari_extra_args} \
            -o truvari/
    >>>

    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File training_vcf_gz = "truvari/tp-comp.vcf.gz"
        File training_vcf_gz_tbi = "truvari/tp-comp.vcf.gz.tbi"
        File? monitoring_log = "monitoring.log"
    }
}


# copy SCORE from INFO field to FORMAT field
task PostprocessVCF {
    input {
        File scored_vcf_gz
        File scored_vcf_gz_tbi
        
        String output_prefix

        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
        set -eou pipefail

        if [ -s ~{monitoring_script} ]; then
          bash ~{monitoring_script} > monitoring.log &
        fi

        # create supplemental header for bcftools annotate
        echo '##FORMAT=<ID=SCORE,Number=1,Type=Float,Description="Score according to the model applied by ScoreVariantAnnotations">' > format.hdr.txt

        # copy SCORE INFO to FORMAT
        bcftools query -f '%CHROM\t%POS\t%ID\t%SCORE\n' ~{scored_vcf_gz} | bgzip -c > format.score.tsv.gz
        tabix -s1 -b2 -e2 format.score.tsv.gz
        bcftools annotate -a format.score.tsv.gz -h format.hdr.txt -c CHROM,POS,ID,FORMAT/SCORE ~{scored_vcf_gz} -Oz -o ~{output_prefix}.score.vcf.gz
        bcftools index -t ~{output_prefix}.score.vcf.gz
    >>>

    runtime {
        docker: "fcunial/callset_integration_phase2"
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File postprocessed_vcf_gz = "~{output_prefix}.score.vcf.gz"
        File postprocessed_vcf_gz_tbi = "~{output_prefix}.score.vcf.gz.tbi"
        File? monitoring_log = "monitoring.log"
    }
}
