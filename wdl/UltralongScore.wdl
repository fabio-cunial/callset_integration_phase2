version 1.0


# For DEL:
#
# annotations_custom = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","BIN_BEFORE_COVERAGE","BIN_LEFT_COVERAGE","BIN_1_COVERAGE","BIN_2_COVERAGE","BIN_3_COVERAGE","BIN_4_COVERAGE","BIN_5_COVERAGE","BIN_6_COVERAGE","BIN_7_COVERAGE","BIN_8_COVERAGE","BIN_9_COVERAGE","BIN_10_COVERAGE","BIN_RIGHT_COVERAGE","BIN_AFTER_COVERAGE","BIN_LEFT_MAPQ","BIN_RIGHT_MAPQ","BIN_LEFT_SECONDARY","BIN_RIGHT_SECONDARY","LL","LR","RL","RR","LL_RL_1","LL_RL_2","LL_RL_3","LL_RL_4","LL_RR_1","LL_RR_2","LL_RR_3","LL_RR_4","LR_RL_1","LR_RL_2","LR_RL_3","LR_RL_4","LR_RR_1","LR_RR_2","LR_RR_3","LR_RR_4"]
# annotations_fex = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK"]
# annotations_sniffles = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","SNIFFLES_GT_COUNT","SNIFFLES_GQ","SNIFFLES_DR","SNIFFLES_DV"]
# annotations_cutefc = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","CUTEFC_GT_COUNT","CUTEFC_GQ","CUTEFC_DR","CUTEFC_DV","CUTEFC_PL_1","CUTEFC_PL_2","CUTEFC_PL_3","CUTEFC_CIPOS_1","CUTEFC_CIPOS_2","CUTEFC_CILEN_1","CUTEFC_CILEN_2","CUTEFC_RE","CUTEFC_STRAND"]
# annotations_lrcaller_left = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","LRCALLER_GTCOUNT1_left","LRCALLER_AD11_left","LRCALLER_AD12_left","LRCALLER_AD13_left","LRCALLER_VA11_left","LRCALLER_VA12_left","LRCALLER_VA13_left","LRCALLER_PL11_left","LRCALLER_PL12_left","LRCALLER_PL13_left","LRCALLER_GTCOUNT2_left","LRCALLER_AD21_left","LRCALLER_AD22_left","LRCALLER_AD23_left","LRCALLER_VA21_left","LRCALLER_VA22_left","LRCALLER_VA23_left","LRCALLER_PL21_left","LRCALLER_PL22_left","LRCALLER_PL23_left","LRCALLER_GTCOUNT3_left","LRCALLER_AD31_left","LRCALLER_AD32_left","LRCALLER_AD33_left","LRCALLER_VA31_left","LRCALLER_VA32_left","LRCALLER_VA33_left","LRCALLER_PL31_left","LRCALLER_PL32_left","LRCALLER_PL33_left","LRCALLER_GTCOUNT4_left","LRCALLER_AD41_left","LRCALLER_AD42_left","LRCALLER_AD43_left","LRCALLER_VA41_left","LRCALLER_VA42_left","LRCALLER_VA43_left","LRCALLER_PL41_left","LRCALLER_PL42_left","LRCALLER_PL43_left","LRCALLER_GTCOUNT5_left","LRCALLER_AD51_left","LRCALLER_AD52_left","LRCALLER_AD53_left","LRCALLER_VA51_left","LRCALLER_VA52_left","LRCALLER_VA53_left","LRCALLER_PL51_left","LRCALLER_PL52_left","LRCALLER_PL53_left"]
# annotations_lrcaller_right = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","LRCALLER_GTCOUNT1_right","LRCALLER_AD11_right","LRCALLER_AD12_right","LRCALLER_AD13_right","LRCALLER_VA11_right","LRCALLER_VA12_right","LRCALLER_VA13_right","LRCALLER_PL11_right","LRCALLER_PL12_right","LRCALLER_PL13_right","LRCALLER_GTCOUNT2_right","LRCALLER_AD21_right","LRCALLER_AD22_right","LRCALLER_AD23_right","LRCALLER_VA21_right","LRCALLER_VA22_right","LRCALLER_VA23_right","LRCALLER_PL21_right","LRCALLER_PL22_right","LRCALLER_PL23_right","LRCALLER_GTCOUNT3_right","LRCALLER_AD31_right","LRCALLER_AD32_right","LRCALLER_AD33_right","LRCALLER_VA31_right","LRCALLER_VA32_right","LRCALLER_VA33_right","LRCALLER_PL31_right","LRCALLER_PL32_right","LRCALLER_PL33_right","LRCALLER_GTCOUNT4_right","LRCALLER_AD41_right","LRCALLER_AD42_right","LRCALLER_AD43_right","LRCALLER_VA41_right","LRCALLER_VA42_right","LRCALLER_VA43_right","LRCALLER_PL41_right","LRCALLER_PL42_right","LRCALLER_PL43_right","LRCALLER_GTCOUNT5_right","LRCALLER_AD51_right","LRCALLER_AD52_right","LRCALLER_AD53_right","LRCALLER_VA51_right","LRCALLER_VA52_right","LRCALLER_VA53_right","LRCALLER_PL51_right","LRCALLER_PL52_right","LRCALLER_PL53_right"]
# annotations_lrcaller_all = [ "SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","LRCALLER_GTCOUNT1_left","LRCALLER_AD11_left","LRCALLER_AD12_left","LRCALLER_AD13_left","LRCALLER_VA11_left","LRCALLER_VA12_left","LRCALLER_VA13_left","LRCALLER_PL11_left","LRCALLER_PL12_left","LRCALLER_PL13_left","LRCALLER_GTCOUNT2_left","LRCALLER_AD21_left","LRCALLER_AD22_left","LRCALLER_AD23_left","LRCALLER_VA21_left","LRCALLER_VA22_left","LRCALLER_VA23_left","LRCALLER_PL21_left","LRCALLER_PL22_left","LRCALLER_PL23_left","LRCALLER_GTCOUNT3_left","LRCALLER_AD31_left","LRCALLER_AD32_left","LRCALLER_AD33_left","LRCALLER_VA31_left","LRCALLER_VA32_left","LRCALLER_VA33_left","LRCALLER_PL31_left","LRCALLER_PL32_left","LRCALLER_PL33_left","LRCALLER_GTCOUNT4_left","LRCALLER_AD41_left","LRCALLER_AD42_left","LRCALLER_AD43_left","LRCALLER_VA41_left","LRCALLER_VA42_left","LRCALLER_VA43_left","LRCALLER_PL41_left","LRCALLER_PL42_left","LRCALLER_PL43_left","LRCALLER_GTCOUNT5_left","LRCALLER_AD51_left","LRCALLER_AD52_left","LRCALLER_AD53_left","LRCALLER_VA51_left","LRCALLER_VA52_left","LRCALLER_VA53_left","LRCALLER_PL51_left","LRCALLER_PL52_left","LRCALLER_PL53_left",
#   "LRCALLER_GTCOUNT1_right","LRCALLER_AD11_right","LRCALLER_AD12_right","LRCALLER_AD13_right","LRCALLER_VA11_right","LRCALLER_VA12_right","LRCALLER_VA13_right","LRCALLER_PL11_right","LRCALLER_PL12_right","LRCALLER_PL13_right","LRCALLER_GTCOUNT2_right","LRCALLER_AD21_right","LRCALLER_AD22_right","LRCALLER_AD23_right","LRCALLER_VA21_right","LRCALLER_VA22_right","LRCALLER_VA23_right","LRCALLER_PL21_right","LRCALLER_PL22_right","LRCALLER_PL23_right","LRCALLER_GTCOUNT3_right","LRCALLER_AD31_right","LRCALLER_AD32_right","LRCALLER_AD33_right","LRCALLER_VA31_right","LRCALLER_VA32_right","LRCALLER_VA33_right","LRCALLER_PL31_right","LRCALLER_PL32_right","LRCALLER_PL33_right","LRCALLER_GTCOUNT4_right","LRCALLER_AD41_right","LRCALLER_AD42_right","LRCALLER_AD43_right","LRCALLER_VA41_right","LRCALLER_VA42_right","LRCALLER_VA43_right","LRCALLER_PL41_right","LRCALLER_PL42_right","LRCALLER_PL43_right","LRCALLER_GTCOUNT5_right","LRCALLER_AD51_right","LRCALLER_AD52_right","LRCALLER_AD53_right","LRCALLER_VA51_right","LRCALLER_VA52_right","LRCALLER_VA53_right","LRCALLER_PL51_right","LRCALLER_PL52_right","LRCALLER_PL53_right"
# ]
# annotations_all = [ "SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","BIN_BEFORE_COVERAGE","BIN_LEFT_COVERAGE","BIN_1_COVERAGE","BIN_2_COVERAGE","BIN_3_COVERAGE","BIN_4_COVERAGE","BIN_5_COVERAGE","BIN_6_COVERAGE","BIN_7_COVERAGE","BIN_8_COVERAGE","BIN_9_COVERAGE","BIN_10_COVERAGE","BIN_RIGHT_COVERAGE","BIN_AFTER_COVERAGE","BIN_LEFT_MAPQ","BIN_RIGHT_MAPQ","BIN_LEFT_SECONDARY","BIN_RIGHT_SECONDARY","LL","LR","RL","RR","LL_RL_1","LL_RL_2","LL_RL_3","LL_RL_4","LL_RR_1","LL_RR_2","LL_RR_3","LL_RR_4","LR_RL_1","LR_RL_2","LR_RL_3","LR_RL_4","LR_RR_1","LR_RR_2","LR_RR_3","LR_RR_4",
#  "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK",
#  "SNIFFLES_GT_COUNT","SNIFFLES_GQ","SNIFFLES_DR","SNIFFLES_DV",
#  "CUTEFC_GT_COUNT","CUTEFC_GQ","CUTEFC_DR","CUTEFC_DV","CUTEFC_PL_1","CUTEFC_PL_2","CUTEFC_PL_3","CUTEFC_CIPOS_1","CUTEFC_CIPOS_2","CUTEFC_CILEN_1","CUTEFC_CILEN_2","CUTEFC_RE","CUTEFC_STRAND",
#  "LRCALLER_GTCOUNT1_left","LRCALLER_AD11_left","LRCALLER_AD12_left","LRCALLER_AD13_left","LRCALLER_VA11_left","LRCALLER_VA12_left","LRCALLER_VA13_left","LRCALLER_PL11_left","LRCALLER_PL12_left","LRCALLER_PL13_left","LRCALLER_GTCOUNT2_left","LRCALLER_AD21_left","LRCALLER_AD22_left","LRCALLER_AD23_left","LRCALLER_VA21_left","LRCALLER_VA22_left","LRCALLER_VA23_left","LRCALLER_PL21_left","LRCALLER_PL22_left","LRCALLER_PL23_left","LRCALLER_GTCOUNT3_left","LRCALLER_AD31_left","LRCALLER_AD32_left","LRCALLER_AD33_left","LRCALLER_VA31_left","LRCALLER_VA32_left","LRCALLER_VA33_left","LRCALLER_PL31_left","LRCALLER_PL32_left","LRCALLER_PL33_left","LRCALLER_GTCOUNT4_left","LRCALLER_AD41_left","LRCALLER_AD42_left","LRCALLER_AD43_left","LRCALLER_VA41_left","LRCALLER_VA42_left","LRCALLER_VA43_left","LRCALLER_PL41_left","LRCALLER_PL42_left","LRCALLER_PL43_left","LRCALLER_GTCOUNT5_left","LRCALLER_AD51_left","LRCALLER_AD52_left","LRCALLER_AD53_left","LRCALLER_VA51_left","LRCALLER_VA52_left","LRCALLER_VA53_left","LRCALLER_PL51_left","LRCALLER_PL52_left","LRCALLER_PL53_left",
#  "LRCALLER_GTCOUNT1_right","LRCALLER_AD11_right","LRCALLER_AD12_right","LRCALLER_AD13_right","LRCALLER_VA11_right","LRCALLER_VA12_right","LRCALLER_VA13_right","LRCALLER_PL11_right","LRCALLER_PL12_right","LRCALLER_PL13_right","LRCALLER_GTCOUNT2_right","LRCALLER_AD21_right","LRCALLER_AD22_right","LRCALLER_AD23_right","LRCALLER_VA21_right","LRCALLER_VA22_right","LRCALLER_VA23_right","LRCALLER_PL21_right","LRCALLER_PL22_right","LRCALLER_PL23_right","LRCALLER_GTCOUNT3_right","LRCALLER_AD31_right","LRCALLER_AD32_right","LRCALLER_AD33_right","LRCALLER_VA31_right","LRCALLER_VA32_right","LRCALLER_VA33_right","LRCALLER_PL31_right","LRCALLER_PL32_right","LRCALLER_PL33_right","LRCALLER_GTCOUNT4_right","LRCALLER_AD41_right","LRCALLER_AD42_right","LRCALLER_AD43_right","LRCALLER_VA41_right","LRCALLER_VA42_right","LRCALLER_VA43_right","LRCALLER_PL41_right","LRCALLER_PL42_right","LRCALLER_PL43_right","LRCALLER_GTCOUNT5_right","LRCALLER_AD51_right","LRCALLER_AD52_right","LRCALLER_AD53_right","LRCALLER_VA51_right","LRCALLER_VA52_right","LRCALLER_VA53_right","LRCALLER_PL51_right","LRCALLER_PL52_right","LRCALLER_PL53_right"
#]
#
workflow UltralongScore {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File resource_vcf_gz
        File resource_vcf_gz_tbi
        String remote_outdir
                
        File? training_resource_bed
        File reference_fa
        File reference_fai

        Array[String] annotations_custom = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","BIN_POS","BIN_POINT_MAPQ","BIN_POINT_SECONDARY","PL","PR","PL_PL_1","PL_PL_2","PL_PL_3","PL_PL_4","PL_PR_1","PL_PR_2","PL_PR_3","PL_PR_4","PR_PR_1","PR_PR_2","PR_PR_3","PR_PR_4"]
        Array[String] annotations_fex = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK"]
        Array[String] annotations_sniffles = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","SNIFFLES_GT_COUNT","SNIFFLES_GQ","SNIFFLES_DR","SNIFFLES_DV"]
        Array[String] annotations_cutefc = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","CUTEFC_GT_COUNT","CUTEFC_GQ","CUTEFC_DR","CUTEFC_DV","CUTEFC_PL_1","CUTEFC_PL_2","CUTEFC_PL_3","CUTEFC_CIPOS_1","CUTEFC_CIPOS_2","CUTEFC_CILEN_1","CUTEFC_CILEN_2","CUTEFC_RE","CUTEFC_STRAND"]
        Array[String] annotations_lrcaller_left = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","LRCALLER_GTCOUNT1_left","LRCALLER_AD11_left","LRCALLER_AD12_left","LRCALLER_AD13_left","LRCALLER_VA11_left","LRCALLER_VA12_left","LRCALLER_VA13_left","LRCALLER_PL11_left","LRCALLER_PL12_left","LRCALLER_PL13_left","LRCALLER_GTCOUNT2_left","LRCALLER_AD21_left","LRCALLER_AD22_left","LRCALLER_AD23_left","LRCALLER_VA21_left","LRCALLER_VA22_left","LRCALLER_VA23_left","LRCALLER_PL21_left","LRCALLER_PL22_left","LRCALLER_PL23_left","LRCALLER_GTCOUNT3_left","LRCALLER_AD31_left","LRCALLER_AD32_left","LRCALLER_AD33_left","LRCALLER_VA31_left","LRCALLER_VA32_left","LRCALLER_VA33_left","LRCALLER_PL31_left","LRCALLER_PL32_left","LRCALLER_PL33_left","LRCALLER_GTCOUNT4_left","LRCALLER_AD41_left","LRCALLER_AD42_left","LRCALLER_AD43_left","LRCALLER_VA41_left","LRCALLER_VA42_left","LRCALLER_VA43_left","LRCALLER_PL41_left","LRCALLER_PL42_left","LRCALLER_PL43_left","LRCALLER_GTCOUNT5_left","LRCALLER_AD51_left","LRCALLER_AD52_left","LRCALLER_AD53_left","LRCALLER_VA51_left","LRCALLER_VA52_left","LRCALLER_VA53_left","LRCALLER_PL51_left","LRCALLER_PL52_left","LRCALLER_PL53_left"]
        Array[String] annotations_lrcaller_right = ["SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV"]
        Array[String] annotations_lrcaller_all = [ "SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","LRCALLER_GTCOUNT1_left","LRCALLER_AD11_left","LRCALLER_AD12_left","LRCALLER_AD13_left","LRCALLER_VA11_left","LRCALLER_VA12_left","LRCALLER_VA13_left","LRCALLER_PL11_left","LRCALLER_PL12_left","LRCALLER_PL13_left","LRCALLER_GTCOUNT2_left","LRCALLER_AD21_left","LRCALLER_AD22_left","LRCALLER_AD23_left","LRCALLER_VA21_left","LRCALLER_VA22_left","LRCALLER_VA23_left","LRCALLER_PL21_left","LRCALLER_PL22_left","LRCALLER_PL23_left","LRCALLER_GTCOUNT3_left","LRCALLER_AD31_left","LRCALLER_AD32_left","LRCALLER_AD33_left","LRCALLER_VA31_left","LRCALLER_VA32_left","LRCALLER_VA33_left","LRCALLER_PL31_left","LRCALLER_PL32_left","LRCALLER_PL33_left","LRCALLER_GTCOUNT4_left","LRCALLER_AD41_left","LRCALLER_AD42_left","LRCALLER_AD43_left","LRCALLER_VA41_left","LRCALLER_VA42_left","LRCALLER_VA43_left","LRCALLER_PL41_left","LRCALLER_PL42_left","LRCALLER_PL43_left","LRCALLER_GTCOUNT5_left","LRCALLER_AD51_left","LRCALLER_AD52_left","LRCALLER_AD53_left","LRCALLER_VA51_left","LRCALLER_VA52_left","LRCALLER_VA53_left","LRCALLER_PL51_left","LRCALLER_PL52_left","LRCALLER_PL53_left" ]
        Array[String] annotations_all = [ "SVLEN","SUPP_SNIFFLES","SUPP_PBSV","SUPP_PAV","BIN_POS","BIN_POINT_MAPQ","BIN_POINT_SECONDARY","PL","PR","PL_PL_1","PL_PL_2","PL_PL_3","PL_PL_4","PL_PR_1","PL_PR_2","PL_PR_3","PL_PR_4","PR_PR_1","PR_PR_2","PR_PR_3","PR_PR_4",
                                          "FEX_DEPTH_RATIO","FEX_DEPTH_MAD","FEX_AB","FEX_CN_SLOP","FEX_MQ_DROP","FEX_CLIP_FRAC","FEX_SPLIT_READS","FEX_READ_LEN_MED","FEX_STRAND_BIAS","FEX_GC_FRAC","FEX_HOMOPOLYMER_MAX","FEX_LCR_MASK",
                                          "SNIFFLES_GT_COUNT","SNIFFLES_GQ","SNIFFLES_DR","SNIFFLES_DV",
                                          "CUTEFC_GT_COUNT","CUTEFC_GQ","CUTEFC_DR","CUTEFC_DV","CUTEFC_PL_1","CUTEFC_PL_2","CUTEFC_PL_3","CUTEFC_CIPOS_1","CUTEFC_CIPOS_2","CUTEFC_CILEN_1","CUTEFC_CILEN_2","CUTEFC_RE","CUTEFC_STRAND",
                                          "LRCALLER_GTCOUNT1_left","LRCALLER_AD11_left","LRCALLER_AD12_left","LRCALLER_AD13_left","LRCALLER_VA11_left","LRCALLER_VA12_left","LRCALLER_VA13_left","LRCALLER_PL11_left","LRCALLER_PL12_left","LRCALLER_PL13_left","LRCALLER_GTCOUNT2_left","LRCALLER_AD21_left","LRCALLER_AD22_left","LRCALLER_AD23_left","LRCALLER_VA21_left","LRCALLER_VA22_left","LRCALLER_VA23_left","LRCALLER_PL21_left","LRCALLER_PL22_left","LRCALLER_PL23_left","LRCALLER_GTCOUNT3_left","LRCALLER_AD31_left","LRCALLER_AD32_left","LRCALLER_AD33_left","LRCALLER_VA31_left","LRCALLER_VA32_left","LRCALLER_VA33_left","LRCALLER_PL31_left","LRCALLER_PL32_left","LRCALLER_PL33_left","LRCALLER_GTCOUNT4_left","LRCALLER_AD41_left","LRCALLER_AD42_left","LRCALLER_AD43_left","LRCALLER_VA41_left","LRCALLER_VA42_left","LRCALLER_VA43_left","LRCALLER_PL41_left","LRCALLER_PL42_left","LRCALLER_PL43_left","LRCALLER_GTCOUNT5_left","LRCALLER_AD51_left","LRCALLER_AD52_left","LRCALLER_AD53_left","LRCALLER_VA51_left","LRCALLER_VA52_left","LRCALLER_VA53_left","LRCALLER_PL51_left","LRCALLER_PL52_left","LRCALLER_PL53_left"
                                        ]
        
        File training_python_script
        File scoring_python_script
        File hyperparameters_json
        
        String docker_image = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb"
    }
    parameter_meta {
        input_vcf_gz: "Assumed to contain all the annotations used in this workflow."
        remote_outdir: "Without final slash"
        training_resource_bed: "The same BED used in `SV_Integration_Workpackage1`."
        hyperparameters_json: "Parameters for `gatk TrainVariantAnnotationsModel`."
    }
    
    call Score as score_custom {
        input:
            id = "custom",
            annotations = annotations_custom,
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            resource_vcf_gz = resource_vcf_gz,
            resource_vcf_gz_tbi = resource_vcf_gz_tbi,
            remote_outdir = remote_outdir,
            training_resource_bed = training_resource_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
    }
    call Score as score_fex {
        input:
            id = "fex",
            annotations = annotations_fex,
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            resource_vcf_gz = resource_vcf_gz,
            resource_vcf_gz_tbi = resource_vcf_gz_tbi,
            remote_outdir = remote_outdir,
            training_resource_bed = training_resource_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
    }
    call Score as score_sniffles {
        input:
            id = "sniffles",
            annotations = annotations_sniffles,
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            resource_vcf_gz = resource_vcf_gz,
            resource_vcf_gz_tbi = resource_vcf_gz_tbi,
            remote_outdir = remote_outdir,
            training_resource_bed = training_resource_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
    }
    call Score as score_cutefc {
        input:
            id = "cutefc",
            annotations = annotations_cutefc,
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            resource_vcf_gz = resource_vcf_gz,
            resource_vcf_gz_tbi = resource_vcf_gz_tbi,
            remote_outdir = remote_outdir,
            training_resource_bed = training_resource_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
    }
    call Score as score_lrcaller_left {
        input:
            id = "lrcaller_left",
            annotations = annotations_lrcaller_left,
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            resource_vcf_gz = resource_vcf_gz,
            resource_vcf_gz_tbi = resource_vcf_gz_tbi,
            remote_outdir = remote_outdir,
            training_resource_bed = training_resource_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
    }
    call Score as score_lrcaller_right {
        input:
            id = "lrcaller_right",
            annotations = annotations_lrcaller_right,
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            resource_vcf_gz = resource_vcf_gz,
            resource_vcf_gz_tbi = resource_vcf_gz_tbi,
            remote_outdir = remote_outdir,
            training_resource_bed = training_resource_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
    }
    call Score as score_lrcaller_all {
        input:
            id = "lrcaller_all",
            annotations = annotations_lrcaller_all,
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            resource_vcf_gz = resource_vcf_gz,
            resource_vcf_gz_tbi = resource_vcf_gz_tbi,
            remote_outdir = remote_outdir,
            training_resource_bed = training_resource_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
    }
    call Score as score_all {
        input:
            id = "all",
            annotations = annotations_all,
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            resource_vcf_gz = resource_vcf_gz,
            resource_vcf_gz_tbi = resource_vcf_gz_tbi,
            remote_outdir = remote_outdir,
            training_resource_bed = training_resource_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            training_python_script = training_python_script,
            scoring_python_script = scoring_python_script,
            hyperparameters_json = hyperparameters_json,
            docker_image = docker_image
    }
    
    output {
    }
}


# Performane on a 2-core, 16GB VM, all HPRC+HGSVC samples:
#
# TOOL                                  CPU            RAM              TIME
# ExtractVariantAnnotations                             8G               20s
# TrainVariantAnnotationsModel                        200M               10s
# ScoreVariantAnnotations                             800M               20s
#
task Score {
    input {
        String id
        
        File input_vcf_gz
        File input_vcf_gz_tbi
        File resource_vcf_gz
        File resource_vcf_gz_tbi
        String remote_outdir
        
        File? training_resource_bed
        File reference_fa
        File reference_fai

        Array[String] annotations
        File training_python_script
        File scoring_python_script
        File hyperparameters_json
        
        String docker_image
        Int n_cpu = 2
        Int ram_size_gb = 10
        Int disk_size_gb = 20
    }
    parameter_meta {
    }
    
    String docker_dir = "/root"
    
    command <<<
        set -euxo pipefail
        
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))
        export GATK_LOCAL_JAR="/root/gatk.jar"
        
        # 1. Ensuring that the input VCFs have the correct format
        bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z ~{input_vcf_gz} --output input_cleaned.vcf.gz
        bcftools index -f -t input_cleaned.vcf.gz
        rm -f ~{input_vcf_gz}
    
        bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z ~{resource_vcf_gz} --output resource_cleaned.vcf.gz
        bcftools index -f -t resource_cleaned.vcf.gz
        rm -f ~{resource_vcf_gz}
        
        # 2. Scoring
        EXCLUDE_CHROMOSOMES="-XL chr1 -XL chr2 -XL chr3 -XL chr4 -XL chr5"
        if ~{defined(training_resource_bed)}
        then
            BED_FLAG="-L ~{training_resource_bed}"
        else
            BED_FLAG=""
        fi
        gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ExtractVariantAnnotations -V input_cleaned.vcf.gz ${EXCLUDE_CHROMOSOMES} -O extract -A ~{sep=" -A " annotations} --resource:resource,training=true,calibration=true resource_cleaned.vcf.gz --maximum-number-of-unlabeled-variants 1000000000 --mode INDEL --mnp-type INDEL ${BED_FLAG}
        ls -laht
        # Output:
        # extract.annot.hdf5
        # extract.unlabeled.annot.hdf5
        # extract.vcf.gz
        # extract.vcf.gz.tbi
        gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" TrainVariantAnnotationsModel --annotations-hdf5 extract.annot.hdf5 --unlabeled-annotations-hdf5 extract.unlabeled.annot.hdf5 --model-backend PYTHON_SCRIPT --python-script ~{training_python_script} --hyperparameters-json ~{hyperparameters_json} -O train.train --mode INDEL --verbosity DEBUG
        ls -laht
        # Output: 
        # train.train.*
        gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ScoreVariantAnnotations -V input_cleaned.vcf.gz -O score -A ~{sep=" -A " annotations} --resource:resource,training=true,calibration=true resource_cleaned.vcf.gz --resource:extracted,extracted=true extract.vcf.gz --model-prefix train.train --model-backend PYTHON_SCRIPT --python-script ~{scoring_python_script} --mode INDEL --mnp-type INDEL --ignore-all-filters --verbosity DEBUG
        ls -laht
        # Output:
        # score.vcf.gz
        # score.vcf.gz.tbi
        # score.annot.hdf5
        # score.scores.hdf5
        
        gsutil -m mv score.vcf.gz ~{remote_outdir}/~{id}_score.vcf.gz
        gsutil -m mv score.vcf.gz.tbi ~{remote_outdir}/~{id}_score.vcf.gz.tbi
    >>>
    
    output {
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
}
