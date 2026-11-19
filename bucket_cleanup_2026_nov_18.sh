#!/bin/bash
#
# November 18, 2026 bucket cleanup.
#
UW_BUCKET="..."
HA_BUCKET="..."
BCM_BUCKET="..."
BI_BUCKET="..."
HPRC_Y2_BUCKET="..."

set -euxo pipefail




# ------------------------ SV integration filesystem ---------------------------
# For each workspace we show the directory structure before the cleanup.

## AoU_DRC_WGS_LongReads_UW_CCS
#
# /cunial_intesample_vcf/
#   /reference/
#   /training_resources/                    # v1 run, 12x and 25x.
#   /truvari_intersample/                   # v1 run, with 07 and 09 and chm13 versions.
#   /truvari_intersample_regenotyped/       # v1 run, with 07 and 09 and chm13 versions.
#   /v2/                                    # v2 run. Files not used by v3.
#       /training_resource/
#       /workpackage_1/
#       /workpackage_2/
#       /workpackage_3/
#       /workpackage_4/
#       /workpackage_9/
#   /v3/
#   /workpackages/                          # Only chunk files
# /cunial_svqc                              # Only small files
#
gcloud storage rm -r ${UW_BUCKET}/scratch/cunial_intersample_vcf/training_resources/
gcloud storage rm -r ${UW_BUCKET}/scratch/cunial_intersample_vcf/truvari_intersample/
gcloud storage rm -r ${UW_BUCKET}/scratch/cunial_intersample_vcf/truvari_intersample_regenotyped/
gcloud storage rm -r ${UW_BUCKET}/scratch/cunial_intersample_vcf/v2/


## AoU_DRC_WGS_LongReads_HA_v9 copy
#
# /cunial_intesample_vcf/
#   /reference/
#   /training_resources/                    # v1 run, 12x and 25x.
#   /truvari_intersample/                   # v1 run, with 07 and 09 versions, but only hg38.
#   /truvari_intersample_regenotyped/       # v1 run, with 07 and 09 versions, but only hg38.
#   /v3/
#   /workpackages/                          # Only chunk files
# /cunial_svqc                              # Only small files
# 
gcloud storage rm -r ${HA_BUCKET}/scratch/cunial_intersample_vcf/training_resources/
gcloud storage rm -r ${HA_BUCKET}/scratch/cunial_intersample_vcf/truvari_intersample/
gcloud storage rm -r ${HA_BUCKET}/scratch/cunial_intersample_vcf/truvari_intersample_regenotyped/


## AoU_DRC_WGS_LongReads_BCM_CCS
#
# /cunial_intesample_vcf/
#   /reference/
#   /training_resources/                    # v1 run, 12x and 25x.
#   /truvari_intersample/                   # v1 run, with 07 and 09 and chm13 versions.
#   /truvari_intersample_regenotyped/       # v1 run, with 07 and 09 and chm13 versions.
#   /v2/                                    # v2 run. Files not used by v3.
#       /training_resource/
#       /workpackage_1/
#       /workpackage_2/
#       /workpackage_3/
#       /workpackage_4/
#       /workpackage_9/
#   /v3/
#   /workpackages/                          # Only small files
# /cunial_svqc                              # Only small files
#
gcloud storage rm -r ${BCM_BUCKET}/scratch/cunial_intersample_vcf/training_resources/
gcloud storage rm -r ${BCM_BUCKET}/scratch/cunial_intersample_vcf/truvari_intersample/
gcloud storage rm -r ${BCM_BUCKET}/scratch/cunial_intersample_vcf/truvari_intersample_regenotyped/
gcloud storage rm -r ${BCM_BUCKET}/scratch/cunial_intersample_vcf/v2/


## AoU_DRC_WGS_LongReads_BI_CCS
#
# /cunial-read-length-distribution/                     # Old analysis of read length, not used any more.
# /cunial_intesample_vcf/
#   /04_scored_vcf/                                     # v1 run, large files.
#   /05_truvari_collapse_output_concatenated_07/        # v1 run
#   /05_truvari_collapse_output_concatenated_09/        # v1 run
#   /08_final_07/                                       # v1 run, large files.
#   /08_final_09/                                       # v1 run, large files.
#   /ancestries/                                        # Only small files
#   /kanpig_inter_tests/                                # Probably all the cohort re-genotyping experiments?
#       /ab/                                            # Large files
#       /ancestry_only/                                 # Large files
#       /beta_binomial/                                 # Large files
#       ... many more subdirs ...
#       /tmp/
#   /reference/
#   /tmp/                                               # PNGs and other tmp files
#   /training_resources/                                # v1 run, 12x and 25x.
#   /v2/
#       ... many subdirs e.g.: ...                      # Probably many analyses in the github discussions
#       /regenotyping_analysis_*/
#       /shapeit4_analysis/
#       /kanpig_hap_analysis_*/
#       /kanpig_trgt_analysis_*/
#       /phab_analysis/
#       /workpackage_*/
#   /v3/
#   /workpackages/                                      # Only small files
#   /workpackages_amr/                                  # Only small files
#   /workpackages_hprc/                                 # Only small files
#   /workpackages_mid/                                  # Only small files
# /cunial_svqc                                          # Only small files
#
gcloud storage rm -r ${BI_BUCKET}/scratch/cunial-read-length-distribution/
gcloud storage rm -r ${BI_BUCKET}/scratch/cunial_intersample_vcf/04_scored_vcf/
gcloud storage rm -r ${BI_BUCKET}/scratch/cunial_intersample_vcf/05_truvari_collapse_output_concatenated_07/
gcloud storage rm -r ${BI_BUCKET}/scratch/cunial_intersample_vcf/05_truvari_collapse_output_concatenated_09/
gcloud storage rm -r ${BI_BUCKET}/scratch/cunial_intersample_vcf/08_final_07/
gcloud storage rm -r ${BI_BUCKET}/scratch/cunial_intersample_vcf/08_final_09/
gcloud storage rm -r ${BI_BUCKET}/scratch/cunial_intersample_vcf/kanpig_inter_tests/
gcloud storage rm -r ${BI_BUCKET}/scratch/cunial_intersample_vcf/tmp/
gcloud storage rm -r ${BI_BUCKET}/scratch/cunial_intersample_vcf/training_resources/
gcloud storage rm -r ${BI_BUCKET}/scratch/cunial_intersample_vcf/v2/


## HPRC Y2
#
# /numerical_id/gsutil/tmp/parallel_composite_uploads/for_details_see/gsuitl_help_cp/hash_id        # Large files
# /tmp/                                                                                             # Mostly small files, not worth removing.
# /v2/
#   /15x/
#   /30x/
#   /training_resource/
# /v3/
# /workpackages/                                                                                    # Only small files
#
gcloud storage rm -r ${HPRC_Y2_BUCKET}'/*/gsutil/tmp/*'
gcloud storage rm -r ${HPRC_Y2_BUCKET}/v2/


# ------------------- Temporary files used for computing R^2 -------------------
# 160 GB in `rsquare_analysis/matrices`
# 130 GB in `rsquare_analysis/bubble_representation/matrices`
RSQUARE_DIR_1="${BI_BUCKET}/scratch/cunial_intersample_vcf/v3/rsquare_analysis/matrices"
RSQUARE_DIR_2="${BI_BUCKET}/scratch/cunial_intersample_vcf/v3/rsquare_analysis/bubble_representation/matrices"

gcloud storage rm ${RSQUARE_DIR_1}/'*.tsv.gz' 
gcloud storage rm ${RSQUARE_DIR_2}/'*.tsv.gz' 


# ------------ Temporary files used for computing smallest distance ------------
# Not worth deleting:
# "${BI_BUCKET}/scratch/cunial_intersample_vcf/v3/smallest_distance_analysis/matrices"  # 40 MB
# "${BI_BUCKET}/scratch/cunial_intersample_vcf/v3/smallest_distance_analysis_minaf_0_minlength_0/matrices"  # 40 MB
# "${BI_BUCKET}/scratch/cunial_intersample_vcf/v3/smallest_distance_analysis_minaf_0_minlength_50/matrices"  # 40 MB
# "${BI_BUCKET}/scratch/cunial_intersample_vcf/v3/smallest_distance_analysis_minaf_0_minlength_0/before_regenotyping/matrices"  # 40 MB
# "${BI_BUCKET}/scratch/cunial_intersample_vcf/v3/smallest_distance_analysis_minaf_0_minlength_50/before_regenotyping/matrices"  # 40 MB


# --------------- Temporary files used for computing numneighbors --------------
# Not worth deleting:
# "${BI_BUCKET}/scratch/cunial_intersample_vcf/v3/numneighbors_analysis_minlength_0"  # 150 MB
# "${BI_BUCKET}/scratch/cunial_intersample_vcf/v3/numneighbors_analysis_minlength_50"  # 75 MB
# "${BI_BUCKET}/scratch/cunial_intersample_vcf/v3/numneighbors_analysis_minlength_0/before_regenotyping"  # 170 MB
# "${BI_BUCKET}/scratch/cunial_intersample_vcf/v3/numneighbors_analysis_minlength_50/before_regenotyping" # 170 MB


# --------------------- Large UltralongScore WDL run outputs -------------------
# Manually removed the output files of every run of UltralongScore WDL in the HPRC_Y2 workspace.
