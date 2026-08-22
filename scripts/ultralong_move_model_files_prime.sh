#!/bin/bash
#
REMOTE_BUCKET_URI="gs://fc-7f861a33-ddb4-4b2f-8d10-5679c9df6108"
REMOTE_DESTINATION_URI_15x="${REMOTE_BUCKET_URI}/v3/15x/workpackage_1_ultralong_annotated_fourth_attempt/merged/scored_final_model"
REMOTE_DESTINATION_URI_30x="${REMOTE_BUCKET_URI}/v3/30x/workpackage_1_ultralong_annotated_fourth_attempt/merged/scored_final_model"

# SVTYPES=("del" "dup" "insdup" "inv" "ins" "bnd")
# SUBMISSION_IDS_15x=("985bdce3-43d2-4402-bd71-fd2a389c6cfb" "3dee98f5-d4ed-4adf-87a7-060f061e973b" "99bfa26a-fbc9-4887-b849-4cfa2a23bb3d" "626fefd4-4c46-4f28-9e89-a1ecd0b82ae1" "cba78ff7-5920-4127-91e3-6dec5d6a4a0f" "05c84a0f-013e-4657-b5ac-d7b9e21ac7b4")
# SUBMISSION_IDS_30x=("d0ba589e-7819-467b-bdac-33ebd2d45144" "1e4eaea8-ee25-4d66-b3b2-56ab86b9472b" "12009e2b-13c8-4f2e-9cc8-b5fd228c2be9" "49a21dbc-5e11-4e45-bf34-988c90e78383" "1be023c6-2346-4aeb-a0fd-0ac0ffdb1715" "3fb02236-56fd-4dc1-bf5d-9590df1d059b")

# Temporary code to move just BND models after realizing missing C1L_C2R_3_0 annotation
SVTYPES=("bnd")
N_TYPES=1
SUBMISSION_IDS_15x=("7c77b13d-1059-40db-adb7-69b09155d08c")
SUBMISSION_IDS_30x=("06a874c3-dfa7-47c1-b029-195f2b24522c")

set -euxo pipefail


for i in $(seq 0 $(( ${N_TYPES} - 1 ))); do
    SOURCE_URI="${REMOTE_BUCKET_URI}/submissions/${SUBMISSION_IDS_15x[${i}]}/SV_Integration_UltralongScore"
    gcloud storage cp ${SOURCE_URI}/'*'/call-score_all/train.train.indel.calibrationScores.hdf5 ${REMOTE_DESTINATION_URI_15x}/${SVTYPES[${i}]}.indel.calibrationScores.hdf5
    gcloud storage cp ${SOURCE_URI}/'*'/call-score_all/train.train.indel.scorer.pkl ${REMOTE_DESTINATION_URI_15x}/${SVTYPES[${i}]}.indel.scorer.pkl
done
for i in $(seq 0 $(( ${N_TYPES} - 1 ))); do
    SOURCE_URI="${REMOTE_BUCKET_URI}/submissions/${SUBMISSION_IDS_30x[${i}]}/SV_Integration_UltralongScore"
    gcloud storage cp ${SOURCE_URI}/'*'/call-score_all/train.train.indel.calibrationScores.hdf5 ${REMOTE_DESTINATION_URI_30x}/${SVTYPES[${i}]}.indel.calibrationScores.hdf5
    gcloud storage cp ${SOURCE_URI}/'*'/call-score_all/train.train.indel.scorer.pkl ${REMOTE_DESTINATION_URI_30x}/${SVTYPES[${i}]}.indel.scorer.pkl
done
