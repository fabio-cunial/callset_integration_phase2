#!/bin/bash
#
BUCKET="???"

set -euxo pipefail

rm -f tasks.csv

echo "15x,del,38f74024-5300-49c1-8dbb-bf80e49873fe" >> tasks.csv
echo "15x,ins,a80b08ad-95fd-467e-88de-d68ad962e426" >> tasks.csv
echo "15x,insdup,a5eb15bb-56d5-4459-9f1a-5ff83781c666" >> tasks.csv
echo "15x,dup,003156c8-eebf-4d74-9366-8601c5fcfb13" >> tasks.csv
echo "15x,inv,13dc0b19-4804-435e-b428-00e802331e24" >> tasks.csv

echo "30x,del,490148d7-6bea-4724-84dc-2c5d9c500b15" >> tasks.csv
echo "30x,ins,768785b4-5bf9-450a-8d66-ca21b84651da" >> tasks.csv
echo "30x,insdup,6e060f58-f70a-4b88-a003-7dd5260f8a21" >> tasks.csv
echo "30x,dup,79cfed04-c941-4d40-909b-d18d20ba894f" >> tasks.csv
echo "30x,inv,6defd175-a8ad-4ad3-ab97-ddcdf0be787f" >> tasks.csv

while read LINE; do
    COVERAGE=$( echo ${LINE} | cut -d ',' -f 1 )
    SVTYPE=$( echo ${LINE} | cut -d ',' -f 2 )
    SUBMISSION_ID=$( echo ${LINE} | cut -d ',' -f 3 )
    gcloud storage cp ${BUCKET}/submissions/${SUBMISSION_ID}/SV_Integration_UltralongScore/'*'/call-score_all/train.train.indel.scorer.pkl "${BUCKET}/v3/${COVERAGE}/workpackage_1_ultralong_annotated/merged/scored_final_model"/${SVTYPE}.indel.scorer.pkl
    gcloud storage cp ${BUCKET}/submissions/${SUBMISSION_ID}/SV_Integration_UltralongScore/'*'/call-score_all/train.train.indel.calibrationScores.hdf5 "${BUCKET}/v3/${COVERAGE}/workpackage_1_ultralong_annotated/merged/scored_final_model"/${SVTYPE}.indel.calibrationScores.hdf5
done < tasks.csv
rm -f tasks.csv
gcloud storage ls -l "${BUCKET}/v3/15x/workpackage_1_ultralong_annotated/merged/scored_final_model"
gcloud storage ls -l "${BUCKET}/v3/30x/workpackage_1_ultralong_annotated/merged/scored_final_model"
