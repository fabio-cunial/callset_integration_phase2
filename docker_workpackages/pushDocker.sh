#!/bin/bash
#
set -euxo pipefail

if [ $# -eq 0 ]; then
  echo "Need to provide tag for the docker image."
  exit 1
fi

TAG="$1"

cp ../scripts/*.java .
cp ../scripts/*.r .
docker build --progress=plain -t fcunial/callset_integration_phase2_workpackages .
docker tag fcunial/callset_integration_phase2_workpackages fcunial/callset_integration_phase2_workpackages:${TAG}
docker push fcunial/callset_integration_phase2_workpackages:${TAG}
rm -f *.java *.class *.r
