#!/bin/bash
#
set -euxo pipefail

if [ $# -eq 0 ]; then
  echo "Need to provide tag for the docker image."
  exit 1
fi

TAG="$1"

cp ../scripts/*.java .
docker build --progress=plain -t fcunial/callset_integration_phase2 .
docker tag fcunial/callset_integration_phase2 fcunial/callset_integration_phase2:${TAG}
docker push fcunial/callset_integration_phase2:${TAG}
rm -f *.java *.class