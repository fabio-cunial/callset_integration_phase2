#!/bin/bash
#
set -euxo pipefail

if [ $# -eq 0 ]; then
  echo "Need to provide tag for the docker image."
  exit 1
fi

TAG="$1"

docker build --progress=plain -t fcunial/callset_integration_phase2_resolve .
docker tag fcunial/callset_integration_phase2_resolve fcunial/callset_integration_phase2_resolve:${TAG}
docker push fcunial/callset_integration_phase2_resolve:${TAG}
