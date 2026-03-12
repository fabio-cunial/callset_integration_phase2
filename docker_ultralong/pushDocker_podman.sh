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
cp ../scripts/*.py .
podman build -t us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong .
podman tag us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:${TAG}
podman push us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong:${TAG}
rm -f *.java *.class *.r *.py
