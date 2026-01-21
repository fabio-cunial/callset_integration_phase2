#!/bin/bash
#
set -euxo pipefail

if [ $# -eq 0 ]; then
  echo "Need to provide tag for the docker image."
  exit 1
fi

TAG="$1"

podman build -t fcunial/callset_integration_phase2_minimap2r10 .
podman tag fcunial/callset_integration_phase2_minimap2r10 fcunial/callset_integration_phase2_minimap2r10:${TAG}
podman push fcunial/callset_integration_phase2_minimap2r10:${TAG}
