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
podman build -t fcunial/kanpig_subintv .
podman tag fcunial/kanpig_subintv fcunial/kanpig_subintv:${TAG}
podman push fcunial/kanpig_subintv:${TAG}
rm -f *.java *.class *.r
