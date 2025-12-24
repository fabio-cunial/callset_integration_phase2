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
docker build --progress=plain -t fcunial/kanpig_simplified .
docker tag fcunial/kanpig_simplified fcunial/kanpig_simplified:${TAG}
docker push fcunial/kanpig_simplified:${TAG}
rm -f *.java *.class *.r
