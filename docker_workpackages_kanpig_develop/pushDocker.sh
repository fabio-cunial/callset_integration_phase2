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
docker build --progress=plain -t fcunial/kanpig_develop .
docker tag fcunial/kanpig_develop fcunial/kanpig_develop:${TAG}
docker push fcunial/kanpig_develop:${TAG}
rm -f *.java *.class *.r *.py
