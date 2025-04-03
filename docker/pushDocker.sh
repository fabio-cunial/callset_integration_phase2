#!/bin/bash
#
TAG=""
cp ../scripts/*.java .
docker build --progress=plain -t fcunial/callset_integration_phase2 .
docker push fcunial/callset_integration_phase2${TAG}
rm -f *.java *.class