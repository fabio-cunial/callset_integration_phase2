#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l InterCenterBench.wdl
java -jar ${WOMTOOL_PATH} validate -l InterCenterMerge.wdl
java -jar ${WOMTOOL_PATH} validate -l PasteGTs.wdl
java -jar ${WOMTOOL_PATH} validate -l KanpigMerged.wdl
java -jar ${WOMTOOL_PATH} validate -l RemoveSamples.wdl
java -jar ${WOMTOOL_PATH} validate -l TruvariIntersamplePhase2.wdl
java -jar ${WOMTOOL_PATH} validate -l HGSVC3Dipcall2BAMs.wdl
java -jar ${WOMTOOL_PATH} validate -l Split.wdl
java -jar ${WOMTOOL_PATH} validate -l HGSVC3ExtractHapsFromAssemblies.wdl
java -jar ${WOMTOOL_PATH} validate -l DipcallPhase2.wdl
java -jar ${WOMTOOL_PATH} validate -l FilterIntrasampleDevPhase2.wdl
java -jar ${WOMTOOL_PATH} validate -l Kanpig.wdl
java -jar ${WOMTOOL_PATH} validate -l TruvariIntrasample.wdl
java -jar ${WOMTOOL_PATH} validate -l Resolve.wdl
java -jar ${WOMTOOL_PATH} validate -l PAV2SVs.wdl
