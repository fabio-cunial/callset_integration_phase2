#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l MapR10Phase2Scattered.wdl
java -jar ${WOMTOOL_PATH} validate -l AddReadGroup.wdl
java -jar ${WOMTOOL_PATH} validate -l MapCCSPhase2Prime.wdl
java -jar ${WOMTOOL_PATH} validate -l MapR10Phase2.wdl
java -jar ${WOMTOOL_PATH} validate -l DownloadAssembly.wdl
java -jar ${WOMTOOL_PATH} validate -l SubsampleSimple.wdl
java -jar ${WOMTOOL_PATH} validate -l ReadLengthDistribution.wdl
java -jar ${WOMTOOL_PATH} validate -l GetLongCalls.wdl
java -jar ${WOMTOOL_PATH} validate -l CheckDeNovo.wdl
java -jar ${WOMTOOL_PATH} validate -l BenchHprcSamples.wdl
java -jar ${WOMTOOL_PATH} validate -l GetPresentCalls.wdl
java -jar ${WOMTOOL_PATH} validate -l BenchCohortSamples.wdl
java -jar ${WOMTOOL_PATH} validate -l InterCenterMerge.wdl
java -jar ${WOMTOOL_PATH} validate -l CheckMendelian.wdl
java -jar ${WOMTOOL_PATH} validate -l PlotHwe.wdl
java -jar ${WOMTOOL_PATH} validate -l QcPlots.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage13.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage12.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage11.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage10.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage9.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage8.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage7.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage6.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage5.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage4.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage3.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage2.wdl
java -jar ${WOMTOOL_PATH} validate -l Workpackage1.wdl
java -jar ${WOMTOOL_PATH} validate -l InterCenterBench.wdl
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
