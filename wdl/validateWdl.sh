#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l TruvariIntrasample.wdl
java -jar ${WOMTOOL_PATH} validate -l Resolve.wdl
java -jar ${WOMTOOL_PATH} validate -l PAV2SVs.wdl
