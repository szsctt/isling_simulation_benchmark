#!/bin/bash
set -e

CONFIG="../../config/experiment1_OTC_chr1/analysis-conditions.yml"
CLUSTER="../../config/experiment1_OTC_chr1/cluster.json"
OUTPATH="../../out/experiment1_OTC_chr1/"
NAME="analysis-conditions"

bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

