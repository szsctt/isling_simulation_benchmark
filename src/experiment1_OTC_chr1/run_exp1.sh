#!/bin/bash
set -e


CLUSTER="../../config/experiment1_OTC_chr1/cluster.json"
OUTPATH="../../out/experiment1_OTC_chr1/"

CONFIG="../../config/experiment1_OTC_chr1/AAV-OTC_parameter-optimise.yml"
NAME="easier-harder"

echo "running with config $CONFIG"
bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/AAV-OTC.yml"
NAME="easier-harder"

echo "running with config $CONFIG"
bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/OTC-condition-breakdown.yml"
NAME="condition-breakdown_OTC-harder"

echo "running with config $CONFIG"
bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

