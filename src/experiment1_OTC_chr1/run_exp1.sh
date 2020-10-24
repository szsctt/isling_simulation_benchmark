#!/bin/bash
set -e

CONFIG="../../config/experiment1_OTC_chr1/easier-harder.yml"
CLUSTER="../../config/experiment1_OTC_chr1/cluster.json"
OUTPATH="../../out/experiment1_OTC_chr1/"
NAME="easier-harder"

bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/read-properties.yml"
CLUSTER="../../config/experiment1_OTC_chr1/cluster.json"
OUTPATH="../../out/experiment1_OTC_chr1/"
NAME="read-properties"

bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/condition-breakdown.yml"
CLUSTER="../../config/experiment1_OTC_chr1/cluster.json"
OUTPATH="../../out/experiment1_OTC_chr1/"
NAME="condition-breakdown"

bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

