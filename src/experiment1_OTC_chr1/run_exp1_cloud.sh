#!/bin/bash
set -e

CONFIG="../../config/experiment1_OTC_chr1/analysis-conditions.yml"
CLUSTER="../../config/experiment1_OTC_chr1/cluster.json"
OUTPATH="../../out/experiment1_OTC_chr1/"
NAME="analysis-conditions"

bash ./run-cloud.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/easier-harder.yml"
NAME="easier-harder"

bash ./run-cloud.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/read-properties.yml"
NAME="read-properties"

bash ./run-cloud.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/condition-breakdown.yml"
NAME="condition-breakdown"

bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/chromosomes.yml"
NAME="chromosomes"

bash ./run-cloud.sh $CONFIG $CLUSTER $OUTPATH $NAME


