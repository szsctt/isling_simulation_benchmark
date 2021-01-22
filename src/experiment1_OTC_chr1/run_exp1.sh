#!/bin/bash
set -e


CLUSTER="../../config/experiment1_OTC_chr1/cluster.json"
OUTPATH="../../out/experiment1_OTC_chr1/"

CONFIG="../../config/experiment1_OTC_chr1/easier-harder.yml"
NAME="easier-harder"

echo "running with config $CONFIG"
bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/condition-breakdown_OTC-harder-3.yml"
NAME="condition-breakdown_OTC-harder"

echo "running with config $CONFIG"
bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/condition-breakdown_OTC-harder-2.yml"
NAME="condition-breakdown_OTC-harder"

#echo "running with config $CONFIG"
#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME


CONFIG="../../config/experiment1_OTC_chr1/condition-breakdown_OTC-harder.yml"
NAME="condition-breakdown_OTC-harder"

#echo "running with config $CONFIG"
#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME


CONFIG="../../config/experiment1_OTC_chr1/condition-breakdown_2_mapq-15.yml"
NAME="condition-breakdown-2"

#echo "running with config $CONFIG"
#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME 

CONFIG="../../config/experiment1_OTC_chr1/condition-breakdown_1_mapq-15.yml"
NAME="condition-breakdown-1"

#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME 




CONFIG="../../config/experiment1_OTC_chr1/condition-breakdown_2.yml"
NAME="condition-breakdown-2"

#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME &
#PID_2=$!

CONFIG="../../config/experiment1_OTC_chr1/condition-breakdown_1.yml"
NAME="condition-breakdown-1"

#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME &
#PID_1=$!

#wait $PID_1 $PID_2

CONFIG="../../config/experiment1_OTC_chr1/condition-breakdown_OTC-harder.yml"
NAME="condition-breakdown_OTC-harder"

#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/analysis-conditions.yml"
NAME="analysis-conditions"

#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/read-properties.yml"
NAME="read-properties"

#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME


