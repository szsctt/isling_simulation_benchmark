#!/bin/bash
set -e


CLUSTER="../../config/experiment1_OTC_chr1/cluster.json"
OUTPATH="../../out/experiment1_OTC_chr1/"

CONFIG="../../config/experiment1_OTC_chr1/AAV-OTC_isling-optimise.yml"
NAME="AAV-OTC_parameters"

#echo "running with config $CONFIG"
#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/OTC_isling-optimise_filtering.yml"
NAME="OTC_parameters"

#echo "running with config $CONFIG"
#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/AAV_isling-optimise_filtering.yml"
NAME="AAV_parameters"

#echo "running with config $CONFIG"
#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/AAV-OTC_isling-optimise_filtering.yml"
NAME="OTC-AAV_parameters"

#echo "running with config $CONFIG"
#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/AAV-OTC_parameter-optimise.yml"
NAME="AAV-OTC_parameters"

#echo "running with config $CONFIG"
#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/AAV-OTC.yml"
NAME="AAV-OTC"

echo "running with config $CONFIG"
bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/OTC-condition-breakdown.yml"
NAME="condition-breakdown_OTC"

echo "running with config $CONFIG"
bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME


CONFIG="../../config/experiment1_OTC_chr1/AAV-OTC_no-ints.yml"
NAME="AAV-OTC"

echo "running with config $CONFIG"
bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

CONFIG="../../config/experiment1_OTC_chr1/AAV-OTC_conditions.yml"
NAME="AAV-OTC_conditions"

#echo "running with config $CONFIG"
#bash ./run.sh $CONFIG $CLUSTER $OUTPATH $NAME

