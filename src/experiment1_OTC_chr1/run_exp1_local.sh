#!/bin/bash
set -e

#CONFIG="../../config/experiment1_OTC_chr1/analysis-conditions.yml"

#bash ./run-local.sh $CONFIG

CONFIG="../../config/experiment1_OTC_chr1/easier-harder.yml"

bash ./run-local.sh $CONFIG

CONFIG="../../config/experiment1_OTC_chr1/read-properties.yml"

bash ./run-local.sh $CONFIG

CONFIG="../../config/experiment1_OTC_chr1/condition-breakdown.yml"

bash ./run-local.sh $CONFIG

CONFIG="../../config/experiment1_OTC_chr1/chromosomes.yml"
NAME="chromosomes"

bash ./run-local.sh $CONFIG


