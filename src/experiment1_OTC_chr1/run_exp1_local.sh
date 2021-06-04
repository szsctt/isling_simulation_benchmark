#!/bin/bash
set -e



CONFIG="../../config/experiment1_OTC_chr1/AAV-OTC.yml"
bash ./run-local.sh $CONFIG

CONFIG="../../config/experiment1_OTC_chr1/OTC-condition-breakdown.yml"

bash ./run-local.sh $CONFIG


