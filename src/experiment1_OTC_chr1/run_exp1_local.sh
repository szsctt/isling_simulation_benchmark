#!/bin/bash
set -e

# run AAV and OTC experiment
CONFIG="../../config/experiment1_OTC_chr1/AAV-OTC.yml"
bash src/experiment1_OTC_chr1/run-local.sh $CONFIG

# run OTC experiment
CONFIG="../../config/experiment1_OTC_chr1/OTC-condition-breakdown.yml"
bash src/experiment1_OTC_chr1/run-local.sh $CONFIG


