#!/bin/bash
set -e

CONFIG="../config/experiment2_time/test_sim.yml"
CLUSTER="../config/experiment2_time/cluster.json"
ISLING="docker://szsctt/isling:latest"

# simulate data
bash run_sim.sh ${CONFIG}

# bwa index host and virus for isling
bash run_isling_index.sh ${CONFIG} ${CLUSTER}

# run isling
bash run_isling.sh ${CONFIG} ${ISLING}

# run other tools

