#!/bin/bash
set -e

#CONFIG="../config/experiment2_time/simulation.yml"
#CONFIG="../config/experiment2_time/simulation_100x.yml"
CONFIG="../config/experiment2_time/test_sim.yml"
CLUSTER="../config/experiment2_time/cluster.json"
VIFI_REPO="../data/references/data_repo"

ISLING="docker://szsctt/isling:runtime"
POLYIDUS="docker://szsctt/polyidus:3"
SEEKSV="docker://szsctt/seeksv:1"
VIFI="docker://szsctt/vifi:1"
VSEQ="docker://szsctt/vseq:1"



# simulate data
bash run_sim.sh ${CONFIG} local

# bwa index host and virus for tools
bash run_index.sh ${CONFIG} ${CLUSTER} ${VIFI_REPO} local

# run tools
bash run_tools.sh ${CONFIG} ${ISLING} ${SEEKSV} ${POLYIDUS} ${VIFI} ${VIFI_REPO} ${VSEQ} 1 local


