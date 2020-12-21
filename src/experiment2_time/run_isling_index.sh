#!/bin/bash
set -e

# ./run_sim.sh <sim config yml>  <cluster_config>

CONFIG=$1
CLUSTER=$2

WD=$(pwd)

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../../intvi_pipeline

# make config file from simulation config
CONFIGFOLDER=$(dirname $CONFIG)
IDXCONFIG="${CONFIGFOLDER}/auto/isling_index.yml"
mkdir -p $(dirname $IDXCONFIG)
python3 ${WD}/make_isling_index_config.py ${CONFIG} ${IDXCONFIG}

# run snakemake until index rule
snakemake \
 --configfile "${IDXCONFIG}" \
 --jobs 5 \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --profile slurm \
 --rerun-incomplete \
 --latency-wait 120 \
 --until index \
 --cluster-config ${CLUSTER}


