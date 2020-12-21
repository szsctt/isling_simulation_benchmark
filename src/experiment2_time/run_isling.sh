#!/bin/bash
set -e

# ./run_isling.sh <sim config yml> <docker string>

CONFIG=$1
DOCKER=$2
OUT=$3

CONTAINER="isling.sif"
WD=$(pwd)

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../../intvi_pipeline

if [ ! -e ${CONTAINER} ]; then
	singularity pull --name ${CONTAINER} --force ${DOCKER}
fi

python3 ${WD}/launch_isling_jobs.py ${CONFIG} ${CONTAINER} "${WD}/make_isling_config.py"

