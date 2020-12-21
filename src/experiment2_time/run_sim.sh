#!/bin/bash

# ./run_sim.sh <config yml> 

CONFIG=$1

WD=$(pwd)

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../../intvi_simulation

snakemake \
 --configfile ${CONFIG} \
 --jobs 100 \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --profile slurm \
 --rerun-incomplete \
 --latency-wait 120 \
 --notemp \
 --until art

