#!/bin/bash

# ./run_sim.sh <config yml> <local>

CONFIG=$1

WD=$(pwd)

if [ -z "$2" ]
	then
		eval "$(conda shell.bash hook)"
		conda activate snakemake
		module load singularity
		SNAKEARGS="--profile slurm --latency-wait 120 --jobs 100"
else
	SNAKEARGS="--resources mem_mb=60000 --cores 15"
fi

cd ../../intvi_simulation

snakemake \
 --configfile ${CONFIG} \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --rerun-incomplete \
 --notemp \
 --until art ${SNAKEARGS}
 
 snakemake \
 --configfile ${CONFIG} \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --rerun-incomplete \
 --notemp \
 --until write_summary ${SNAKEARGS}

