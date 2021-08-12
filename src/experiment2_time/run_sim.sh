#!/bin/bash

# ./run_sim.sh <config yml> <local>

CONFIG=$1
CORES=$2

WD=$(pwd)

if [ -z "$3" ]
	then
		eval "$(conda shell.bash hook)"
		conda activate snakemake
		module load singularity
		SNAKEARGS="--profile slurm --latency-wait 120 --jobs ${CORES}"
else
	SNAKEARGS="--resources mem_mb=20000 --cores ${CORES}"
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

