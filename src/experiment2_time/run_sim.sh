#!/bin/bash

# ./run_sim.sh <config yml> <local>

CONFIG=$1

WD=$(pwd)

if [ -z "$2" ]
	then
		eval "$(conda shell.bash hook)"
		conda activate snakemake
		module load singularity
fi

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
 
 snakemake \
 --configfile ${CONFIG} \
 --jobs 100 \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --profile slurm \
 --rerun-incomplete \
 --latency-wait 120 \
 --notemp \
 --until write_summary

