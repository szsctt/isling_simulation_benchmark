#!/bin/bash

# ./run-cloud.sh <config yml> <cores>

CONFIG=$1
CORES=$2

cd src/snakemake_sim_analysis

snakemake \
 --snakefile combined_snakefile \
 --configfile ${CONFIG} \
 --cores ${CORES} \
 --jobs 1 \
 --restart-times 3 \
 --keep-going \
 --use-singularity \
 --scheduler greedy \
 --singularity-args '-B $(realpath ../../..)' \
 --rerun-incomplete


