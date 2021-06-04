#!/bin/bash

# ./run-cloud.sh <config yml>

CONFIG=$1

cd ../snakemake_sim_analysis 

snakemake \
 --snakefile combined_snakefile \
 --configfile ${CONFIG} \
 --cores 1 \
 --jobs 1 \
 --restart-times 3 \
 --keep-going \
 --resources mem_mb=20000 \
 --use-singularity \
 --singularity-args '-B $(realpath ../..)' \
 --rerun-incomplete 



