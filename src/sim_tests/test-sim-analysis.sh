#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../snakemake_sim_analysis 
snakemake --configfile ../../config/test/sim_and_detect.yml --cores 1 --use-singularity
