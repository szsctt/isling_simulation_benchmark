#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../../intvi_simulation 
snakemake --configfile ../config/test/simulation.yml --cores 1 --use-singularity
