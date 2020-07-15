#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate snakemake
module load singularity

cd ../../intvi_simulation 
snakemake --configfile ../config/experiment1/simulation.yml --jobs 50 --use-singularity --profile slurm
