#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../../intvi_simulation 
snakemake --configfile ../config/test/simulation.yml --jobs 10 --use-singularity --profile slurm --rerun-incomplete

