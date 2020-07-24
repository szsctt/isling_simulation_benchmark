#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../../intvi_pipeline 
snakemake --configfile ../config/test/detection.yml --cores 1 --use-singularity
