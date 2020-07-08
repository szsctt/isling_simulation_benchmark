#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate snakemake

cd ../../intvi_pipeline
snakemake --configfile ../config/experiment1/detection.yml --cores 1 --use-conda