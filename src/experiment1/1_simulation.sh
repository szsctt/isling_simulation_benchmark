#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate simvi

cd ../../intvi_simulation 
snakemake --configfile ../config/experiment1/simulation.yml --cores 1