#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate simvi

cd ../../intvi_simulation 
snakemake --configfile ../config/test-sim.yml --cores 1