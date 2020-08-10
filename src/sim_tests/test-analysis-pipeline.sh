#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../../intvi_pipeline 
#snakemake --configfile ../config/test/detection.yml --cores 1 --use-singularity --dag | dot -Tsvg > ../out/test_dag.svg
snakemake --configfile ../config/test/detection.yml --jobs 50 --use-singularity  --profile slurm --rerun-incomplete
