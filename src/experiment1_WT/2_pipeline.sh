#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate snakemake
module load singularity

cd ../../intvi_pipeline
snakemake \
--configfile ../config/experiment1/detection.yml \
 --jobs 50 \
 --use-singularity \
 --profile slurm \
 --cluster-config ../config/experiment1/cluster.json
