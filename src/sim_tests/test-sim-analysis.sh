#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

#cd ../snakemake_sim_analysis 
#snakemake --configfile ../../config/test/sim_and_detect.yml --cores 1 --use-singularity

cd ../../intvi_simulation 
echo "running simulation..."
snakemake --configfile ../config/test/simulation.yml --jobs 50 --use-singularity --profile slurm --rerun-incomplete

cd ../intvi_pipeline 
echo ""
echo "running analysis..."
snakemake --configfile ../config/test/detection.yml --jobs 50 --use-singularity --profile slurm --rerun-incomplete --cluster-config ../config/test/detection-cluster.json
