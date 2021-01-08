#!/bin/bash
set -e


eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../snakemake_sim_analysis 

mkdir -p ../../out/test
#snakemake --configfile ../../config/test/sim_and_detect.yml --snakefile combined_snakefile --dag | dot -Tsvg > ../../out/test/test_combined.dag.svg

#snakemake \
# --snakefile combined_snakefile \
# --configfile ../../config/test/sim_and_detect.yml\
# --jobs 50 \
# --use-singularity \
# --profile slurm \
# --rerun-incomplete
 
# for local execution
#snakemake \
# --snakefile combined_snakefile \
# --configfile ../../config/test/sim_and_detect.yml\
# --jobs 1 \
# --use-conda \
# --conda-frontend mamba \
# --rerun-incomplete

#srun -c5 --time 2:00:00 --mem 50gb \
snakemake \
 --snakefile combined_snakefile \
 --configfile ../../config/test/sim_and_detect.yml \
 --cores 1 \
 --jobs 1 \
 --use-singularity \
 --singularity-args '-B $(realpath ../..)' \
 --rerun-incomplete --until run_discordant --forceall -np

