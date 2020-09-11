#!/bin/bash
set -e

CONFIG="../../config/experiment0_prelim/virus.yml"
CLUSTER="../../config/experiment2_AAV-OTC/cluster.json"
OUTPATH="../../out/experiment0_prelim/"
NAME="viruses"

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../snakemake_sim_analysis 

mkdir -p ${OUTPATH}
snakemake --configfile ${CONFIG} --snakefile combined_snakefile --dag | dot -Tsvg > ${OUTPATH}/${NAME}.dag.svg

snakemake \
 --snakefile combined_snakefile \
 --configfile ${CONFIG}\
 --jobs 50 \
 --use-singularity \
 --profile slurm \
 --cluster-config ${CLUSTER}
 

