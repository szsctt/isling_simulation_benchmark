#!/bin/bash
set -e

CONFIG="../../config/experiment0_prelim/condition-breakdown.yml"
CLUSTER="../../config/experiment0_prelim/cluster.json"
OUTPATH="../../out/experiment0_prelim/"
NAME="condition-breakdown"

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../snakemake_sim_analysis 

mkdir -p ${OUTPATH}
snakemake --configfile ${CONFIG} --snakefile combined_snakefile --dag | dot -Tsvg > ${OUTPATH}/${NAME}.dag.svg

snakemake \
 --snakefile combined_snakefile \
 --configfile ${CONFIG}\
 --jobs 1 \
 --use-singularity \
 --profile slurm \
 --cluster-config ${CLUSTER}  
 

