#!/bin/bash
set -e

NAME="min-sep"
CONFIG="../../config/experiment0_short-refs/${NAME}.yml"
CLUSTER="../../config/experiment0_short-refs/cluster.json"
OUTPATH="../../out/experiment0_short-refs/"


eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../snakemake_sim_analysis 

mkdir -p ${OUTPATH}
snakemake --configfile ${CONFIG} --snakefile combined_snakefile --dag | dot -Tsvg > ${OUTPATH}/${NAME}.dag.svg

snakemake \
 --snakefile combined_snakefile \
 --configfile ${CONFIG}\
 --jobs 100 \
 --use-singularity \
 --profile slurm \
 --cluster-config ${CLUSTER} \
 --rerun-incomplete
 


