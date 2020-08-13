#!/bin/bash
set -e

CONFIG="../../config/experiment1_WT/chromosomes.yml"
CLUSTER="../../config/experiment1_WT/cluster.json"
OUTPATH="../../out/experiment1_WT/"
NAME="chomosomes"

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../snakemake_sim_analysis 

mkdir -p ${OUTPATH}
snakemake --configfile ${CONFIG} --snakefile combined_snakefile --dag | dot -Tsvg > ${OUTPATH}/${name}.dag.svg

snakemake \
 --snakefile combined_snakefile \
 --configfile ${CONFIG}\
 --jobs 50 \
 --use-singularity \
 --profile slurm \
 --cluster-config ${CLUSTER}
 

