#!/bin/bash
set -e

CONFIG="../../config/experiment2_AAV-OTC/episomes.yml"
CLUSTER="../../config/experiment2_AAV-OTC/cluster.json"
OUTPATH="../../out/experiment2_AAV-OTC/"
NAME="chomosomes"

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
 --cluster-config ${CLUSTER}
 

