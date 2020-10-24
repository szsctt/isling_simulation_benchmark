#!/bin/bash
set -e

CONFIG="../../config/experiment1_OTC_chr1/read-properties.yml"
CLUSTER="../../config/experiment1_OTC_chr1/cluster.json"
OUTPATH="../../out/experiment1_OTC_chr1/"
NAME="read-properties"

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
 --cluster-config ${CLUSTER} \
 --rerun-incomplete -n
 

