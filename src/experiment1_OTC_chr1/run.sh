#!/bin/bash

# ./run.sh <config yml> <cluster json> <dag outpath> <dag name>

CONFIG=$1
CLUSTER=$2
OUTPATH=$3
NAME=$4

WD=$(pwd)

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../snakemake_sim_analysis 

#mkdir -p ${OUTPATH}
#snakemake --configfile ${CONFIG} --snakefile combined_snakefile --dag --until combine_ints | dot -Tsvg > ${OUTPATH}/${NAME}.dag.svg

snakemake \
 --snakefile combined_snakefile \
 --configfile ${CONFIG}\
 --jobs 100 \
 --use-singularity \
 --singularity-args "-B $(realpath ../../)" \
 --profile slurm \
 --rerun-incomplete \
 --latency-wait 120 \
 --scheduler greedy
 
# --cluster-config ${CLUSTER} \

