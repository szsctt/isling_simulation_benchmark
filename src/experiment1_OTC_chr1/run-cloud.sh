#!/bin/bash
set -e

# ./run.sh <config yml> <cluster json> <dag outpath> <dag name>

CONFIG=$1
CLUSTER=$2
OUTPATH=$3
NAME=$4

cd ../snakemake_sim_analysis 

mkdir -p ${OUTPATH}

# need to strip "../.." from input and output file paths in config file
filename="${$CONFIG%.*}"
NEW="${filename}.cloud.yml"
sed 's%../../data%data%g' $CONFIG > $NEW
sed -i 's%../../out%out%g' $NEW

snakemake \
 --snakefile combined_snakefile \
 --configfile ${NEW}\
 --cores 32 \
 --use-singularity \
 --rerun-incomplete \
 --default-remote-provider "GS" \
 --default-remote-prefix "sjs_snakemake_test"


