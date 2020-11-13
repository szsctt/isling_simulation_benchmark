#!/bin/bash

# ./run-cloud.sh <config yml>

CONFIG=$1

cd ../snakemake_sim_analysis 

# need to strip "../.." from input and output file paths in config file
filename="$(dirname $CONFIG)/$(basename $CONFIG .yml)"
NEW="${filename}.cloud.yml"
sed 's%../../data%data%g' $CONFIG > $NEW
sed -i 's%../../out%out%g' $NEW

snakemake \
 --snakefile combined_snakefile \
 --configfile ${CONFIG} \
 --cores 8 \
 --jobs 3 \
 --restart-times 3 \
 --keep-going \
 --resources mem_mb=28000 \
 --use-singularity \
 --singularity-args '-B $(realpath ../..)' \
 --rerun-incomplete \
 --until vifi

# --configfile ${NEW}\
# --default-remote-provider "GS" \
# --default-remote-prefix sjs_snakemake_test


