#!/bin/bash
set -e

# ./run-cloud.sh <config yml>

CONFIG=$1

cd ../snakemake_sim_analysis 

# need to strip "../.." from input and output file paths in config file
filename="$(dirname $CONFIG)/$(basename $CONFIG .yml)"
NEW="${filename}.cloud.yml"
echo $NEW
sed 's%../../data%data%g' $CONFIG > $NEW
sed -i 's%../../out%out%g' $NEW

snakemake \
 --snakefile combined_snakefile \
 --configfile ${NEW}\
 --cores 8 \
 --use-singularity \
 --rerun-incomplete \
 --default-remote-provider "GS" \
 --default-remote-prefix sjs_snakemake_test


