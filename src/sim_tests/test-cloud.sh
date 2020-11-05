#!/bin/bash
set -e

cd ../snakemake_sim_analysis 

snakemake \
 --snakefile combined_snakefile \
 --configfile ../../config/test/cloud_test.yml\
 --cores 2 \
 --use-singularity \
 --rerun-incomplete \
 --default-remote-provider "GS" \
 --default-remote-prefix "sjs_snakemake_test"

