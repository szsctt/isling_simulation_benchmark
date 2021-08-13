#!/bin/bash

set -e

CORES=$1

cd ../../

snakemake \
--configfile benchmark/simulated_data/config/experiment1_OTC_chr1/example_rmd.yml \
--jobs ${CORES} \
--use-singularity \
--rerun-incomplete
