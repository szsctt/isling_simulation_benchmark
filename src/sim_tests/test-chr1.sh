#!/bin/bash
set -e

cd ../snakemake_sim_analysis 

snakemake \
 --snakefile combined_snakefile \
 --configfile ../../config/test/test_chr1.yml\
 --jobs 100 \
 --use-singularity \
 --singularity-args "-B $(realpath ../..)" \
 --rerun-incomplete \
 --notemp --profile slurm --until vseq_toolkit
