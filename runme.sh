#!/bin/bash

set -e

# set number of cores to use
CORES="15"

# dependencies for running
conda list -n sim_isling || mamba create -n sim_isling -c bioconda -c conda-forge entrez-direct>=11 snakemake=6.1 -y

eval "$(conda shell.bash hook)"
conda activate sim_isling

# download references
echo "downloading references"
bash src/references/make_refs.sh

# link necessary scripts in combined snakefile directory
cd src/snakemake_sim_analysis
bash link_scripts.sh
cd ../..

# simulate data and analyse for AAV and OTC
echo "simulating and analysing data"
bash src/experiment1_OTC_chr1/run_exp1_local.sh $CORES

# generate example rmd report
bash src/experiment1_OTC_chr1/example_rmd.sh

# generate figures and tables
echo "generating tables and figures"
bash src/experiment1_OTC_chr1/make_figures.sh

# do runtime stuff
cd src/experiment2_time
echo "checking runtime"
./run_time_local.sh $CORES
cd ../..

# make runtime figures
bash src/experiment2_time/make_figures.sh
