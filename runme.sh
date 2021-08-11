#!/bin/bash

set -euo pipefail

# dependencies for running
conda list -n sim_isling || conda create -n sim_isling -c bioconda -c conda-forge entrez-direct=13.9 snakemake=6.6




# download references
bash src/references/make_refs.sh

# link necessary scripts in combined snakefile directory
cd src/snakemake_sim_analysis
bash link_scripts.sh
cd ../..

# simulate data and analyse for AAV and OTC
bash src/experiment1_OTC_chr1/run_exp1_local.sh

# generate figures and tables
