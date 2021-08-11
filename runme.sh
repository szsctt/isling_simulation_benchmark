#!/bin/bash

set -euo pipefail

# dependencies for running
conda list -n sim_isling || conda create -n sim_isling -c bioconda -c conda-forge entrez-direct=13.9 snakemake=6.6

# download references
bash src/references/make_refs.sh


