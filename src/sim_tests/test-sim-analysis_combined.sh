#!/bin/bash
set -e

SIM="../../intvi_simulation"
ANALYSIS="../../intvi_pipeline"

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

cd ../snakemake_sim_analysis 

mkdir -p snakemake_rules
ln -sf $(realpath "$SIM/snakemake_rules/"*smk) $(realpath snakemake_rules)
ln -sf "$SIM/snakemake_rules/make_df.py" snakemake_rules
ln -sf "$SIM/snakemake_rules/make_df.py" snakemake_rules
ln -sf "$SIM/scripts" .

ln -sf $(realpath "$ANALYSIS/snakemake_rules/"*smk) $(realpath snakemake_rules)
ln -sf $ANALYSIS/*.pl .
ln -sf $ANALYSIS/*.R .
ln -sf $ANALYSIS/*.sh .
ln -sf $ANALYSIS/post .


# python stuff
mkdir -p python_scripts
ln -sf $(realpath "$SIM/snakemake_rules/parse_config.py") $(realpath python_scripts)
ln -sf $(realpath "$ANALYSIS/snakemake_rules/make_df.py") $(realpath python_scripts)

snakemake --snakefile combined_snakefile --configfile ../../config/test/sim_and_detect.yml --cores 1 --use-singularity --profile slurm


