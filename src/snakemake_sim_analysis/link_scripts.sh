#!/bin/bash

# make symbolic links 

SIM="../../intvi_simulation"
ANALYSIS="../../intvi_pipeline"

# simulation stuff
mkdir -p snakemake_rules
ln -sf $(realpath "$SIM/snakemake_rules/"*smk) $(realpath snakemake_rules)
ln -sf "$SIM/snakemake_rules/make_df.py" snakemake_rules
ln -sf "$SIM/snakemake_rules/make_df.py" snakemake_rules
ln -sf "$SIM/scripts" .

# analysis stuff
ln -sf $(realpath "$ANALYSIS/snakemake_rules/"*smk) $(realpath snakemake_rules)
ln -sf $ANALYSIS/*.pl .
ln -sf $ANALYSIS/*.pm .
ln -sf $ANALYSIS/*.R .
ln -sf $ANALYSIS/*.sh .
ln -sf $ANALYSIS/post .

# python stuff
mkdir -p python_scripts
ln -sf $(realpath "$SIM/snakemake_rules/parse_config.py") $(realpath python_scripts)
ln -sf $(realpath "$ANALYSIS/snakemake_rules/make_df.py") $(realpath python_scripts)
