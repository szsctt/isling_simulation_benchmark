#!/bin/bash
set -e

# make symbolic links 

SIM="../../intvi_simulation"
ANALYSIS="../../intvi_pipeline"

# simulation stuff
mkdir -p snakemake_rules
find "$SIM/snakemake_rules" -name '*smk' -exec ln -sf  "$(realpath '{}')" "$(realpath snakemake_rules)" \;
ln -sf "$SIM/snakemake_rules/make_df.py" "$(realpath snakemake_rules)"
ln -sf "$SIM/scripts" .

# analysis stuff
find "$ANALYSIS/snakemake_rules" -name '*smk' -exec ln -sf "$(realpath '{}')" "$(realpath snakemake_rules)" \;
ln -sf $ANALYSIS/*.pl .
ln -sf $ANALYSIS/*.pm .
ln -sf $ANALYSIS/*.R .
ln -sf $ANALYSIS/*.sh .
ln -sf $ANALYSIS/post .

# python stuff
mkdir -p python_scripts
find "$SIM/snakemake_rules" -name 'parse_config.py' -exec ln -sf "$(realpath '{}')" "$(realpath python_scripts)" \;
find "$ANALYSIS/snakemake_rules" -name 'make_df.py' -exec ln -sf "$(realpath '{}')" "$(realpath python_scripts)" \;

# conda envs
mkdir -p envs
find "$SIM/envs" -name '*.yml' -exec ln -sf "$(realpath '{}')" "$(realpath envs)" \;
find "$ANALYSIS/envs" -name '*.yml' -exec ln -sf "$(realpath '{}')" "$(realpath envs)" \;
