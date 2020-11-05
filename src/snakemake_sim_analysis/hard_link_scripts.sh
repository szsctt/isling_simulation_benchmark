#!/bin/bash
set -e

# make symbolic links 

SIM=../../intvi_simulation
ANALYSIS=../../intvi_pipeline

# simulation stuff
mkdir -p snakemake_rules
mkdir -p scripts
find $SIM/snakemake_rules -name "*smk" -exec ln  $(realpath {}) $(realpath snakemake_rules) \;
find $SIM/scripts -name "*" -exec ln $(realpath {}) $(realpath scripts) \;

# analysis stuff
mkdir -p scripts/post
find $ANALYSIS/snakemake_rules -name "*smk" -exec ln  $(realpath {}) $(realpath snakemake_rules) \;
find $ANALYSIS/scripts -name "*.pl" -exec ln  $(realpath {}) $(realpath scripts) \;
find $ANALYSIS/scripts -name "*.pm" -exec ln  $(realpath {}) $(realpath scripts) \;
find $ANALYSIS/scripts -name "*.R" -exec ln  $(realpath {}) $(realpath scripts) \;
find $ANALYSIS/scripts -name "*.sh" -exec ln  $(realpath {}) $(realpath scripts) \;
find $ANALYSIS/scripts/post -name "*.R" -exec ln  $(realpath {}) $(realpath scripts/post) \;

# python stuff
mkdir -p python_scripts
find $SIM/snakemake_rules -name "parse_config.py" -exec ln  $(realpath {}) $(realpath python_scripts) \;
find $ANALYSIS/snakemake_rules -name "make_df.py" -exec ln  $(realpath {}) $(realpath python_scripts) \;

# conda envs
mkdir -p envs
find $SIM/envs -name "*.yml" -exec ln  $(realpath {}) $(realpath envs) \;
find $ANALYSIS/envs -name "*.yml" -exec ln  $(realpath {}) $(realpath envs) \;
