#!/bin/bash
set -e

# make symbolic links 

SIM=../../intvi_simulation
ANALYSIS=../../intvi_pipeline
TOOLS=../../intvi_other-tools

LNK="-sf"

# simulation stuff
mkdir -p snakemake_rules
mkdir -p scripts/post

find $SIM/snakemake_rules -name "*smk" -exec ln $LNK  $(realpath {}) $(realpath snakemake_rules/) \;
find $SIM/scripts -name "*" -exec ln $LNK $(realpath {}) $(realpath scripts/) \;


# analysis stuff
find $ANALYSIS/snakemake_rules -name "*smk" -exec ln $LNK $(realpath {}) $(realpath snakemake_rules/) \;
find $ANALYSIS/scripts -name "*.pl" -exec ln $LNK $(realpath {}) $(realpath scripts/) \;
find $ANALYSIS/scripts -name "*.pm" -exec ln $LNK $(realpath {}) $(realpath scripts/) \;
find $ANALYSIS/scripts -name "*.R" -exec ln $LNK $(realpath {}) $(realpath scripts/) \;
find $ANALYSIS/scripts -name "*.sh" -exec ln $LNK $(realpath {}) $(realpath scripts/) \;
find $ANALYSIS/scripts/post -name "*.R" -exec ln $LNK $(realpath {}) $(realpath scripts/post/) \;

# other tools stuff
find $TOOLS/scripts -name "*py" -exec ln $LNK $(realpath {}) $(realpath scripts/) \;
find $TOOLS/snakemake_rules -name "*smk" -exec ln $LNK $(realpath {}) $(realpath snakemake_rules/) \;

# python stuff
mkdir -p python_scripts
find $SIM/snakemake_rules -name "parse_config.py" -exec ln $LNK $(realpath {}) $(realpath python_scripts/) \;
find $ANALYSIS/snakemake_rules -name "make_df.py" -exec ln $LNK $(realpath {}) $(realpath python_scripts/) \;

# conda envs
mkdir -p envs
find $SIM/envs -name "*.yml" -exec ln $LNK $(realpath {}) $(realpath envs/) \;
find $ANALYSIS/envs -name "*.yml" -exec ln $LNK $(realpath {}) $(realpath envs/) \;
