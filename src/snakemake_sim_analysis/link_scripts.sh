#!/bin/bash
set -e

# make symbolic links 

SIM="../../intvi_simulation"
ANALYSIS="../../intvi_pipeline"

# simulation stuff
mkdir -p snakemake_rules
mkdir -p scripts
find "$SIM/snakemake_rules" -name '*smk' -exec ln -sf  "$(realpath '{}')" "$(realpath snakemake_rules)" \;
ln -sf "$SIM/snakemake_rules/make_df.py" "$(realpath snakemake_rules)"
find "$SIM/scripts" -name '*' -exec ln -sf "$(realpath '{}')" "$(realpath scripts)" \;

# analysis stuff
find "$ANALYSIS/snakemake_rules" -name '*smk' -exec ln -sf "$(realpath '{}')" "$(realpath snakemake_rules)" \;
find "$SIM/scripts" -name '*.pl' -exec ln -sf "$(realpath '{}')" "$(realpath scripts)" \;
find "$SIM/scripts" -name '*.pm' -exec ln -sf "$(realpath '{}')" "$(realpath scripts)" \;
find "$SIM/scripts" -name '*.R' -exec ln -sf "$(realpath '{}')" "$(realpath scripts)" \;
find "$SIM/scripts" -name '*.sh' -exec ln -sf "$(realpath '{}')" "$(realpath scripts)" \;
find "$SIM/scripts" -name 'post' -exec ln -sf "$(realpath '{}')" "$(realpath scripts)" \;

# python stuff
mkdir -p python_scripts
find "$SIM/snakemake_rules" -name 'parse_config.py' -exec ln -sf "$(realpath '{}')" "$(realpath python_scripts)" \;
find "$ANALYSIS/snakemake_rules" -name 'make_df.py' -exec ln -sf "$(realpath '{}')" "$(realpath python_scripts)" \;

# conda envs
mkdir -p envs
find "$SIM/envs" -name '*.yml' -exec ln -sf "$(realpath '{}')" "$(realpath envs)" \;
find "$ANALYSIS/envs" -name '*.yml' -exec ln -sf "$(realpath '{}')" "$(realpath envs)" \;
