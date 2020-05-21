#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate simvi

SCRIPT="../intvi_simulation/insert_virus.py"
HOST="../data/references/test_host.fa"
VIRUS="../data/references/test_virus.fa"

OUTDIR="../out/test/"
mkdir -p ${OUTDIR}
INTS="${OUTDIR}integrations.fa"
LOCS="${OUTDIR}locs.txt"
HOSTINTS="${OUTDIR}host_locs.txt"

srun -N 1 python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS} \
 --locs ${LOCS} \
 --host_ints test_host_ints.csv \
 --int_num 1 