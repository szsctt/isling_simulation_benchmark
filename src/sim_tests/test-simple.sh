#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate simvi

SCRIPT="../../intvi_simulation/insert_virus_simple.py"
HOST="../../data/references/test_host.fa"
VIRUS="../../data/references/test_virus.fa"

OUTDIR="../../out/test/"
mkdir -p ${OUTDIR}


i=0
INTS="${OUTDIR}integrations_${i}.fa"
LOCS="${OUTDIR}ints_info_${i}.txt"
EPI="${OUTDIR}epi_info_${i}.txt"

python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS} \
 --int_info ${LOCS} \
 --int_num 10 \
 --seed 1 \
 --epi_num 5 \
 --epi_info ${EPI} \
 --verbose


