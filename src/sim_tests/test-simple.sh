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
 --p_host_deletion 0.9 \
 --lambda_host_deletion 2 \
 --verbose
 
TESTHOSTSCRIPT="../../intvi_simulation/reconstruct_host.py"
python3 ${TESTHOSTSCRIPT} \
--int-fa ${INTS} \
--int-info ${LOCS} \
--host-fa ${HOST} \

printf "\n\n"

HOST="../../data/references/test_human.fa"
VIRUS="../../data/references/test_AAV.fa"

i=1
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
 --p_host_deletion 0.9 \
 --lambda_host_deletion 2 \
 --verbose
 
python3 ${TESTHOSTSCRIPT} \
--int-fa ${INTS} \
--int-info ${LOCS} \
--host-fa ${HOST} \