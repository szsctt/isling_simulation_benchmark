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
--host-fa ${HOST} 

TESTVIRUSCRIPT="../../intvi_simulation/reconstruct_virus.py"
python3 ${TESTVIRUSCRIPT} \
--int-fa ${INTS} \
--int-info ${LOCS} \
--virus-fa ${VIRUS} 



printf "\n\n"

HOST="../../data/references/test_human.fa"
VIRUS="../../data/references/test_AAV.fa"

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
--host-fa ${HOST} 

python3 ${TESTVIRUSCRIPT} \
--int-fa ${INTS} \
--int-info ${LOCS} \
--virus-fa ${VIRUS} 

# change seed
printf "\n\n"

python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS} \
 --int_info ${LOCS} \
 --int_num 10 \
 --seed 12345 \
 --epi_num 5 \
 --epi_info ${EPI} \
 --p_host_deletion 0.9 \
 --lambda_host_deletion 2 \
 --verbose
 
python3 ${TESTHOSTSCRIPT} \
--int-fa ${INTS} \
--int-info ${LOCS} \
--host-fa ${HOST} 

python3 ${TESTVIRUSCRIPT} \
--int-fa ${INTS} \
--int-info ${LOCS} \
--virus-fa ${VIRUS} 

# check that two runs with the same seed produce the same results
printf "\n\n"

i=1
INTS1="${OUTDIR}integrations_${i}.fa"
LOCS1="${OUTDIR}ints_info_${i}.txt"
EPI1="${OUTDIR}epi_info_${i}.txt"

python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS1} \
 --int_info ${LOCS1} \
 --int_num 10 \
 --seed 12345 \
 --epi_num 5 \
 --epi_info ${EPI1} \
 --p_host_deletion 0.9 \
 --lambda_host_deletion 2 \
 --verbose
 
python3 ${TESTHOSTSCRIPT} \
--int-fa ${INTS1} \
--int-info ${LOCS1} \
--host-fa ${HOST} 

python3 ${TESTVIRUSCRIPT} \
--int-fa ${INTS1} \
--int-info ${LOCS1} \
--virus-fa ${VIRUS} 

diff ${INTS} ${INTS1}
diff ${LOCS} ${LOCS1}
diff ${EPI} ${EPI1}
