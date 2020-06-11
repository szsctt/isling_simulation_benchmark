#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate simvi

SCRIPT="../../intvi_simulation/insert_virus.py"
HOST="../../data/references/test_host.fa"
VIRUS="../../data/references/test_virus.fa"

OUTDIR="../../out/test/"
mkdir -p ${OUTDIR}

echo "*************************"
echo "**** starting test 1 ****"
echo "*************************"

i=0
INTS="${OUTDIR}integrations_${i}.fa"
LOCS="${OUTDIR}locs_${i}.txt"
HOSTINTS="${OUTDIR}host_locs_${i}.txt"

python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS} \
 --locs ${LOCS} \
 --locs_host ${HOSTINTS} \
 --int_num 1 \
 --seed 1 \
 &
 
((i=i+1))
INTS="${OUTDIR}integrations_${i}.fa"
LOCS="${OUTDIR}locs_${i}.txt"
HOSTINTS="${OUTDIR}host_locs_${i}.txt"

python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS} \
 --locs ${LOCS} \
 --locs_host ${HOSTINTS} \
 --int_num 1 \
 --seed 1 \
 &
 
wait
 
 # make sure that files are the same with the same seed
SEED_TEST="${OUTDIR}seed_test.txt"
diff "${OUTDIR}integrations_0.fa" "${OUTDIR}integrations_1.fa" | wc -l > ${SEED_TEST}
diff "${OUTDIR}locs_0.txt" "${OUTDIR}locs_1.txt" | wc -l >> ${SEED_TEST}
diff "${OUTDIR}host_locs_0.txt" "${OUTDIR}host_locs_1.txt" | wc -l >> ${SEED_TEST}
awk '($0 !~ /0/){exit 1}' ${SEED_TEST}
rm "${OUTDIR}integrations_0.fa" "${OUTDIR}integrations_1.fa" "${OUTDIR}locs_0.txt" "${OUTDIR}locs_1.txt"  "${OUTDIR}host_locs_0.txt" "${OUTDIR}host_locs_1.txt" ${SEED_TEST}
echo "*************************"
echo "***** passed test 1 *****"
echo "*************************"


echo "*************************"
echo "**** starting test 2 ****"
echo "*************************"

i=0
INTS="${OUTDIR}integrations_${i}.fa"
LOCS="${OUTDIR}locs_${i}.txt"
HOSTINTS="${OUTDIR}host_locs_${i}.txt"

python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS} \
 --locs ${LOCS} \
 --locs_host ${HOSTINTS} \
 --int_num 1 \
 --set_len 10 \
 --seed 1 \
 &
 
((i=i+1))
INTS="${OUTDIR}integrations_${i}.fa"
LOCS="${OUTDIR}locs_${i}.txt"
HOSTINTS="${OUTDIR}host_locs_${i}.txt"

python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS} \
 --locs ${LOCS} \
 --locs_host ${HOSTINTS} \
 --int_num 1 \
 --set_len 10 \
 --seed 1 \
 &
 
wait
 
 # make sure that files are the same with the same seed
SEED_TEST="${OUTDIR}seed_test.txt"
diff "${OUTDIR}integrations_0.fa" "${OUTDIR}integrations_1.fa" | wc -l > ${SEED_TEST}
diff "${OUTDIR}locs_0.txt" "${OUTDIR}locs_1.txt" | wc -l >> ${SEED_TEST}
diff "${OUTDIR}host_locs_0.txt" "${OUTDIR}host_locs_1.txt" | wc -l >> ${SEED_TEST}
awk '($0 !~ /0/){exit 1}' ${SEED_TEST}
rm "${OUTDIR}integrations_0.fa" "${OUTDIR}integrations_1.fa" "${OUTDIR}locs_0.txt" "${OUTDIR}locs_1.txt"  "${OUTDIR}host_locs_0.txt" "${OUTDIR}host_locs_1.txt" ${SEED_TEST}
echo "*************************"
echo "***** passed test 2 *****"
echo "*************************"


echo "*************************"
echo "**** starting test 3 ****"
echo "*************************"

i=0
INTS="${OUTDIR}integrations_${i}.fa"
LOCS="${OUTDIR}locs_${i}.txt"
HOSTINTS="${OUTDIR}host_locs_${i}.txt"

python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS} \
 --locs ${LOCS} \
 --locs_host ${HOSTINTS} \
 --sep 5 \
 --int_portion "whole" \
 --int_deletion "none" \
 --int_rearrange "none" \
 --set_junc "clean" \
 --int_num 3 \
 --seed 1
 
((i=i+1))
INTS="${OUTDIR}integrations_${i}.fa"
LOCS="${OUTDIR}locs_${i}.txt"
HOSTINTS="${OUTDIR}host_locs_${i}.txt"

python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS} \
 --locs ${LOCS} \
 --locs_host ${HOSTINTS} \
 --sep 5 \
 --int_portion "whole" \
 --int_deletion "none" \
 --int_rearrange "none" \
 --set_junc "clean" \
 --int_num 3 \
 --seed 1
 
wait
 
 # make sure that files are the same with the same seed
SEED_TEST="${OUTDIR}seed_test.txt"
diff "${OUTDIR}integrations_0.fa" "${OUTDIR}integrations_1.fa" | wc -l > ${SEED_TEST}
diff "${OUTDIR}locs_0.txt" "${OUTDIR}locs_1.txt" | wc -l >> ${SEED_TEST}
diff "${OUTDIR}host_locs_0.txt" "${OUTDIR}host_locs_1.txt" | wc -l >> ${SEED_TEST}
awk '($0 !~ /0/){exit 1}' ${SEED_TEST}
rm "${OUTDIR}integrations_0.fa" "${OUTDIR}integrations_1.fa" "${OUTDIR}locs_0.txt" "${OUTDIR}locs_1.txt"  "${OUTDIR}host_locs_0.txt" "${OUTDIR}host_locs_1.txt" ${SEED_TEST}
echo "*************************"
echo "***** passed test 3 *****"
echo "*************************"

