#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate simvi

SCRIPT="../../intvi_simulation/insert_virus_simple.py"
HOST="../../data/references/test_host.fa"
VIRUS="../../data/references/test_virus.fa"

OUTDIR="../../out/test/"
mkdir -p ${OUTDIR}

echo "*************************"
echo "**** starting test 0 ****"
echo "*************************"

i=0
INTS="${OUTDIR}integrations_${i}.fa"
LOCS="${OUTDIR}info_${i}.txt"

python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS} \
 --info ${LOCS} \
 --int_num 10 \
 --seed 1 \
 --verbose

 
((i=i+1))
INTS="${OUTDIR}integrations_${i}.fa"
LOCS="${OUTDIR}info_${i}.txt"

python3 ${SCRIPT} \
 --host ${HOST} \
 --virus ${VIRUS} \
 --ints ${INTS} \
 --info ${LOCS} \
 --int_num 10 \
 --seed 1 \
 --verbose

 
wait
 
 # make sure that files are the same with the same seed
SEED_TEST="${OUTDIR}seed_test.txt"
diff "${OUTDIR}integrations_0.fa" "${OUTDIR}integrations_1.fa" | wc -l > ${SEED_TEST}
diff "${OUTDIR}info_0.txt" "${OUTDIR}info_1.txt" | wc -l >> ${SEED_TEST}
awk '($0 !~ /0/){exit 1}' ${SEED_TEST}
#rm "${OUTDIR}integrations_0.fa" "${OUTDIR}integrations_1.fa" "${OUTDIR}locs_0.txt" "${OUTDIR}locs_1.txt"  "${OUTDIR}host_locs_0.txt" "${OUTDIR}host_locs_1.txt" ${SEED_TEST}
echo "*************************"
echo "***** passed test 0 *****"
echo "*************************"
