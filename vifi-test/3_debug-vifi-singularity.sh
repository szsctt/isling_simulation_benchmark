#!/bin/bash

#https://github.com/namphuon/ViFi


module load singularity

PWD=$(pwd)
IMAGE=$PWD/"vifi_latest.sif"
if [ ! $IMAGE ]; then
	singularity pull --name $IMAGE docker://namphuon/vifi:latest
fi

# data downloaded from google drive lines above
AA_DATA_REPO=$PWD/data_repo
REFERENCE_REPO=$PWD/data


#mkdir -p $REFERENCE_REPO/otc
#ln -s $(realpath ../data/references/OTC-vec_rAAV-genome-only.fa) $(realpath $REFERENCE_REPO/otc)


# make bwa indices of virus + host
#cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/otc/OTC-vec_rAAV-genome-only.fa > $REFERENCE_REPO/otc/hg19_otc.fas

#srun --time 12:00:00 --mem 30gb \
#singularity exec -B $REFERENCE_REPO/otc:/home/repo/data $IMAGE \
#bwa index /home/repo/data/hg19_otc.fas


INPUT_DIR="/scratch1/sco305/intvi_simulation-experiments/out/test/test-easier/sim_reads"
READ1="cond0.rep01.fq"
READ2="cond0.rep02.fq"
VIFI_DIR="/scratch1/sco305/intvi_simulation-experiments/ViFi"
OUTPUT_DIR=$PWD/test_output_3
mkdir -p $OUTPUT_DIR
CPUS="1"

cd ..
srun --time 12:00:00 --mem 30gb \
singularity exec \
-B $REFERENCE_REPO:/home/repo/data \
-B $INPUT_DIR:/usr/share/ \
-B $AA_DATA_REPO:/home/data_repo/ \
-B $OUTPUT_DIR:/home/output/ \
$IMAGE \
/bin/bash

# run inside container
VIFI_DIR="/scratch1/sco305/intvi_simulation-experiments/ViFi"
CPUS="1"
READ1="cond0.rep01.fq"
READ2="cond0.rep02.fq"
#python /home/scripts/run_vifi.py -c ${CPUS} -f /usr/share/${READ1} -r /usr/share/${READ2} -v otc -o /home/output/ -d True




