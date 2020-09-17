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

mkdir -p $REFERENCE_REPO/rep68
# make bwa indices of virus + host
if [ ! $REFERENCE_REPO/rep68/rep68.fa ]; then

	ln -s $(realpath ../data/references/test_AAV.fa) $(realpath $REFERENCE_REPO/rep68/rep68.fa)
	
	cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/rep68/rep68.fa > $REFERENCE_REPO/rep68/hg19_rep68.fas
	
	srun --time 12:00:00 --mem 30gb \
		singularity exec -B $REFERENCE_REPO/rep68:/home/repo/data $IMAGE \
		bwa index /home/repo/data/hg19_rep68.fas
fi
cd ../..

INPUT_DIR="/scratch1/sco305/intvi_simulation-experiments/out/test/test-easier/sim_reads"
OUTPUT_DIR=$PWD/test_output_2
mkdir -p $OUTPUT_DIR
export SINGULARITYENV_VIFI_DIR="/scratch1/sco305/intvi_simulation-experiments/ViFi"
CPUS="1"
READ1="cond0.rep01.fq"
READ2="cond0.rep02.fq"


srun --time 2:00:00 --mem 30gb \
singularity exec \
-B $REFERENCE_REPO:/home/repo/data \
-B $INPUT_DIR:/usr/share/ \
-B $AA_DATA_REPO:/home/data_repo/ \
-B $OUTPUT_DIR:/home/output/ \
$IMAGE \
python $SINGULARITYENV_VIFI_DIR/scripts/run_vifi.py -c ${CPUS} -f /usr/share/${READ1} -r /usr/share/${READ2} -v rep68 -o /home/output/ -d True





