#!/bin/bash
set -e

#https://github.com/namphuon/ViFi

# download example data: https://drive.google.com/open?id=0ByYcg0axX7udUDRxcTdZZkg0X1k

# download HMM models for HPV and HBV: https://drive.google.com/open?id=0Bzp6XgpBhhghSTNMd3RWS2VsVXM

module load singularity

PWD=$(pwd)
IMAGE=$PWD/"vifi_latest.sif"
if [ ! $IMAGE ]; then
	singularity pull --name $IMAGE docker://namphuon/vifi:latest
fi

# data downloaded from google drive lines above
AA_DATA_REPO=$PWD/data_repo
REFERENCE_REPO=$PWD/data


# make bwa indices of virus + host
if [ ! $REFERENCE_REPO/hg19_hbv.fas ]; then
	cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/hpv/hpv.unaligned.fas > $REFERENCE_REPO/hpv/hg19_hpv.fas

	srun --time 12:00:00 --mem 30gb \
	singularity exec -B $REFERENCE_REPO/hpv:/home/repo/data $IMAGE \
	bwa index /home/repo/data/hg19_hpv.fas
fi

INPUT_DIR="/scratch1/sco305/intvi_simulation-experiments/polyidus/data/fastqfiles"
READ1="SiHa_R1.fastq.gz"
READ2="SiHa_R2.fastq.gz"
OUTPUT_DIR=$PWD/test_output_0
mkdir -p $OUTPUT_DIR
CPUS="1"
export SINGULARITYENV_VIFI_DIR="/scratch1/sco305/intvi_simulation-experiments/ViFi"

# vifi command from /home/vifi/sh (inside container)
#docker run -e CPUS=$CPUS -v $REFERENCE_REPO:/home/repo/data -v $INPUT_DIR:/home/fastq/ -e READ1=$READ1 -e READ2=$READ2 -v $AA_DATA_REPO:/home/data_repo/ -v $OUTPUT_DIR:/home/output/ vifi:latest sh /home/vifi.sh

# translated to singularity syntax
srun --time 12:00:00 --mem 30gb \
singularity exec \
-B $REFERENCE_REPO:/home/repo/data \
-B $INPUT_DIR:/usr/share/ \
-B $AA_DATA_REPO:/home/data_repo/ \
-B $OUTPUT_DIR:/home/output/ \
$IMAGE \
python /home/scripts/run_vifi.py -c ${CPUS} -f /usr/share/${READ1} -r /usr/share/${READ2} -v hpv -o /home/output/

