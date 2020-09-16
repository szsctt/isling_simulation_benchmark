#!/bin/bash

#https://github.com/namphuon/ViFi

# download example data: https://drive.google.com/open?id=0ByYcg0axX7udUDRxcTdZZkg0X1k

# download HMM models for HPV and HBV: https://drive.google.com/open?id=0Bzp6XgpBhhghSTNMd3RWS2VsVXM

module load singularity

IMAGE="vifi_latest.sif"
#singularity pull --name $IMAGE docker://namphuon/vifi:latest


# data downloaded from google drive lines above
PWD=$(pwd)
AA_DATA_REPO=$PWD/data_repo
REFERENCE_REPO=$PWD/data


# make bwa indices of virus + host
#cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/hpv/hpv.unaligned.fas > $REFERENCE_REPO/hpv/hg19_hpv.fas

#srun --time 12:00:00 --mem 30gb \
#singularity exec -B $REFERENCE_REPO/hpv:/home/repo/data $IMAGE \
#bwa index /home/repo/data/hg19_hpv.fas

#cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/hbv/hbv.unaligned.fas > $REFERENCE_REPO/hbv/hg19_hbv.fas
#srun --time 12:00:00 --mem 30gb \
#singularity exec -B $REFERENCE_REPO/hbv:/home/repo/data $IMAGE \
#bwa index /home/repo/data/hg19_hbv.fas


INPUT_DIR="/scratch1/sco305/intvi_simulation-experiments/polyidus/data/fastqfiles"
READ1="SiHa_R1.fastq.gz"
READ2="SiHa_R2.fastq.gz"
OUTPUT_DIR=$PWD/test_output_2
mkdir -p $OUTPUT_DIR
CPUS="1"

# vifi command from /home/vifi/sh (inside container)
#docker run -e CPUS=$CPUS -v $REFERENCE_REPO:/home/repo/data -v $INPUT_DIR:/home/fastq/ -e READ1=$READ1 -e READ2=$READ2 -v $AA_DATA_REPO:/home/data_repo/ -v $OUTPUT_DIR:/home/output/ vifi:latest sh /home/vifi.sh

# translated to singularity syntax
srun --time 12:00:00 --mem 30gb \
singularity exec \
-B $REFERENCE_REPO:/home/repo/data \
-B $INPUT_DIR:/home/fastq/ \
-B $AA_DATA_REPO:/home/data_repo/ \
-B $OUTPUT_DIR:/home/output/ \
$IMAGE \
python /home/scripts/run_vifi.py -c ${CPUS} -f /home/fastq/${READ1} -r /home/fastq/${READ2} -v hpv -o /home/output/ -d True


