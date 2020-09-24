#!/bin/bash
set -e

cd ..
IMAGE1="polyidus_1.sif"
DOCKER1="docker://szsctt/polyidus:1"
IMAGE2="polyidus_2.sif"
DOCKER2="docker://szsctt/polyidus:2"

module load singularity

if [ ! -e $IMAGE1 ]; then
	echo "pulling container $IMAGE1"
	singularity pull --name $IMAGE1 $DOCKER1
fi

if [ ! -e $IMAGE2 ]; then
	echo "pulling container $IMAGE2"
	singularity pull --name $IMAGE2 $DOCKER2
fi

HOST="polyidus/test/host/test_human_bwt2"
VIRUS="polyidus/test/virus/test_virus_bwt2"
OUT="/polyidus/test/easier"
FASTQ1="out/test/test-easier/sim_reads/cond0.rep01.fq"
FASTQ2="out/test/test-easier/sim_reads/cond0.rep02.fq"

echo "first verion of dockerfile"
singularity exec -B$(pwd) $IMAGE1 python /usr/src/app/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT 

echo
echo "second verion of dockerfile"
singularity exec -B$(pwd) $IMAGE2 python /usr/src/app/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT 

#echo 
#echo "running with modulefile..."
#module load polyidus
#time python polyidus/src/polyidus.py \
#$(pwd)/$HOST $(pwd)/$VIRUS \
#--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
#--outdir $(pwd)/$OUT

