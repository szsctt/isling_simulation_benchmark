#!/bin/bash

cd ..
IMAGE="polyidus.sif"
DOCKER="docker://szsctt/polyidus:3"

module load singularity

if [ ! -e $IMAGE ]; then
	echo "pulling container"
	singularity pull --name $IMAGE $DOCKER 
fi


# make bowtie2 indicies
HOSTFA="data/references/GRCh38.fa"
HOST="polyidus/data/hg38_whole/hg38"

mkdir -p $(dirname $HOST)
if [ ! -e $HOST.1.bt2 ]; then
	srun --time 2:00:00 --mem 50gb \
	singularity exec $IMAGE \
	bowtie2-build $HOSTFA $HOST
fi


VIRUSFA="data/references/NC_001401.2.fa"
VIRUS="polyidus/data/AAV2/AAV2"

mkdir -p $(dirname $VIRUS)
if [ ! -e $VIRUS.1.bt2 ]; then
	srun --time 2:00:00 --mem 50gb \
	singularity exec $IMAGE \
	bowtie2-build $VIRUSFA $VIRUS
fi

# use previously simulated data
FASTQ1="out/experiment0_prelim/AAV2-easier/sim_reads/cond0.rep01.fq"
FASTQ2="out/experiment0_prelim/AAV2-easier/sim_reads/cond0.rep02.fq"

OUT="polyidus/polyidusOutput_singularity_larger"
mkdir -p $OUT
rm -rf $OUT/*


echo "running with singularity image..."
srun --time 2:00:00 --mem 100gb \
time singularity exec -B$(pwd) $IMAGE python3 $(pwd)/polyidus/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT > $OUT/log.log
sort -k1,1 -k2,2n $OUT/results/exactHpvIntegrations.tsv | cut -f1-9  > $OUT/results/exactHpvIntegrations.sorted.tsv

module unload singularity

OUT="polyidus/polyidusOutput_module_larger"
mkdir -p $OUT
rm -rf $OUT/*

echo 
echo "running with modulefile..."
module load polyidus
srun --time 2:00:00 --mem 100gb \
time python3 $(pwd)/polyidus/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT > $OUT/log.log
sort -k1,1 -k2,2n $OUT/results/exactHpvIntegrations.tsv | cut -f1-9 > $OUT/results/exactHpvIntegrations.sorted.tsv

