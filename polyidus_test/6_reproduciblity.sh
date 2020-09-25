#!/bin/bash

cd ..
IMAGE="polyidus.sif"
DOCKER="docker://szsctt/polyidus:3"

module load singularity

if [ ! -e $IMAGE ]; then
	echo "pulling container"
	singularity pull --name $IMAGE $DOCKER 
fi

HOST="polyidus/data/hg38/hg38_bwt2_index"
VIRUS="polyidus/data/hpv16/hpv16_bowt_ind"

FASTQ1="polyidus/data/fastqfiles/SiHa_R1.fastq.gz"
FASTQ2="polyidus/data/fastqfiles/SiHa_R2.fastq.gz"

OUT="polyidus/polyidusOutput_singularity_1"
mkdir -p $OUT
rm -rf $OUT/*


echo "running with singularity image..."
singularity exec -B$(pwd) $IMAGE python3 $(pwd)/polyidus/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT > $OUT/log.log
sort -k1,1 -k2,2n $OUT/results/exactHpvIntegrations.tsv | rev | cut -f3- | rev > $OUT/results/exactHpvIntegrations.sorted.tsv


OUT="polyidus/polyidusOutput_singularity_2"
mkdir -p $OUT
rm -rf $OUT/*

echo
echo "running with singularity image..."
singularity exec -B$(pwd) $IMAGE python3 $(pwd)/polyidus/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT > $OUT/log.log
sort -k1,1 -k2,2n $OUT/results/exactHpvIntegrations.tsv | rev | cut -f3- | rev > $OUT/results/exactHpvIntegrations.sorted.tsv


module unload singularity

echo 
echo "comparing the two singularity runs"
diff polyidus/polyidusOutput_singularity_1/results/exactHpvIntegrations.sorted.tsv \
polyidus/polyidusOutput_singularity_2/results/exactHpvIntegrations.sorted.tsv

eval "$(conda shell.bash hook)"
conda activate polyidus

OUT="polyidus/polyidusOutput_conda_1"
mkdir -p $OUT
rm -rf $OUT/*

echo 
echo "running with conda..."
module load polyidus
python3 $(pwd)/polyidus/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT > $OUT/log.log
sort -k1,1 -k2,2n $OUT/results/exactHpvIntegrations.tsv | rev | cut -f3- | rev > $OUT/results/exactHpvIntegrations.sorted.tsv

OUT="polyidus/polyidusOutput_conda_2"
mkdir -p $OUT
rm -rf $OUT/*

echo 
echo "running with conda..."
module load polyidus
python3 $(pwd)/polyidus/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT > $OUT/log.log
sort -k1,1 -k2,2n $OUT/results/exactHpvIntegrations.tsv | rev | cut -f3- | rev > $OUT/results/exactHpvIntegrations.sorted.tsv

conda deactivate polyidus

echo 
echo "comparing the two conda runs"
diff polyidus/polyidusOutput_conda_1/results/exactHpvIntegrations.sorted.tsv \
polyidus/polyidusOutput_conda_2/results/exactHpvIntegrations.sorted.tsv

OUT="polyidus/polyidusOutput_module_1"
mkdir -p $OUT
rm -rf $OUT/*

echo 
echo "running with modulefile..."
module load polyidus
python3 $(pwd)/polyidus/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT > $OUT/log.log
sort -k1,1 -k2,2n $OUT/results/exactHpvIntegrations.tsv | rev | cut -f3- | rev > $OUT/results/exactHpvIntegrations.sorted.tsv

OUT="polyidus/polyidusOutput_module_2"
mkdir -p $OUT
rm -rf $OUT/*

echo 
echo "running with modulefile..."
python3 $(pwd)/polyidus/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT > $OUT/log.log
sort -k1,1 -k2,2n $OUT/results/exactHpvIntegrations.tsv | rev | cut -f3- | rev > $OUT/results/exactHpvIntegrations.sorted.tsv

module unload polyidus

echo 
echo "comparing the two module runs"
diff polyidus/polyidusOutput_module_1/results/exactHpvIntegrations.sorted.tsv \
polyidus/polyidusOutput_module_2/results/exactHpvIntegrations.sorted.tsv


echo 
echo "comparing module 1 with conda 1"
diff polyidus/polyidusOutput_module_1/results/exactHpvIntegrations.sorted.tsv \
polyidus/polyidusOutput_conda_1/results/exactHpvIntegrations.sorted.tsv

echo 
echo "comparing module 1 with singularity 1"
diff polyidus/polyidusOutput_module_1/results/exactHpvIntegrations.sorted.tsv \
polyidus/polyidusOutput_singularity_1/results/exactHpvIntegrations.sorted.tsv

echo 
echo "comparing conda 1 with singularity 1"
diff polyidus/polyidusOutput_conda_1/results/exactHpvIntegrations.sorted.tsv \
polyidus/polyidusOutput_singularity_1/results/exactHpvIntegrations.sorted.tsv


