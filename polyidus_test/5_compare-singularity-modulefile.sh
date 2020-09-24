#!/bin/bash

cd ..
IMAGE="polyidus.sif"
DOCKER="docker://szsctt/polyidus:2"

module load singularity

if [ ! -e $IMAGE ]; then
	echo "pulling container"
	singularity pull --name $IMAGE $DOCKER 
fi

HOST="polyidus/data/hg38/hg38_bwt2_index"
VIRUS="polyidus/data/hpv16/hpv16_bowt_ind"
OUT="polyidus/polyidusOutput_singularity"
FASTQ1="polyidus/data/fastqfiles/SiHa_R1.fastq.gz"
FASTQ2="polyidus/data/fastqfiles/SiHa_R2.fastq.gz"

mkdir -p $OUT

echo "running with singularity image..."
time singularity exec -B$(pwd) $IMAGE python $(pwd)/polyidus/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT > polyidus/singularity3.log


OUT="polyidus/polyidusOutput_module"
mkdir -p $OUT
echo 
echo "running with modulefile..."
module load polyidus
time python $(pwd)/polyidus/src/polyidus.py \
$(pwd)/$HOST $(pwd)/$VIRUS \
--fastq $(pwd)/$FASTQ1 $(pwd)/$FASTQ2 \
--outdir $(pwd)/$OUT > polyidus/module3.log


awk 'NF{NF-=1};1' polyidus/polyidusOutput_singularity/results/exactHpvIntegrations.tsv | sort -k1,1 -k2,2n > polyidus/polyidusOutput_singularity/results/exactHpvIntegrations.sorted.tsv 

awk 'NF{NF-=1};1' polyidus/polyidusOutput_module/results/exactHpvIntegrations.tsv | sort -k1,1 -k2,2n > polyidus/polyidusOutput_module/results/exactHpvIntegrations.sorted.tsv 

echo "checking for differences between outputs"
diff polyidus/polyidusOutput_singularity/results/exactHpvIntegrations.sorted.tsv polyidus/polyidusOutput_module/results/exactHpvIntegrations.sorted.tsv

# running with singularity image...
## run 1
### real	0m33.845s
### user	0m26.765s
### sys	0m6.447s

## run 2
### real	0m31.471s
### user	0m29.655s
### sys	0m1.531s



# running with modulefile...
## run 1
### real	0m31.638s
### user	0m24.707s
### sys	0m6.249s

## run 2
### real	0m27.746s
### user	0m25.363s
### sys	0m1.393s





