#!/bin/bash
set -e

CONT="szsctt/polyidus:3"

cd ..
HOST="polyidus/data/hg38/hg38_bwt2_index"
VIRUS="polyidus/data/hpv16/hpv16_bowt_ind"

FASTQ1="polyidus/data/fastqfiles/SiHa_R1.fastq.gz"
FASTQ2="polyidus/data/fastqfiles/SiHa_R2.fastq.gz"

OUT="polyidus/polyidusOutput_docker_1"
mkdir -p $OUT
rm -rf $OUT/*


docker run -v$(pwd):/home --rm $CONT \
/opt/conda/bin/python /usr/src/app/src/polyidus.py \
/home/$HOST /home/$VIRUS \
--fastq /home/$FASTQ1 /home/$FASTQ2 \
--outdir /home/$OUT
sort -k1,1 -k2,2n $OUT/results/exactHpvIntegrations.tsv | cut -f1-9 > $OUT/results/exactHpvIntegrations.sorted.tsv


OUT="polyidus/polyidusOutput_docker_2"
mkdir -p $OUT
rm -rf $OUT/*

docker run -v$(pwd):/home --rm $CONT \
/opt/conda/bin/python /usr/src/app/src/polyidus.py \
/home/$HOST /home/$VIRUS \
--fastq /home/$FASTQ1 /home/$FASTQ2 \
--outdir /home/$OUT
sort -k1,1 -k2,2n $OUT/results/exactHpvIntegrations.tsv | cut -f1-9 > $OUT/results/exactHpvIntegrations.sorted.tsv


diff polyidus/polyidusOutput_docker_1/results/exactHpvIntegrations.sorted.tsv \
polyidus/polyidusOutput_docker_2/results/exactHpvIntegrations.sorted.tsv 

