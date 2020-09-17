#!/bin/bash
set -e

eval "$(conda shell.bash hook)"
conda activate polyidus

cd ../polyidus
mkdir -p data/polyidusOutput
cd src

srun --time 2:00:00 --mem 50gb \
python3 polyidus.py \
../data/hg38/hg38_bwt2_index \
../data/hpv16/hpv16_bowt_ind \
--fastq ../data/fastqfiles/SiHa_R1.fastq.gz ../data/fastqfiles/SiHa_R2.fastq.gz \
--outdir ../data/polyidusOutput
