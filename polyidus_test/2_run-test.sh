#!/bin/bash
set -e

eval "$(conda shell.bash hook)"


# simulate data
cd ../src/sim_tests
bash test-sim-analysis_combined.sh

# bowtie2 indicies of host and virus
conda activate polyidus
cd ../../polyidus
mkdir -p test
cd test
mkdir -p host
cd host
HOSTREF="../../../data/references/test_human.fa"
srun --time 2:00:00 --mem 50gb \
bowtie2-build $HOSTREF test_human_bwt2

cd ..
mkdir -p virus
cd virus
VIRUSREF="../../../data/references/test_AAV.fa"
srun --time 2:00:00 --mem 50gb \
bowtie2-build $VIRUSREF test_virus_bwt2

# analyse data
cd ../../src
OUTPUT="../test/easier"
FQ1="../../out/test/test-easier/sim_reads/cond0.rep01.fq"
FQ2="../../out/test/test-easier/sim_reads/cond0.rep02.fq"
mkdir -p $OUTPUT
srun --time 2:00:00 --mem 50gb \
python3 polyidus.py \
../test/host/test_human_bwt2 \
../test/virus/test_virus_bwt2 \
--fastq $FQ1 $FQ2 \
--outdir $OUTPUT

OUTPUT="../test/harder"
FQ1="../../out/test/test-harder/sim_reads/cond0.rep01.fq"
FQ2="../../out/test/test-harder/sim_reads/cond0.rep02.fq"
mkdir -p $OUTPUT
srun --time 2:00:00 --mem 50gb \
python3 polyidus.py \
../test/host/test_human_bwt2 \
../test/virus/test_virus_bwt2 \
--fastq $FQ1 $FQ2 \
--outdir $OUTPUT
