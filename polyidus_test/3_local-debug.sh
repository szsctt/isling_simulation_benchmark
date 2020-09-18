#!/bin/bash
set -e

# index host and viral references
cd ~/Documents/Projects/viInt/experiments/expt6_simulations
mkdir -p polyidus_test/refs/test_human
mkdir -p polyidus_test/refs/test_AAV
mkdir -p polyidus_test/test-easier

docker run --rm -it -v$(pwd):/usr/data szsctt/polyidus:latest \
bowtie2-build /usr/data/data/references/test_human.fa /usr/data/polyidus_test/refs/test_human/test_human

docker run --rm -it -v$(pwd):/usr/data szsctt/polyidus:latest \
bowtie2-build /usr/data/data/references/test_AAV.fa /usr/data/polyidus_test/refs/test_AAV/test_AAV

docker run --rm -it -v$(pwd):/usr/data szsctt/polyidus:latest python /usr/data/polyidus/src/polyidus.py \
/usr/data/polyidus_test/refs/test_human/test_human \
/usr/data/polyidus_test/refs/test_AAV/test_AAV \
--fastq /usr/data/out/test/test-easier/sim_reads/cond0.rep01.fq /usr/data/out/test/test-easier/sim_reads/cond0.rep02.fq \
--outdir /usr/data/polyidus_test/test-easier \
--skip-alignment
