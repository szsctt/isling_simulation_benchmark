#!/bin/bash
set -e

# usage: ./run_seeksv.sh <fq1> <fq2> <bwa_prefix> <threads> <out_dir>

FQ1=$1
FQ2=$2
PRF=$3
THR=$4
OUT=$5

# make output directory
mkdir -p ${OUT}

ALN="${OUT}/aln1.bam"

CLPRF="${OUT}/getclip"
CLP="${CLPRF}.clip.gz"
CFQ="${CLPRF}.clip.fq.gz"
UM1="${CLPRF}.unmapped_1.fq.gz"
UM2="${CLPRF}.unmapped_2.fq.gz"

ALN2="${OUT}/aln2.bam"

SV="${OUT}/output.sv.txt"
UM="${OUT}/output.unmapped.clip.fq.gz"

# align reads
bwa mem -t ${THR} ${PRF} ${FQ1} ${FQ2} | samtools sort -o ${ALN} -
samtools index ${ALN}

# seeksv get clip
/var/work/seeksv/seeksv getclip -o ${CLPRF} ${ALN}

# align clipped reads
bwa mem -t ${THR} ${PRF} ${CFQ} | samtools view  -Sb -o ${ALN2} -

# get output
/var/work/seeksv/seeksv getsv ${ALN2} ${ALN} ${CFQ} ${SV} ${UM}

