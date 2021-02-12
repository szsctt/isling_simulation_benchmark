#!/bin/bash
set -e

cd ../../intvi_pipeline/

CONFIGDIR="../config/experiment2_time/split"

OUTDIR="../out/experiment2_time/split/"
mkdir -p ${OUTDIR}

eval "$(conda shell.bash hook)"
conda activate snakemake
module load singularity
module load parallel

SRUN="srun --time 12:00:00 --mem 60gb -c10"
SINGULARITY="singularity run -B $(realpath ..) isling.sif"
TIME="/usr/bin/time"
SNAKEMAKE="snakemake --cores 10 --resources mem_mb=60000 --forceall"



ls ${CONFIGDIR}/* | parallel basename | parallel --delay 0.1 "date > ${OUTDIR}/$(basename {}).log; ${SRUN} ${TIME} ${SINGULARITY} ${SNAKEMAKE} --configfile ${CONFIGDIR}/{} >> ${OUTDIR}/$(basename {}).log 2>&1; echo ${?} >> ${OUTDIR}/$(basename {}).log; date >> ${OUTDIR}/$(basename {}).log"


