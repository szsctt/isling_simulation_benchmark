#!/bin/bash
set -e

# ./run_tools.sh <sim config yml> <isling docker string> <seeksv docker string> <polyidus docker string> 
#   <vifi docker string> <vifi data repo> <vseq-toolkit docker string> <cores> <local>

CONFIG=$1
ISDOCK=$2
SEDOCK=$3
PODOCK=$4
VIDOCK=$5
VIDATA=$6
VSDOCK=$7
CORES=$8

if [ -z "$9" ]
	then
		eval "$(conda shell.bash hook)"
		conda activate snakemake
		module load singularity
		ARGS="--parallel"
fi

ISCONT="isling.sif"
SECONT="seeksv.sif"
POCONT="polyidus.sif"
VICONT="vifi.sif"
VSCONT="vseq.sif"
WD=$(pwd)

cd ../../intvi_pipeline

# pull docker containers
if [ ! -e ${ISCONT} ]; then
	singularity pull --name ${ISCONT} ${ISDOCK}
fi
if [ ! -e ${SECONT} ]; then
	singularity pull --name ${SECONT} ${SEDOCK}
fi
if [ ! -e ${POCONT} ]; then
	singularity pull --name ${POCONT} ${PODOCK}
fi
if [ ! -e ${VICONT} ]; then
	singularity pull --name ${VICONT} ${VIDOCK}
fi
if [ ! -e ${VSCONT} ]; then
	singularity pull --name ${VSCONT} ${VSDOCK}
fi

python3 ${WD}/run_tools.py \
 --sim-config ${CONFIG} \
 --isling-sif ${ISCONT} \
 --isling-config-script "${WD}/make_isling_config.py" \
 --seeksv-sif ${SECONT} \
 --seeksv-script "${WD}/run_seeksv.sh" \
 --polyidus-sif ${POCONT} \
 --vifi-data-repo ${VIDATA} \
 --vifi-sif ${VICONT} \
 --vseq-toolkit-sif ${VSCONT} \
 --replicates 3 \
 --cores ${CORES} ${ARGS}
