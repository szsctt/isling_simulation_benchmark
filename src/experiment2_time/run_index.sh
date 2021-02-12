#!/bin/bash
set -e

# ./run_sim.sh <sim config yml> <cluster_config> <vifi_data_repo> <local>

CONFIG=$1
CLUSTER=$2
VIFI=$3

if [ -z "$4" ]
	then
		eval "$(conda shell.bash hook)"
		conda activate snakemake
		module load singularity
		SNAKEARGS="--profile slurm --cluster-config ${CLUSTER} --latency-wait 120 --jobs 100"
else
	SNAKEARGS="--resources mem_mb=60000 --cores 15"
fi

WD=$(pwd)

cd ../../intvi_pipeline

# make config files from simulation config
CONFIGFOLDER="$(dirname ${CONFIG})/auto"
ISLCONFIG="${CONFIGFOLDER}/isling_index.yml"
OTHCONFIG="${CONFIGFOLDER}/other_tools_index.yml"
mkdir -p $CONFIGFOLDER
echo "creating config files for indexing...\n"
python3 ${WD}/make_index_configs.py ${CONFIG} ${VIFI} ${ISLCONFIG} ${OTHCONFIG}

# run snakemake until index rule
echo 
echo "indexing references for isling"
snakemake \
 --configfile "${ISLCONFIG}" \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --rerun-incomplete \
 --until index \
 --quiet ${SNAKEARGS}

# also make config file for other tools
cd ../intvi_other-tools
echo
echo "indexing references for polyidus"
snakemake \
 --configfile "${OTHCONFIG}" \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --rerun-incomplete \
 --until bwt2_index ${SNAKEARGS}

echo
echo "indexing references for vifi"
snakemake \
 --configfile "${OTHCONFIG}" \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --rerun-incomplete \
 --until host_virus_index ${SNAKEARGS}
 
echo
echo "indexing references for seeksv"
snakemake \
 --configfile "${OTHCONFIG}" \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --rerun-incomplete \
 --until host_virus_index_seeksv ${SNAKEARGS}

echo
echo "indexing references for vseq-toolkit"
snakemake \
 --configfile "${OTHCONFIG}" \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --rerun-incomplete \
 --until bwa_index ${SNAKEARGS}

echo
echo "done indexing references"
