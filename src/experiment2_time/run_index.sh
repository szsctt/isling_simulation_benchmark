#!/bin/bash
set -e

# ./run_sim.sh <sim config yml> <cluster_config> <vifi_data_repo>

CONFIG=$1
CLUSTER=$2
VIFI=$3

WD=$(pwd)

eval "$(conda shell.bash hook)"
conda activate snakemake

module load singularity

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
 --jobs 5 \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --profile slurm \
 --rerun-incomplete \
 --latency-wait 120 \
 --until index \
 --cluster-config ${CLUSTER} \
 --quiet

# also make config file for other tools
cd ../intvi_other-tools
echo
echo "indexing references for polyidus"
snakemake \
 --configfile "${OTHCONFIG}" \
 --jobs 5 \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --profile slurm \
 --rerun-incomplete \
 --latency-wait 120 \
 --until bwt2_index

echo
echo "indexing references for vifi"
snakemake \
 --configfile "${OTHCONFIG}" \
 --jobs 5 \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --profile slurm \
 --rerun-incomplete \
 --latency-wait 120 \
 --until host_virus_index
 
echo
echo "indexing references for seeksv"
snakemake \
 --configfile "${OTHCONFIG}" \
 --jobs 5 \
 --use-singularity \
 --singularity-args "-B $(realpath ../)" \
 --profile slurm \
 --rerun-incomplete \
 --latency-wait 120 \
 --until host_virus_index_seeksv

echo
echo "done indexing references"
