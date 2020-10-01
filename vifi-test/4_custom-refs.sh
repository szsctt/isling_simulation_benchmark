#!/bin/bash
set -e

#https://github.com/namphuon/ViFi


module load singularity
cd ..

IMAGE="vifi_latest.sif"
if [ ! -e $IMAGE ]; then
	singularity pull --name $IMAGE docker://namphuon/vifi:latest
fi



#### mimic organisation of data downloaded from google drive
AA_DATA_REPO="out/vifi_refs/data_repo" # contains host sequence
REFERENCE_REPO="out/vifi_refs/data" # contains viral sequences and host/virus combined sequences
mkdir -p ${AA_DATA_REPO}
mkdir -p ${REFERENCE_REPO}

VIRUS="data/references/test_AAV.fa"
VIRUSNAME="rep68"

HOST="data/references/test_human.fa"
HOSTNAME="test_human"

# make bwa indices of virus + host
if [ ! -e $REFERENCE_REPO/${VIRUSNAME}/${HOSTNAME}_${VIRUSNAME}.fas ]; then

	mkdir -p ${REFERENCE_REPO}/${VIRUSNAME}
	mkdir -p ${AA_DATA_REPO}/${HOSTNAME}
	
	ln -s $(realpath ${VIRUS}) $(realpath ${REFERENCE_REPO}/${VIRUSNAME}/${VIRUSNAME}.fa)
	ln -s $(realpath ${HOST}) $(realpath ${AA_DATA_REPO}/${HOSTNAME}/${HOSTNAME}.fa)
	
	echo "concatenating virus and host references"
	cat ${AA_DATA_REPO}/${HOSTNAME}/${HOSTNAME}.fa \
	${REFERENCE_REPO}/${VIRUSNAME}/${VIRUSNAME}.fa \
	> $REFERENCE_REPO/${VIRUSNAME}/${HOSTNAME}_${VIRUSNAME}.fas
	
	echo "indexing combined virus and host"
	srun --time 2:00:00 --mem 5 gb \
		singularity exec -B $(realpath $REFERENCE_REPO/${VIRUSNAME}):/home/repo/data $IMAGE \
		bwa index /home/repo/data/${HOSTNAME}_${VIRUSNAME}.fas
		
	echo "extracting names of host chromosomes"
	#this needs to be a file with a single line with the space-delimited list of chromosomes
	srun --time 2:00:00 --mem 5 gb \
		singularity exec -B $(realpath ${AA_DATA_REPO}/${HOSTNAME}):/home/repo/data $IMAGE \
		samtools faidx /home/repo/data/${HOSTNAME}.fa

	awk 'BEGIN { ORS = " " }
				{a[$1]} 
				END {for (i in a) print i}' \
	${AA_DATA_REPO}/${HOSTNAME}/${HOSTNAME}.fa.fai > \
	${REFERENCE_REPO}/${VIRUSNAME}/${HOSTNAME}_chromosome-list.txt
	
	# stuf that's undocumented, but seems to be required
	echo ${HOSTNAME} > ${AA_DATA_REPO}/reference.txt
	
	echo "fa_file 		                ${HOSTNAME}.fa" > ${AA_DATA_REPO}/${HOSTNAME}/file_list.txt
	echo "chrLen_file 		            ${HOSTNAME}.fa.fai" >> ${AA_DATA_REPO}/${HOSTNAME}/file_list.txt
	

	
fi

INPUT_DIR="out/test/test-easier/sim_reads"
OUTPUT_DIR="out/test/test-easier/vifi"
mkdir -p $OUTPUT_DIR
CPUS="1"
READ1="cond0.rep01.fq"
READ2="cond0.rep02.fq"

export SINGULARITYENV_VIFI_DIR="/scratch1/sco305/intvi_simulation-experiments/ViFi"

srun --time 2:00:00 --mem 30gb \
singularity exec \
-B $(realpath $REFERENCE_REPO/${VIRUSNAME}):/home/repo/data/ \
-B $(realpath $INPUT_DIR):/usr/share/ \
-B $(realpath $AA_DATA_REPO):/home/data_repo/ \
-B $(realpath $OUTPUT_DIR):/opt/ \
$IMAGE \
python $SINGULARITYENV_VIFI_DIR/scripts/run_vifi.py \
-c ${CPUS} \
-f /usr/share/${READ1} -r /usr/share/${READ2} \
--reference /home/repo/data/${HOSTNAME}_${VIRUSNAME}.fas \
-v $VIRUSNAME \
-o /opt/ \
-d True \
-C /home/repo/data//${HOSTNAME}_chromosome-list.txt



