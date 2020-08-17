#!/bin/bash
set -e

# align simulated reads to host
# use all reads, also reads that were indicated to be integrations
# this is a nice visual of what the pipeline is doing

module load singularity

BWAPULL="docker://szsctt/bwa:1"
BWAIMG="bwa.sif"
BWA="singularity exec ${BWAIMG} bwa"
SAMTOOLS="singularity exec ${BWAIMG} samtools"

SIMVIPULL="docker://szsctt/simvi:2"
SIMVIIMG="simvi.sif"
PYTHON="singularity exec ${SIMVIIMG} python3"

CWD=$(pwd)

# data
WORK="../../out/test-combined"
SAMPLE="cond0.rep0"
READS="sim_reads"
R1="${SAMPLE}1.fq"
R2="${SAMPLE}2.fq"

# relative to $WORK
INTS="../test-combined_analysis1_human_AAV/ints/${SAMPLE}.human.AAV.integrations.txt"

# host bwa index
HOST="../../out/references/human/human"

# filtering script
SCRIPT="../../intvi_simulation/scripts/filter_host_align.py"

# outputs
OUT="demo-align"
ALL="${SAMPLE}.all.bam"
CHIM="${SAMPLE}.chimeric.bam"
DISC="${SAMPLE}.discordant.bam"
BOTH="${SAMPLE}.ints.bam"

# make output folder
cd ${WORK}
mkdir -p ${OUT}
pwd

# make necessary symlinks
ln -sf $(realpath ${CWD}/${HOST}*) .
ln -sf $(realpath $CWD/${SCRIPT}) .
ln -sf $(realpath $INTS) .

# pull singularity images
if [ ! -f ${BWAIMG} ] ; then
	echo pulling bwa singularity image
	singularity pull --name ${BWAIMG} ${BWAPULL}
fi

if [ ! -f ${SIMVIIMG} ] ; then
	echo pulling simvi singularity image
	singularity pull --name ${SIMVIIMG} ${SIMVIPULL}
fi

# align reads to host
echo aligning all reads to host

${BWA} mem -t 5 $(basename ${HOST}) ${READS}/${R1} ${READS}/${R2} | 
${SAMTOOLS} sort -o ${OUT}/${ALL} -

${SAMTOOLS} index ${OUT}/${ALL}
 
# pull out reads from analysis
echo filtering reads

${PYTHON} filter_host_align.py \
   --analysis-info $(basename $INTS) \
   --sim-bam ${OUT}/${ALL} \
   --output-bam ${OUT}/${CHIM} \
   --chimeric
   

${PYTHON} filter_host_align.py \
   --analysis-info $(basename $INTS) \
   --sim-bam ${OUT}/${ALL} \
   --output-bam ${OUT}/${DISC} \
   --discordant

${PYTHON} filter_host_align.py \
   --analysis-info $(basename $INTS) \
   --sim-bam ${OUT}/${ALL} \
   --output-bam ${OUT}/${BOTH} \
   --discordant \
   --chimeric

${SAMTOOLS} index ${OUT}/${CHIM}
${SAMTOOLS} index ${OUT}/${DISC}
${SAMTOOLS} index ${OUT}/${BOTH}

# clean up
rm $(basename ${HOST})*
rm $(basename ${SCRIPT})
rm $(basename $INTS)


