#!/bin/bash
set -e

#SBATCH --job-name=virusFinderTest   			# Job name
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1                  		# Run on a single CPU
#SBATCH --mem-per-cpu=20gb                     	# Job memory request
#SBATCH --time=24:00:00               		# Time limit hrs:min:sec
#SBATCH --output=../logs/virusFinderTest_%j.log   	# Standard output and error log

module load singularity

HOST="/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/data/references/GRCh38.fa"
VIRUS="/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/data/references/OTC-vec_rAAV-genome-only.fa"
READDIR="/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/verse_test/reads/"
CONFIG="/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/verse_test/config/test-chr1.txt"
OUTDIR="/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/verse_test/test-chr1"
REFDIR="/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/verse_test/references"

mkdir -p ${REFDIR}
cd ${REFDIR}

if [ ! -e verse_1.sif ]; then
	singularity pull docker://szsctt/verse:1
fi

rsync $VIRUS .
rsync $HOST .

if [ ! -e GRCh38.1.bt2 ]; then
	srun --time 24:00:00 --mem 20gb \
	singularity exec -B$(dirname $HOST) verse_1.sif \
	bowtie2-build $HOST GRCh38
fi

if [ ! -e GRCh38.nhr ]; then
srun --time 2:00:00 --mem 20gb \
singularity exec -B$(dirname $HOST) verse_1.sif  \
makeblastdb -in $HOST -dbtype nucl -out GRCh38
fi

if [ ! -e OTC-vec_rAAV-genome-only.nhr ]; then
	srun --time 2:00:00 --mem 20gb \
	singularity exec -B$(dirname $VIRUS) verse_1.sif  \
	makeblastdb -in $VIRUS -dbtype nucl -out OTC-vec_rAAV-genome-only
fi

mkdir -p ${OUTDIR}
cd ${OUTDIR}
rm -rf *

if [ ! -e verse_1.sif ]; then
	singularity pull docker://szsctt/verse:1
fi


srun --time 24:00:00 --mem 200gb -c8 \
singularity exec \
 -B${READDIR} -B$REFDIR \
 verse_1.sif  \
 perl /var/work/VirusFinder2.0/VirusFinder.pl -c ${CONFIG}

