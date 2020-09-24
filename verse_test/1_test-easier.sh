#!/bin/bash
set -e

#SBATCH --job-name=virusFinderTest   			# Job name
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1                  		# Run on a single CPU
#SBATCH --mem-per-cpu=20gb                     	# Job memory request
#SBATCH --time=24:00:00               		# Time limit hrs:min:sec
#SBATCH --output=../logs/virusFinderTest_%j.log   	# Standard output and error log

eval "$(conda shell.bash hook)"
conda activate virusFinder2
module load jdk/1.7.0_65

CONFIG="/scratch1/sco305/intvi_simulation-experiments/verse_test/config/test-easier.txt"
VIRUSFINDER="/scratch1/sco305/intvi_external-data/test-virusfinder2/tools/VirusFinder2.0/VirusFinder.pl"

OUTDIR="/scratch1/sco305/intvi_simulation-experiments/verse_test/out/test-easier"

mkdir -p ${OUTDIR}
cd ${OUTDIR}
rm -rf *

perl ${VIRUSFINDER} -c ${CONFIG}

#### step 1 ####

# run preprocess script
# perl -d -I /scratch1/sco305/intvi_external-data/test-virusfinder2/tools/VirusFinder2.0  /scratch1/sco305/intvi_external-data/test-virusfinder2/tools/VirusFinder2.0/preprocess.pl -c /scratch1/sco305/intvi_external-data/test-virusfinder2/config/pAAV2-test-config.txt -o /scratch1/sco305/intvi_external-data/out-virusFinder2/virusFinder2-test/step1

## bowtie2 command line: this is just the default settings

# /scratch1/sco305/conda/envs/virusFinder2/bin/bowtie2 -p 8 -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x /scratch1/sco305/intvi_external-data/data/references/hg19 -1 /scratch1/sco305/intvi_cmri/data/reads/FRG_OTC/FRG203_C73YB_CGAACTTA_L001_R1.fastq.gz -2 /scratch1/sco305/intvi_cmri/data/reads/FRG_OTC/FRG203_C73YB_CGAACTTA_L001_R2.fastq.gz -S /scratch1/sco305/intvi_external-data/out-virusFinder2/virusFinder2-test/step1/alignment.sam

## samtools to convert sam to bam
# samtools view -bS /scratch1/sco305/intvi_external-data/out-virusFinder2/virusFinder2-test/step1/alignment.sam -o /scratch1/sco305/intvi_external-data/out-virusFinder2/virusFinder2-test/step1/alignment.bam

#### step 2 ####

## run find virus script
# perl -d -I /scratch1/sco305/intvi_external-data/test-virusfinder2/tools/VirusFinder2.0  /scratch1/sco305/intvi_external-data/test-virusfinder2/tools/VirusFinder2.0/detect_virus.pl -c /scratch1/sco305/intvi_external-data/test-virusfinder2/config/pAAV2-test-config.txt -o /scratch1/sco305/intvi_external-data/out-virusFinder2/virusFinder2-test/step2

# /scratch1/sco305/intvi_external-data/test-virusfinder2/tools/VirusFinder2.0/bin/blat /scratch1/sco305/intvi_external-data/data/references/pAAV2-OTC.fa -minIdentity=80 chopped_unmapped.2.fa chopped_unmapped.2.psl

# perl -I /scratch1/sco305/intvi_external-data/test-virusfinder2/tools/VirusFinder2.0  /scratch1/sco305/intvi_external-data/test-virusfinder2/tools/trinityrnaseq_r2013-02-16/Trinity.pl --seqType fa  --JM 4G --single blat_out_candidate_singlelane.fa --min_contig_length 300  --output trinity_output --CPU 8 --bfly_opts "-V 10 --stderr"

