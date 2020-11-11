#!/bin/bash
set -e
CONFIG="$HOME/intvi_simulation-experiments/verse_test/config/sim-chr1.yml"
SIMDIR="$HOME/intvi_simulation-experiments/intvi_simulation"

cd $SIMDIR
snakemake \
 --configfile $CONFIG \
 --cores 8 \
 --jobs 3 \
 --restart-times 3 \
 --keep-going \
 --resources mem_mb=28000 \
 --use-singularity \
 --rerun-incomplete \
 --default-remote-provider "GS" \
 --default-remote-prefix sjs_snakemake_test \
 --notemp

# download simulated data from bucket
DATADIR="$HOME/intvi_simulation-experiments/verse_test/data/test_chr1_cloud/"

mkdir -p $DATADIR
gsutil cp gs://sjs_snakemake_test/verse_test/chr1_sim_data/test/sim_reads/*fq $DATADIR

# index references
REFDIR="$HOME/intvi_simulation-experiments/verse_test/refs/"
mkdir -p $REFDIR

# just get them from bucket
gsutil cp -r gs://sjs_snakemake_test/out/experiment1_OTC_chr1/easier-harder/verse_references/OTC $REFDIR
gsutil cp -r gs://sjs_snakemake_test/out/experiment1_OTC_chr1/easier-harder/verse_references/hg38 $REFDIR


OUTDIR="$HOME/intvi_simulation-experiments/verse_test/test-chr1-cloud/"
mkdir -p ${OUTDIR}
cd ${OUTDIR}
rm -rf *
CONFIG="$HOME/intvi_simulation-experiments/verse_test/config/test-chr1-cloud.txt"

# pull singularity container
if [ ! -e verse_1.sif ]; then
	singularity pull --name verse_1.sif docker://szsctt/verse:1
fi

cp $CONFIG .

# run container
singularity exec \
 -B$(realpath $DATADIR) -B$(realpath $REFDIR) \
 verse_1.sif  \
 perl /var/work/VirusFinder2.0/VirusFinder.pl -c $(basename $CONFIG) &> verse_test.log

