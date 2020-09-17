#!/bin/bash


eval "$(conda shell.bash hook)"

# requirements from git hub:
#bowtie2 (v2.2.6), samtools (v1.9), pysam (v0.8.4), and bedtools (v2.26.0) installed. Also, Polyidus also uses python packages pandas (v0.22.0), numpy (1.18.1), and psutil (v5.4.3). You need to run the following commands:

# note thatusing the following line results in numpy 1.15 rather than 1.18 - using 1.18 results in a conflict
mamba create -n polyidus -c bioconda -c conda-forge bowtie2=2.2.6 samtools=1.9 bedtools=2.26.0 pysam=0.8.4 pandas=0.22.0 psutil=5.4.3 

 




