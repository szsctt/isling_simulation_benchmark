#!/bin/bash

# get data repo for vifi
#https://drive.google.com/drive/folders/0ByYcg0axX7udeGFNVWtaUmxrOFk

set -e

eval "$(conda shell.bash hook)"
#conda activate sim_isling

cd data/references

# get human reference and split into chromosomes
mkdir GRCh38 && cd GRCh38

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# https://stackoverflow.com/questions/11818495/split-a-fasta-file-and-rename-on-the-basis-of-first-line
mkdir chrs
perl -pe 's/^>(\w+?) .+/>$1/' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | awk '/^>chr/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}'
mv chr*.fa chrs

# get human gff from ENCODE
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gff3.gz

#gunzip gencode.v35.annotation.gff3.gz

# filter genes and exons
#awk 'match ($3, /gene/) || match($1, /^#/)' gencode.v35.annotation.gff3 > hg38.gencode.v35.annotation.genes.gff3

#awk 'match ($3, /exon/) || match($1, /^#/)' gencode.v35.annotation.gff3 > hg38.gencode.v35.annotation.exons.gff3

cd ..
# get AAV2 sequence (NCBI Reference Sequence: NC_001401.2)
#conda activate eutils
esearch -db nucleotide -query NC_001401.2 | efetch -format fasta > NC_001401.2.fa

# get HBV sequence (NCBI Refernece sequence: NC_003977.2)
esearch -db nucleotide -query NC_003977.2 | efetch -format fasta > NC_003977.2.fa

# get HBV sequence (NCBI Refernece sequence: NC_027779.1)
esearch -db nucleotide -query NC_027779.1 | efetch -format fasta > NC_027779.1.fa

