
mkdir -p ../../data/references
cd ../../data/references

# get human reference and split into chromosomes
mkdir GRCh38 && cd GRCh38

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
 

# https://stackoverflow.com/questions/11818495/split-a-fasta-file-and-rename-on-the-basis-of-first-line
mkdir chrs
perl -pe 's/^>(\w+?) .+/>$1/' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | awk '/^>chr/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}'
mv chr* chrs

# get human gff from ENCODE
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gff3.gz

gunzip gencode.v35.annotation.gff3.gz

# filter genes and exons
awk 'match ($3, /gene/) || match($1, /^#/)' gencode.v35.annotation.gff3 > hg38.gencode.v35.annotation.genes.gff3

awk 'match ($3, /exon/) || match($1, /^#/)' gencode.v35.annotation.gff3 > hg38.gencode.v35.annotation.exons.gff3


# get data repo for vifi
#https://drive.google.com/drive/folders/0ByYcg0axX7udeGFNVWtaUmxrOFk

# get mappability/uniqueness
#https://wiki.bits.vib.be/index.php/Create_a_mappability_track
GEM="gem.sif"
module load singularity
if [ ! -e $GEM ]; then
	singularity pull --name $GEM "docker://szsctt/gem:1"
fi

mkdir -p hg38_mappability/idx
srun --time 2:00:00 -c8 --mem 30gb \
singularity exec $GEM \
gem-indexer -T 8 -c dna -i GRCh38.fa -o hg38_mappability/idx/hg38 

srun --time 2:00:00 -c8 \
singularity exec $GEM \
gem-mappability -T 8 -I hg38_mappability/idx/hg38.gem -l 35 -o hg38_mappability/hg38_35
 
singularity exec $GEM \
gem-2-wig -I hg38_mappability/idx/hg38.gem -i hg38_35.mappability -o hg38_mappability/hg38_35

singularity exec $GEM \
wigToBigWig hg38_mappability/hg38_35.wig hg38_mappability/hg38.sizes hg38_mappability/hg38_35.bw

# ENCODE blacklisted regions - use this instead of Duke Excluded regions
# https://genome.sph.umich.edu/wiki/MappabilityScores
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
gunzip ENCFF356LFX.bed.gz

# get segmental duplication
# https://humanparalogy.gs.washington.edu/code/WGAC_HOWTO.pdf
#wget https://humanparalogy.gs.washington.edu/build38/data/GRCh38GenomicSuperDup.tab
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
gunzip genomicSuperDups.txt.gz
cut -f2- genomicSuperDups.txt > genomicSuperDups.bed

# hg38 centromeres
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
gunzip centromeres.txt.gz
cut -f2- centromeres.txt > centromeres.bed

## Amplicon archtect fork for hg38
#https://github.com/jluebeck/AmpliconArchitect
#https://drive.google.com/drive/folders/18T83A12CfipB0pnGsTs3s-Qqji6sSvCu

# conserved regions in amplicon architect - not sure what this file is exactly...
wget https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed


cd ..
mkdir AAV2 && cd AAV2
# get AAV2 sequence (NCBI Reference Sequence: NC_001401.2)
conda activate eutils
esearch -db nucleotide -query NC_001401.2 | efetch -format fasta > NC_001401.2.fa

