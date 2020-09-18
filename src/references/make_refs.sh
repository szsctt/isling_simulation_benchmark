
mkdir -p ../../data/references
cd ../../data/references

# get human reference and split into chromosomes

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38.fa

# https://stackoverflow.com/questions/11818495/split-a-fasta-file-and-rename-on-the-basis-of-first-line
perl -pe 's/^>(\w+?) .+/>$1/' GRCh38.fa | awk '/^>chr/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}'

# get AAV2 sequence (NCBI Reference Sequence: NC_001401.2)
conda activate eutils
esearch -db nucleotide -query NC_001401.2 | efetch -format fasta > NC_001401.2.fa

