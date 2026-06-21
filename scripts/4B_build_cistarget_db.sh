#!/bin/bash
########## BUILD CUSTOM DB: topics + DARs ##########
set -euo pipefail

########################
# Parameters (passed in by Snakemake, in this order)
# These replace the original "EDIT THESE PER PROJECT" block.
########################
CTX_DB_BASE="$1"        # base cistarget folder (cbust, motif collection, helper repo)
DB_DIR="$2"             # where to build this DB
DB_PREFIX="$3"          # name prefix for the DB files
TOPICS_TOP_DIR="$4"     # region sets: Topics_top3k
TOPICS_OTSU_DIR="$5"    # region sets: Topics_otsu
DARS_DIR="$6"           # region sets: DAR BEDs
GENOME_FA="$7"          # genome FASTA
CHROMSIZES="$8"         # chrom sizes
N_THREADS="$9"          # number of CPUs

create_cistarget_databases_dir="$CTX_DB_BASE/create_cisTarget_databases"
MOTIF_DIR="$CTX_DB_BASE/aertslab_motif_collection/v10nr_clust_public/singletons"
CBUST="$CTX_DB_BASE/cbust"
chmod u+x "$CBUST"

########################
# BUILD (no edits below)
########################

mkdir -p "$DB_DIR"
cd "$DB_DIR"

echo "1) Make consensus BED from topics + DARs..."
cat "$TOPICS_TOP_DIR"/*.bed "$TOPICS_OTSU_DIR"/*.bed "$DARS_DIR"/*.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > consensus_regions.bed

echo "Consensus regions:"
wc -l consensus_regions.bed

echo "Filtering consensus_regions.bed to chromosomes present in chromsizes..."
# consensus_regions.bed: chr / start / end     CHROMSIZES: chr / size
awk 'NR==FNR { ok[$1]; next } ($1 in ok)' "$CHROMSIZES" consensus_regions.bed \
  > consensus_regions.filtered.bed
mv consensus_regions.filtered.bed consensus_regions.bed
echo "Filtered consensus regions:"
wc -l consensus_regions.bed

echo "2) Pad by 1kb..."
bedtools slop -i consensus_regions.bed -g "$CHROMSIZES" -b 1000 \
  > consensus_regions.pad1kb.bed

chmod u+x "$create_cistarget_databases_dir"/*.sh "$create_cistarget_databases_dir"/*.py

echo "3) Build FASTA with 1kb background padding..."
"$create_cistarget_databases_dir/create_fasta_with_padded_bg_from_bed.sh" \
  "$GENOME_FA" \
  "$CHROMSIZES" \
  "consensus_regions.bed" \
  "genome.consensus.with_1kb_bg_padding.fa" \
  1000 \
  yes

echo "4) List motifs..."
ls "$MOTIF_DIR" > motifs.txt

echo "5) Score motifs..."
"$create_cistarget_databases_dir/create_cistarget_motif_databases.py" \
  -f "genome.consensus.with_1kb_bg_padding.fa" \
  -M "$MOTIF_DIR" \
  -m "motifs.txt" \
  -o "$DB_PREFIX" \
  --bgpadding 1000 \
  -t "$N_THREADS" \
  -c "$CBUST"

echo "Custom DB built in: $DB_DIR"
