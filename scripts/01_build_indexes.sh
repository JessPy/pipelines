#!/usr/bin/env bash
# Constrói os índices: Bowtie1 (hg38 + dm3 combinados) e HISAT2 (GRCh38 Ensembl).
# Precisa de ~16 GB RAM (Bowtie) e ~8 GB (HISAT2 sem SNPs/splice sites) ; ~2-3 h.
set -euo pipefail
THREADS=${THREADS:-8}
cd ref

# --- Genoma combinado humano + Drosophila (spike-in) ---------------------
# Prefixo "dm3_" nos cromossomos de mosca para separá-los depois. [INFERIDO]
if [[ ! -s hg38_dm3.fa ]]; then
  zcat hg38.fa.gz  | awk '/^>/{keep = ($1 !~ /_/)} keep' > hg38.main.fa   # só chr1-22,X,Y,M
  zcat dm3.fa.gz   | sed 's/^>/>dm3_/' > dm3.pref.fa
  cat hg38.main.fa dm3.pref.fa > hg38_dm3.fa
  samtools faidx hg38_dm3.fa
  cut -f1,2 hg38_dm3.fa.fai | grep -v '^dm3_' > hg38.chrom.sizes
fi
mkdir -p bowtie1
[[ -s bowtie1/hg38_dm3.1.ebwtl || -s bowtie1/hg38_dm3.1.ebwt ]] || \
  bowtie-build --threads "$THREADS" hg38_dm3.fa bowtie1/hg38_dm3

# --- HISAT2 (Ensembl GRCh38 primário) -------------------------------------
mkdir -p hisat2
if [[ ! -s hisat2/GRCh38.1.ht2 ]]; then
  zcat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > GRCh38.fa
  hisat2-build -p "$THREADS" GRCh38.fa hisat2/GRCh38
fi
[[ -s Homo_sapiens.GRCh38.108.gtf ]] || gunzip -k Homo_sapiens.GRCh38.108.gtf.gz

# --- GTF com prefixo "chr" para usar com coordenadas UCSC (metaplots) ----
[[ -s genes.hg38.chr.bed ]] || \
  awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && $9 ~ /gene_biotype "protein_coding"/ {
       match($9,/gene_id "[^"]+"/); id=substr($9,RSTART+9,RLENGTH-10);
       c=($1=="MT")?"chrM":"chr"$1; print c,$4-1,$5,id,".",$7}' Homo_sapiens.GRCh38.108.gtf \
  | sort -k1,1 -k2,2n > genes.hg38.chr.bed
echo "Índices prontos."
