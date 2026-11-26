#!/usr/bin/env bash
# Baixa FASTQ do SRA a partir de uma lista de SRX/SRR/GSM e as referências.
# Uso: bash scripts/00_download.sh config/samples_chor.tsv
#      (tsv: nome<TAB>acesso SRX ou SRR)
set -euo pipefail
SHEET=${1:?"informe o samples.tsv"}
THREADS=${THREADS:-8}
mkdir -p raw ref

# 1) FASTQ  ----------------------------------------------------------------
while IFS=$'\t' read -r name acc; do
  [[ -z "$name" || "$name" == \#* ]] && continue
  # um SRX pode ter vários SRR (réplicas fundidas pelos autores)
  for srr in $(curl -s "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=${acc}" | awk -F, 'NR>1{print $1}'); do
    if ls raw/${name}__${srr}*.fastq.gz >/dev/null 2>&1; then continue; fi
    prefetch -O raw/sra "$srr"
    fasterq-dump --split-files -e "$THREADS" -O raw/tmp "raw/sra/$srr"
    for f in raw/tmp/${srr}*.fastq; do
      b=$(basename "$f"); pigz -p "$THREADS" -c "$f" > "raw/${name}__${b}.gz"; rm "$f"
    done
  done
done < "$SHEET"

# 2) Referências -----------------------------------------------------------
cd ref
# hg38 (UCSC, nomes chr1..) e dm3 (spike-in Drosophila, versão usada pelos autores)
[[ -s hg38.fa.gz ]] || wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
[[ -s dm3.fa.gz  ]] || wget -q https://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/dm3.fa.gz
# Ensembl 108 (anotação citada no GSE281650) - genoma primário + GTF
[[ -s Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ]] || \
  wget -q https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
[[ -s Homo_sapiens.GRCh38.108.gtf.gz ]] || \
  wget -q https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
echo "Download concluído."
