#!/usr/bin/env bash
# ChIP-seq / ChOR-seq / EdU-seq com spike-in de Drosophila
# FASTQ -> Bowtie1 (-m 1 --best) em hg38+dm3 -> dedup -> separa humano/mosca
#       -> fatores CPM e RRPM -> bedGraph/bigWig (bedtools genomecov)
# Uso: bash scripts/20_chorseq.sh config/samples_chor.tsv
#   variáveis: THREADS=8  PE=1 (usar pares; 0 = só R1)  FRAGLEN=250 (extensão em single-end)
set -euo pipefail
SHEET=${1:?"samples.tsv"}
THREADS=${THREADS:-8}; PE=${PE:-1}; FRAGLEN=${FRAGLEN:-250}
IDX=ref/bowtie1/hg38_dm3
SIZES=ref/hg38.chrom.sizes
mkdir -p qc/bowtie bam tracks

HUMAN_CHR=$(cut -f1 "$SIZES" | tr '\n' ' ')

for name in $(grep -v '^#' "$SHEET" | cut -f1); do
  parts=()
  for r1 in raw/${name}__*.fastq.gz; do
    [[ "$r1" == *_2.fastq.gz ]] && continue
    id=$(basename "$r1" .fastq.gz); id=${id%_1}
    out=bam/$id.raw.bam; parts+=("$out")
    [[ -s $out ]] && continue
    r2=${r1/_1.fastq.gz/_2.fastq.gz}
    # [DESCRITO] Bowtie 1, parâmetros "-m1 --best"
    if [[ $PE == 1 && -s $r2 && "$r1" == *_1.fastq.gz ]]; then
      # [INFERIDO] modo pareado; -X 1000 aceita fragmentos de 250-500 bp (sonicação descrita)
      bowtie -p "$THREADS" -m 1 --best -X 1000 -S "$IDX" -1 <(zcat "$r1") -2 <(zcat "$r2") \
        2> qc/bowtie/$id.log | samtools view -b -F 4 - > $out.tmp
    else
      bowtie -p "$THREADS" -m 1 --best -S "$IDX" <(zcat "$r1") \
        2> qc/bowtie/$id.log | samtools view -b -F 4 - > $out.tmp
    fi
    mv $out.tmp $out
  done

  # junta lanes/réplicas do mesmo GSM [DESCRITO: "Merge BAM Files"]
  if [[ ! -s bam/$name.dedup.bam ]]; then
    samtools merge -f -@ "$THREADS" -n bam/$name.merged.bam "${parts[@]}"
    # [INFERIDO] remoção de duplicatas: os autores usaram RmDup (Galaxy) ou UMI-tools.
    # Os UMIs ficam no nome da leitura na saída do BCL Convert, mas o SRA normalmente
    # substitui o nome -> use samtools markdup -r (equivalente moderno do rmdup).
    samtools fixmate -@ "$THREADS" -m bam/$name.merged.bam - \
      | samtools sort -@ "$THREADS" - \
      | samtools markdup -@ "$THREADS" -r -s - bam/$name.dedup.bam 2> qc/bowtie/$name.markdup.txt
    samtools index bam/$name.dedup.bam
    rm bam/$name.merged.bam
  fi

  # separa humano e Drosophila
  samtools view -@ "$THREADS" -b bam/$name.dedup.bam $HUMAN_CHR > bam/$name.hg38.bam
  samtools index bam/$name.hg38.bam
  # leituras (ou pares) únicas na mosca: base do fator de spike-in
  if samtools view -c -f 1 bam/$name.dedup.bam | grep -qv '^0$'; then FLAG="-f 64"; else FLAG=""; fi
  n_dm=$(samtools idxstats bam/$name.dedup.bam | awk '$1 ~ /^dm3_/ {s+=$3} END{print s+0}')
  n_hs=$(samtools view -c $FLAG bam/$name.hg38.bam)
  [[ -n $FLAG ]] && n_dm=$(( n_dm / 2 ))
  echo -e "$name\t$n_hs\t$n_dm" >> qc/spikein_counts.tsv

  # [DESCRITO] CPM: bedtools genomecov 2.30.0, escala 1e6 / leituras humanas
  cpm=$(awk -v n=$n_hs 'BEGIN{print 1e6/n}')
  # [DESCRITO 2018] RRPM (Orlando et al. 2014): escala 1e6 / leituras únicas de Drosophila
  rrpm=$(awk -v n=$n_dm 'BEGIN{ if(n>0) print 1e6/n; else print 0 }')
  if [[ -n $FLAG ]]; then
    samtools sort -n -@ "$THREADS" bam/$name.hg38.bam -o bam/$name.hg38.nsort.bam
    GC_IN="-ibam bam/$name.hg38.nsort.bam -pc"          # cobre o fragmento inteiro
  else
    GC_IN="-ibam bam/$name.hg38.bam -fs $FRAGLEN"       # estende leituras (2018: 250/500 bp)
  fi
  for norm in cpm rrpm; do
    sc=$([[ $norm == cpm ]] && echo $cpm || echo $rrpm)
    [[ $sc == 0 ]] && continue
    bedtools genomecov $GC_IN -bg -scale "$sc" | LC_ALL=C sort -k1,1 -k2,2n > tracks/$name.$norm.bedGraph
    bedGraphToBigWig tracks/$name.$norm.bedGraph "$SIZES" tracks/$name.$norm.bw
  done
  rm -f bam/$name.hg38.nsort.bam
done

column -t qc/spikein_counts.tsv
echo "Tracks em tracks/*.cpm.bw e tracks/*.rrpm.bw"
