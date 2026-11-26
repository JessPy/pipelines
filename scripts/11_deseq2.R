#!/usr/bin/env Rscript
# [SUGERIDO] Expressão diferencial com DESeq2 a partir do featureCounts.
# Exemplo: GSE281650 (RPE-1, DMSO x etoposídeo, cinética de triptolide 0/20/40 min, 3 réplicas)
# Uso: Rscript scripts/11_deseq2.R counts/featureCounts_genes.txt results/deseq2
suppressPackageStartupMessages({ library(DESeq2); library(tidyverse); library(pheatmap) })
args <- commandArgs(TRUE)
fc_file <- ifelse(length(args) > 0, args[1], "counts/featureCounts_genes.txt")
outdir  <- ifelse(length(args) > 1, args[2], "results/deseq2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fc  <- read.delim(fc_file, comment.char = "#", check.names = FALSE)
cts <- as.matrix(fc[, -(1:6)]); rownames(cts) <- fc$Geneid
colnames(cts) <- sub("^bam/", "", sub("\\.bam$", "", colnames(cts)))

# nomes no padrão R1_RNAseq_DMSO_T0__SRRxxx -> réplica, tratamento, tempo
meta <- tibble(sample = colnames(cts)) |>
  mutate(base = sub("__.*", "", sample)) |>
  separate(base, c("rep", "assay", "treat", "time"), sep = "_", remove = FALSE) |>
  mutate(treat = factor(treat, c("DMSO", "ETOP")),
         time  = factor(time, c("T0", "T20", "T40")),
         rep   = factor(rep)) |>
  column_to_rownames("sample")
# se um GSM tiver vários SRR (lanes), some as colunas antes:
cts <- t(rowsum(t(cts), group = meta$base)); meta <- meta[!duplicated(meta$base), ]
rownames(meta) <- meta$base; cts <- cts[, rownames(meta)]

# 1) Efeito do etoposídeo no estado estacionário (T0)
m0  <- meta[meta$time == "T0", ]
dds <- DESeqDataSetFromMatrix(cts[, rownames(m0)], m0, ~ rep + treat)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "treat_ETOP_vs_DMSO", type = "apeglm")
write.csv(as.data.frame(res[order(res$padj), ]), file.path(outdir, "T0_ETOP_vs_DMSO.csv"))

# QC: PCA e correlação entre amostras
vsd <- vst(dds, blind = TRUE)
pdf(file.path(outdir, "QC_PCA_T0.pdf")); print(plotPCA(vsd, intgroup = c("treat", "rep"))); dev.off()
pdf(file.path(outdir, "QC_cor_T0.pdf")); pheatmap(cor(assay(vsd))); dev.off()
pdf(file.path(outdir, "MA_T0.pdf")); plotMA(res, ylim = c(-3, 3)); dev.off()

# 2) Decaimento após triptolide diferente entre tratamentos: interação tempo:tratamento (LRT)
ddsI <- DESeqDataSetFromMatrix(cts, meta, ~ rep + treat + time + treat:time)
ddsI <- ddsI[rowSums(counts(ddsI) >= 10) >= 3, ]
ddsI <- DESeq(ddsI, test = "LRT", reduced = ~ rep + treat + time)
write.csv(as.data.frame(results(ddsI)[order(results(ddsI)$padj), ]),
          file.path(outdir, "LRT_treat_x_time.csv"))
cat("DESeq2 concluído em", outdir, "\n")
