suppressPackageStartupMessages({library(data.table); library(DESeq2); library(ggplot2); library(ggrepel)})
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)

counts <- as.data.frame(fread("data/GSE159699/GSE159699_summary_count.star.txt.gz"))
rownames(counts) <- make.unique(as.character(counts[[1]])); counts <- counts[, -1]
counts <- counts[, sapply(counts, is.numeric)]
info <- data.frame(sample = colnames(counts), row.names = colnames(counts))
info$group <- factor(sub(".*-", "", info$sample))

# White-matter score (myelin genes minus neuron genes), used as a covariate in the sensitivity analysis
dds0 <- DESeqDataSetFromMatrix(round(as.matrix(counts)), colData = info, design = ~ 1)
dds0 <- dds0[rowSums(counts(dds0) >= 10) >= 5, ]
v <- assay(vst(dds0, blind = TRUE))
myelin <- c("MBP", "PLP1", "MOG", "MAG", "MOBP", "CLDN11")
neuron <- c("SNAP25", "RBFOX3", "SYT1", "SLC17A7", "CAMK2A")
info$wm_score <- as.numeric(scale(colMeans(v[intersect(myelin, rownames(v)), ]) -
                                  colMeans(v[intersect(neuron, rownames(v)), ])))

padj_cut <- 0.05; lfc_cut <- 0.58

run_de <- function(samples, design, label) {
  ci <- info[samples, ]; ci$group <- relevel(droplevels(ci$group), ref = "Old")
  dds <- DESeqDataSetFromMatrix(round(as.matrix(counts[, samples])), colData = ci, design = design)
  dds <- dds[rowSums(counts(dds) >= 10) >= min(table(ci$group)), ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", "AD", "Old"))
  shr <- lfcShrink(dds, coef = "group_AD_vs_Old", type = "apeglm", quiet = TRUE)
  df <- data.frame(gene = rownames(res), baseMean = res$baseMean, log2FC = res$log2FoldChange,
                   log2FC_shrunk = shr$log2FoldChange, stat = res$stat, pvalue = res$pvalue, padj = res$padj)
  df <- df[order(df$padj), ]
  df$status <- ifelse(!is.na(df$padj) & df$padj < padj_cut & df$log2FC >  lfc_cut, "Up in AD",
               ifelse(!is.na(df$padj) & df$padj < padj_cut & df$log2FC < -lfc_cut, "Down in AD", "NS"))
  write.csv(df, paste0("results/tables/DE_", label, ".csv"), row.names = FALSE)
  d <- df[!is.na(df$padj), ]; top <- head(d[d$status != "NS", ], 15)
  ggsave(paste0("results/figures/volcano_", label, ".png"),
    ggplot(d, aes(log2FC, -log10(padj), colour = status)) + geom_point(alpha = 0.6, size = 1.2) +
      scale_colour_manual(values = c("Down in AD" = "#2166ac", "NS" = "grey75", "Up in AD" = "#b2182b")) +
      geom_text_repel(data = top, aes(label = gene), size = 3, max.overlaps = 30, colour = "black") +
      geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
      geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
      theme_bw() + ggtitle(paste("AD vs Old:", label), subtitle = paste(table(ci$group)[["AD"]], "AD vs", table(ci$group)[["Old"]], "Old")),
    width = 7, height = 6, dpi = 300)
  cat("\n=====", label, "=====\n"); print(table(df$status))
  df
}

exclude <- c("22-2T-AD", "25-5T-AD")
ad_old  <- rownames(info)[info$group %in% c("AD", "Old")]

primary   <- run_de(setdiff(ad_old, exclude), ~ group,            "primary_excl2")
sensitiv  <- run_de(ad_old,                   ~ wm_score + group, "sensitivity_all_adj")

cat("\n--- Top 15 genes, primary analysis ---\n")
print(head(primary[, c("gene", "log2FC", "padj", "status")], 15), row.names = FALSE)

sig <- function(x) x$gene[x$status != "NS"]
cat("\n--- Overlap of significant genes ---\n")
cat("Primary:", length(sig(primary)), " Sensitivity:", length(sig(sensitiv)),
    " Shared:", length(intersect(sig(primary), sig(sensitiv))), "\n")
