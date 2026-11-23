suppressPackageStartupMessages({library(data.table); library(DESeq2); library(ggplot2); library(ggrepel); library(pheatmap)})
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

counts <- as.data.frame(fread("data/GSE159699/GSE159699_summary_count.star.txt.gz"))
rownames(counts) <- make.unique(as.character(counts[[1]])); counts <- counts[, -1]
counts <- counts[, sapply(counts, is.numeric)]

info <- data.frame(sample = colnames(counts), row.names = colnames(counts))
info$group  <- sub(".*-", "", info$sample)
info$suffix <- sub("^[0-9]+-[0-9]+([AT])-.*$", "\\1", info$sample)
info$donor  <- sub("^[0-9]+-([0-9]+)[AT]-.*$", "\\1", info$sample)
print(table(info$group, info$suffix))

dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), colData = info, design = ~ 1)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
vsd <- vst(dds, blind = TRUE)

pca <- plotPCA(vsd, intgroup = c("group", "suffix"), returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
ggsave("results/figures/00_PCA_explore.png",
  ggplot(pca, aes(PC1, PC2, colour = group, shape = suffix, label = name)) +
    geom_point(size = 3) + geom_text_repel(size = 2.5, max.overlaps = 40) +
    xlab(paste0("PC1: ", pv[1], "%")) + ylab(paste0("PC2: ", pv[2], "%")) + theme_bw(),
  width = 8, height = 6, dpi = 300)

cm <- cor(assay(vsd))
pheatmap(cm, annotation_col = info[, c("group", "suffix")], fontsize = 6,
         filename = "results/figures/00_sample_correlation.png", width = 9, height = 8)

diag(cm) <- NA
nn <- apply(cm, 1, function(x) names(which.max(x)))
cat("\nMost similar sample to each sample:\n")
print(data.frame(sample = names(nn), nearest = nn, row.names = NULL))
