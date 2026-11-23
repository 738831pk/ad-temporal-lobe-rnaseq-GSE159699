suppressPackageStartupMessages({
  library(data.table); library(DESeq2); library(ggplot2); library(ggrepel)
})
dir.create("results/figures/final", recursive = TRUE, showWarnings = FALSE)
wrap_txt <- function(x, w = 42) sapply(x, function(s) paste(strwrap(s, width = w), collapse = "\n"))

# ---- Fig 1: PCA QC (all 30 samples, highlighting white-matter-rich cluster) ----
counts <- as.data.frame(fread("data/GSE159699/GSE159699_summary_count.star.txt.gz"))
rownames(counts) <- make.unique(as.character(counts[[1]])); counts <- counts[, -1]
counts <- counts[, sapply(counts, is.numeric)]
info <- data.frame(sample = colnames(counts), row.names = colnames(counts))
info$condition <- factor(sub(".*-", "", info$sample), levels = c("Young", "Old", "AD"))
wm_samples <- c("22-2T-AD", "25-5T-AD", "8-15A-Young", "3-17T-Young", "5-18T-Young")
info$qc <- ifelse(info$sample %in% wm_samples, "White-matter-rich", "Grey-matter-rich")

dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), colData = info, design = ~ 1)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
pca <- plotPCA(vst(dds, blind = TRUE), intgroup = c("condition", "qc"), returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))

p1 <- ggplot(pca, aes(PC1, PC2, colour = condition, shape = qc)) +
  geom_point(size = 3.2, alpha = 0.9) +
  geom_text_repel(data = pca[pca$qc == "White-matter-rich", ], aes(label = name),
                  size = 3, show.legend = FALSE, box.padding = 0.4) +
  scale_colour_manual(values = c("Young" = "#1b9e77", "Old" = "#377eb8", "AD" = "#d95f02")) +
  scale_shape_manual(values = c("Grey-matter-rich" = 16, "White-matter-rich" = 17)) +
  labs(title = "PCA of GSE159699 Bulk RNA-seq (n = 30)",
       subtitle = "PC1 separates 5 white-matter-rich dissections (2 AD excluded from primary DE)",
       x = paste0("PC1 (", pv[1], "% variance)"), y = paste0("PC2 (", pv[2], "% variance)"),
       colour = "Group", shape = "Tissue Profile") +
  theme_bw(base_size = 12) + theme(legend.position = "right")
ggsave("results/figures/final/fig1_pca_qc.png", p1, width = 8, height = 5.5, dpi = 300)

# ---- Fig 2: Volcano plot (Primary AD vs Old, 10 vs 10) ----
de <- read.csv("results/tables/DE_primary_excl2.csv")
de <- de[!is.na(de$padj), ]
top_lab <- head(de[de$status != "NS", ], 12)
n_up <- sum(de$status == "Up in AD"); n_dn <- sum(de$status == "Down in AD")

p2 <- ggplot(de, aes(log2FC, -log10(padj), colour = status)) +
  geom_point(alpha = 0.6, size = 1.4) +
  scale_colour_manual(values = c("Down in AD" = "#2166ac", "NS" = "grey78", "Up in AD" = "#b2182b")) +
  geom_text_repel(data = top_lab, aes(label = gene), size = 3.2, colour = "black",
                  max.overlaps = 30, box.padding = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", colour = "grey40") +
  labs(title = "Differential Expression: AD vs Old Controls (Primary, 10 vs 10)",
       subtitle = paste0("padj < 0.05 & |log2FC| > 0.58  |  Up: ", n_up, "  |  Down: ", n_dn),
       x = "log2 Fold Change (AD vs Old)", y = "-log10(adjusted p-value)", colour = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave("results/figures/final/fig2_volcano.png", p2, width = 7.5, height = 6, dpi = 300)

# ---- Fig 3: GO BP ORA (Downregulated in AD) ----
ora <- read.csv("results/tables/GO_ORA_down_in_AD.csv")
ora_top <- head(ora[order(ora$p.adjust), ], 15)
ora_top$GeneRatioNum <- sapply(strsplit(ora_top$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
ora_top$Label <- factor(wrap_txt(ora_top$Description, 44), levels = rev(wrap_txt(ora_top$Description, 44)))

p3 <- ggplot(ora_top, aes(GeneRatioNum, Label, size = Count, colour = p.adjust)) +
  geom_point() +
  scale_colour_gradient(low = "#b2182b", high = "#2166ac", trans = "log10") +
  labs(title = "GO Biological Process ORA: Downregulated in AD",
       subtitle = "Top 15 of 44 enriched terms (0 terms enriched for upregulated genes)",
       x = "Gene Ratio", y = NULL, colour = "adj. p-value", size = "Gene Count") +
  theme_bw(base_size = 11)
ggsave("results/figures/final/fig3_go_ora.png", p3, width = 8.5, height = 6.5, dpi = 300)

# ---- Fig 4: GSEA GO BP (Top Activated & Suppressed) ----
gsea <- read.csv("results/tables/GSEA_GO_BP.csv")
top_act <- head(gsea[gsea$NES > 0, ][order(gsea$p.adjust[gsea$NES > 0]), ], 10)
top_sup <- head(gsea[gsea$NES < 0, ][order(gsea$p.adjust[gsea$NES < 0]), ], 10)
gplot <- rbind(top_act, top_sup)
gplot$Direction <- ifelse(gplot$NES > 0, "Activated in AD", "Suppressed in AD")
gplot <- gplot[order(gplot$NES), ]
gplot$Label <- factor(wrap_txt(gplot$Description, 44), levels = wrap_txt(gplot$Description, 44))

p4 <- ggplot(gplot, aes(NES, Label, size = setSize, colour = p.adjust)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point() +
  scale_colour_gradient(low = "#b2182b", high = "#2166ac", trans = "log10") +
  labs(title = "GSEA GO Biological Process (Primary AD vs Old)",
       subtitle = paste0("Top 10 activated and top 10 suppressed terms (", nrow(gsea), " total significant)"),
       x = "Normalized Enrichment Score (NES)", y = NULL, colour = "adj. p-value", size = "Set Size") +
  theme_bw(base_size = 11)
ggsave("results/figures/final/fig4_gsea.png", p4, width = 9, height = 7, dpi = 300)

# ---- Fig 5: Targeted Neuroinflammation & Glial Markers ----
panel <- data.frame(
  gene = c("TREM2","TYROBP","CD68","AIF1","CX3CR1","P2RY12",
           "C1QA","C1QB","C1QC","C3",
           "HLA-DRA","HLA-DRB1","B2M","CD74",
           "GFAP","CHI3L1","IL1B","APOE","CLU"),
  category = c(rep("Microglia / DAM", 6), rep("Complement", 4),
               rep("MHC / Antigen", 4), rep("Astrocyte / Cytokine", 3), rep("Lipid / AD Risk", 2))
)
imm <- merge(panel, de, by = "gene")
imm$sig <- ifelse(imm$padj < 0.05, "padj < 0.05 (below LFC cutoff)", "Not significant (padj >= 0.05)")
imm <- imm[order(imm$log2FC), ]
imm$gene <- factor(imm$gene, levels = imm$gene)

p5 <- ggplot(imm, aes(log2FC, gene, colour = sig)) +
  geom_vline(xintercept = 0, colour = "grey30") +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", colour = "grey65") +
  geom_segment(aes(x = 0, xend = log2FC, y = gene, yend = gene), linewidth = 0.7) +
  geom_point(size = 3) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("padj < 0.05 (below LFC cutoff)" = "#d95f02",
                                 "Not significant (padj >= 0.05)" = "#4575b4")) +
  labs(title = "Canonical Neuroinflammation & Glial Markers (AD vs Old)",
       subtitle = "No significant upregulation detected in bulk cortex (TNF & IL6 filtered out by low counts)",
       x = "log2 Fold Change (AD vs Old)", y = NULL, colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top", strip.text.y = element_text(angle = 0))
ggsave("results/figures/final/fig5_immune_markers.png", p5, width = 8, height = 7, dpi = 300)

cat("Done! Generated files:\n")
print(list.files("results/figures/final", full.names = TRUE))
