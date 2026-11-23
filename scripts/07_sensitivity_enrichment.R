suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db)})
set.seed(42)

df <- read.csv("results/tables/DE_sensitivity_all_adj.csv")
df <- df[!is.na(df$padj), ]
map <- suppressWarnings(bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
map <- map[!duplicated(map$SYMBOL), ]
df$entrez <- map$ENTREZID[match(df$gene, map$SYMBOL)]
universe <- unique(na.omit(df$entrez))

# ---- ORA (GO BP) ----
run_ora <- function(genes, label) {
  genes <- unique(na.omit(genes))
  ego <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = universe,
                  pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
  out <- as.data.frame(ego)
  cat("\n=====", label, ": genes =", length(genes), " enriched terms =", nrow(out), "=====\n")
  if (nrow(out) > 0) {
    write.csv(out, paste0("results/tables/GO_ORA_", label, ".csv"), row.names = FALSE)
    print(head(out[, c("Description", "Count", "p.adjust")], 6), row.names = FALSE)
  }
}
run_ora(df$entrez[df$status == "Up in AD"],   "sens_up_in_AD")
run_ora(df$entrez[df$status == "Down in AD"], "sens_down_in_AD")

# ---- GSEA (GO BP) ----
rk <- df[!is.na(df$entrez) & !is.na(df$stat), ]
rk <- rk[order(-abs(rk$stat)), ]; rk <- rk[!duplicated(rk$entrez), ]
glist <- sort(setNames(rk$stat, rk$entrez), decreasing = TRUE)
cat("\nRunning sensitivity GSEA (1 to 5 minutes)...\n")
g <- gseGO(glist, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
           minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0, verbose = FALSE)
g <- setReadable(g, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gd <- as.data.frame(g)
write.csv(gd, "results/tables/GSEA_GO_BP_sens.csv", row.names = FALSE)

cat("\n===== Sensitivity GSEA: significant terms =", nrow(gd), "=====\n")
cat("\nTop ACTIVATED:\n");  print(head(gd[gd$NES > 0, c("Description", "NES", "p.adjust")], 8), row.names = FALSE)
cat("\nTop SUPPRESSED:\n"); print(head(gd[gd$NES < 0, c("Description", "NES", "p.adjust")], 8), row.names = FALSE)

# ---- Compare with primary GSEA ----
prim <- read.csv("results/tables/GSEA_GO_BP.csv")
m <- merge(prim[, c("ID", "NES")], gd[, c("ID", "NES")], by = "ID", suffixes = c("_prim", "_sens"))
cat("\n===== GSEA Overlap: Primary =", nrow(prim), " Sensitivity =", nrow(gd), " Shared =", nrow(m), "=====\n")
if (nrow(m) > 2) cat("NES correlation on shared terms:", round(cor(m$NES_prim, m$NES_sens), 3),
                     "| Same direction:", sum(sign(m$NES_prim) == sign(m$NES_sens)), "of", nrow(m), "\n")

# ---- Immune GSEA terms & key immune genes ----
imm <- gd[grepl("immun|inflamm|cytokine|complement|microgli|interferon|leukocyte|macrophage|phagocyt",
                gd$Description, ignore.case = TRUE), ]
cat("\n===== Sensitivity immune GSEA terms:", nrow(imm), "=====\n")
if (nrow(imm) > 0) print(head(imm[order(imm$p.adjust), c("Description", "NES", "p.adjust")], 15), row.names = FALSE)

genes <- c("TREM2","TYROBP","C1QA","C1QB","C1QC","C3","CD68","AIF1","CX3CR1","P2RY12",
           "GFAP","CHI3L1","IL1B","TNF","IL6","HLA-DRA","HLA-DRB1","B2M","CD74","APOE","CLU")
sub <- df[df$gene %in% genes, c("gene", "log2FC", "padj")]
cat("\n===== Key immune genes (sensitivity, padj < 0.25) =====\n")
print(sub[order(sub$padj), ][sub[order(sub$padj), "padj"] < 0.25, ], row.names = FALSE)
