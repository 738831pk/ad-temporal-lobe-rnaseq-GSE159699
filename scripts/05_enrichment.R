suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot); library(ggplot2)})
set.seed(42)

df <- read.csv("results/tables/DE_primary_excl2.csv")
df <- df[!is.na(df$padj), ]
map <- suppressWarnings(bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
map <- map[!duplicated(map$SYMBOL), ]
df$entrez <- map$ENTREZID[match(df$gene, map$SYMBOL)]
cat("Genes mapped to Entrez:", sum(!is.na(df$entrez)), "of", nrow(df), "\n")
universe <- unique(na.omit(df$entrez))

# ---- Over-representation analysis (GO Biological Process) ----
run_ora <- function(genes, label) {
  genes <- unique(na.omit(genes))
  ego <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = universe,
                  pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
  out <- as.data.frame(ego)
  cat("\n=====", label, ": genes =", length(genes), " enriched terms =", nrow(out), "=====\n")
  if (nrow(out) > 0) {
    write.csv(out, paste0("results/tables/GO_ORA_", label, ".csv"), row.names = FALSE)
    print(head(out[, c("Description", "Count", "p.adjust")], 8), row.names = FALSE)
    ggsave(paste0("results/figures/GO_ORA_", label, ".png"),
           dotplot(ego, showCategory = 15) + ggtitle(paste("GO BP:", label)), width = 9, height = 7, dpi = 300)
  }
}
run_ora(df$entrez[df$status == "Up in AD"],   "up_in_AD")
run_ora(df$entrez[df$status == "Down in AD"], "down_in_AD")

# ---- GSEA on all genes ranked by the DESeq2 statistic ----
rk <- df[!is.na(df$entrez) & !is.na(df$stat), ]
rk <- rk[order(-abs(rk$stat)), ]; rk <- rk[!duplicated(rk$entrez), ]
glist <- sort(setNames(rk$stat, rk$entrez), decreasing = TRUE)
cat("\nRunning GSEA (takes 1 to 5 minutes)...\n")
g <- gseGO(glist, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
           minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0, verbose = FALSE)
g <- setReadable(g, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gd <- as.data.frame(g)
write.csv(gd, "results/tables/GSEA_GO_BP.csv", row.names = FALSE)
cat("\n===== GSEA: significant terms =", nrow(gd), "=====\n")
cat("\nTop terms ACTIVATED in AD:\n")
print(head(gd[gd$NES > 0, c("Description", "NES", "p.adjust")], 8), row.names = FALSE)
cat("\nTop terms SUPPRESSED in AD:\n")
print(head(gd[gd$NES < 0, c("Description", "NES", "p.adjust")], 8), row.names = FALSE)
if (nrow(gd) > 0)
  ggsave("results/figures/GSEA_dotplot.png", dotplot(g, showCategory = 12, split = ".sign") + facet_grid(. ~ .sign),
         width = 11, height = 8, dpi = 300)

# ---- Immune / inflammation focus ----
imm <- gd[grepl("immun|inflamm|cytokine|complement|microgli|interferon|leukocyte|macrophage|phagocyt",
                gd$Description, ignore.case = TRUE), ]
write.csv(imm, "results/tables/GSEA_immune_terms.csv", row.names = FALSE)
cat("\n===== Immune-related GSEA terms:", nrow(imm), "=====\n")
if (nrow(imm) > 0) print(head(imm[order(imm$p.adjust), c("Description", "NES", "p.adjust")], 15), row.names = FALSE)
