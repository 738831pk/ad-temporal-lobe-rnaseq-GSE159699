suppressPackageStartupMessages(library(dplyr))
g <- read.csv("results/tables/GSEA_immune_terms.csv")
cat("--- Leading-edge genes per immune term ---\n")
for (i in seq_len(nrow(g))) cat(g$Description[i], ":", g$core_enrichment[i], "\n\n")

d <- read.csv("results/tables/DE_primary_excl2.csv")
genes <- c("TREM2","TYROBP","C1QA","C1QB","C1QC","C3","CD68","AIF1","CX3CR1","P2RY12",
           "GFAP","CHI3L1","IL1B","TNF","IL6","HLA-DRA","HLA-DRB1","B2M","CD74","APOE","CLU")
cat("--- Key microglia/complement/inflammation genes ---\n")
print(d[d$gene %in% genes, c("gene","log2FC","padj")] %>% arrange(padj), row.names = FALSE)
