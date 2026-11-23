suppressPackageStartupMessages({library(data.table); library(DESeq2)})

counts <- as.data.frame(fread("data/GSE159699/GSE159699_summary_count.star.txt.gz"))
rownames(counts) <- make.unique(as.character(counts[[1]])); counts <- counts[, -1]
counts <- counts[, sapply(counts, is.numeric)]
info <- data.frame(sample = colnames(counts), row.names = colnames(counts))
info$group <- sub(".*-", "", info$sample)

# 1. Sequencing depth and mitochondrial fraction
lib <- colSums(counts)
mt  <- colSums(counts[grepl("^MT-", rownames(counts)), , drop = FALSE]) / lib
cat("\n--- Library size and mito fraction ---\n")
print(data.frame(sample = names(lib), group = info$group,
                 million_counts = round(lib / 1e6, 2), mito_pct = round(100 * mt, 1), row.names = NULL))

# 2. PCA and cluster membership
dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), colData = info, design = ~ 1)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
mat <- assay(vst(dds, blind = TRUE))
pc  <- prcomp(t(mat))
info$PC1 <- pc$x[, 1]
right <- info$PC1 > 25
cat("\n--- Samples in the right cluster ---\n"); print(info$sample[right])

# 3. Genes driving PC1
ld <- pc$rotation[, 1]
cat("\n--- 20 genes HIGHER in right cluster ---\n"); print(names(sort(ld, decreasing = TRUE))[1:20])
cat("\n--- 20 genes LOWER in right cluster ---\n");  print(names(sort(ld))[1:20])

# 4. Brain cell-type markers: right cluster vs rest (log-scale values)
markers <- c(neuron = "SNAP25", neuron = "RBFOX3", neuron = "SYT1", neuron = "SLC17A7", neuron = "GAD1",
             astrocyte = "GFAP", astrocyte = "AQP4", microglia = "AIF1", microglia = "P2RY12",
             microglia = "CX3CR1", oligo = "MBP", oligo = "PLP1", oligo = "MOG", blood = "HBB", blood = "HBA1")
present <- markers[markers %in% rownames(mat)]
tab <- data.frame(celltype = names(present), gene = unname(present),
                  right_cluster = round(rowMeans(mat[present, right, drop = FALSE]), 2),
                  rest = round(rowMeans(mat[present, !right, drop = FALSE]), 2), row.names = NULL)
tab$difference <- tab$right_cluster - tab$rest
cat("\n--- Marker genes: right cluster vs rest ---\n"); print(tab)
