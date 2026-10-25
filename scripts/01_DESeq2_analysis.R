# ================== STEP 0: SETUP ==================
setwd("//wsl.localhost/Ubuntu-22.04/home/mda/T2D_Biomarker_discovery/counts")

counts_raw <- read.table("gene_counts.txt", header=TRUE, row.names=1, comment.char="#")
counts_data <- counts_raw[, 6:ncol(counts_raw), drop=FALSE]
colnames(counts_data) <- gsub(".*(SRR[0-9]+).*", "\\1", colnames(counts_data))

n_samples <- ncol(counts_data)
print(paste("Total samples found:", n_samples))

counts_data <- counts_data[, order(colnames(counts_data))]

samples <- colnames(counts_data)
condition_map <- c(
  "SRR31131588" = "Control",
  "SRR31131589" = "Control",
  "SRR31131590" = "Control",
  "SRR31131592" = "Treated",
  "SRR31131595" = "Treated",
  "SRR31131596" = "Treated"
)
condition <- factor(condition_map[samples], levels = c("Control", "Treated"))
colData <- data.frame(condition = condition, row.names = samples)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "Treated", "Control"))

# ================== STEP 1: ADDING PCA FOR QC ==================
library(ggplot2)
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA - QC Check: Control vs Treated")

print(pca_plot)
ggsave("PCA_QC_plot.png", pca_plot, width=8, height=6, dpi=300)

# ================== STEP 2: ANNOTATION ==================
library(org.Hs.eg.db)
res_df <- as.data.frame(res)
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res_df),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

PADJ_CUTOFF <- 0.07
LFC_CUTOFF <- 1.0

write.csv(res_df, "DESeq2_full_results.csv")

sig_genes <- subset(res_df, padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF)
print(paste("Total significant genes:", nrow(sig_genes)))

write.csv(assay(vsd), "vst_counts_for_ML.csv")
