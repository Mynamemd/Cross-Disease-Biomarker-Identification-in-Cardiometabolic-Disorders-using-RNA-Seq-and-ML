setwd("//wsl.localhost/Ubuntu-22.04/home/mda/Pan_Metabolics_diseases_Biomarkers")

library(ggplot2)
library(org.Hs.eg.db)

pan <- read.csv("merged_results/Final_28_with_padj.csv")
pan$gene <- gsub("\\..*", "", pan$gene)  # ENSG00000162366.1 -> ENSG00000162366

# ENSG se SYMBOL me convert
pan$symbol <- mapIds(org.Hs.eg.db, 
                     keys = pan$gene, 
                     column = "SYMBOL", 
                     keytype = "ENSEMBL", 
                     multiVals = "first")

# whose don't map remain in ENSG form
pan$symbol[is.na(pan$symbol)] <- pan$gene[is.na(pan$symbol)]
print(pan) # Now it will gene symbol name

# Highlight - BSG gene who common in both (ACS+SA
pan$type <- ifelse(pan$symbol == "BSG", "Pan (ACS+SA)", "Unique")

png("merged_results/DotPlot_Pan_SYMBOL_final.png", width=2800, height=2000, res=300)
ggplot(pan, aes(x=disease, y=reorder(symbol, disease), color=type)) + 
  geom_point(size=5) +
  scale_color_manual(values=c("Unique"="steelblue", "Pan (ACS+SA)"="red")) +
  theme_bw() +
  theme(axis.text.y = element_text(size=10, face="italic")) +
  labs(title="Distribution of Pan-Metabolic Biomarkers",
       x="Disease", y="Gene Symbol", color="") +
  geom_text(aes(label=symbol), hjust=-0.3, size=3, show.legend = FALSE)
dev.off()

# UpSet plot
library(UpSetR)
gene_list_symbol <- split(pan$symbol, pan$disease)

png("merged_results/UpSet_plot_SYMBOL.png", width=2500, height=1500, res=300)
upset(fromList(gene_list_symbol), nsets=length(gene_list_symbol), 
      order.by="freq", main.bar.color="steelblue")
dev.off()
