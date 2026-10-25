setwd("//wsl.localhost/Ubuntu-22.04/home/mda/Pan_Metabolics_diseases_Biomarkers")

files <- list.files("DESeq2_results", pattern="*.csv", full.names=TRUE)
all <- do.call(rbind, lapply(files, function(f){
  df <- read.csv(f)
  gcol <- grep("gene|ENSG|X", colnames(df), value=TRUE)[1]
  if(is.na(gcol)) df$gene <- rownames(df) else df$gene <- df[[gcol]]
  df$gene <- gsub("\\..*", "", df$gene)
  return(df)
}))
final_table <- merge(pan, all[,c("gene","log2FoldChange","pvalue","padj")], by="gene", all.x=TRUE)
write.csv(final_table, "merged_results/Final_34_with_padj.csv", row.names=FALSE)
