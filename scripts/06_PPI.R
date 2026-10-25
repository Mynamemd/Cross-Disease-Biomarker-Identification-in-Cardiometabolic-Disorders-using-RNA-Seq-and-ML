# PPI Network
library(org.Hs.eg.db)
setwd("//wsl.localhost/Ubuntu-22.04/home/mda/Pan_Metabolics_diseases_Biomarkers")

pan <- read.csv("merged_results/pan_all_biomarkers.csv")
pan$gene <- gsub("\\..*", "", pan$gene)
pan$symbol <- mapIds(org.Hs.eg.db, keys=pan$gene, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
pan$symbol[is.na(pan$symbol)] <- pan$gene[is.na(pan$symbol)]

pan_clean <- pan[!grepl("^ENSG", pan$symbol), ]

# Saving file
write.csv(pan_clean, "merged_results/pan_clean_SYMBOL.csv", row.names=FALSE)
print(pan_clean)

# Now PPi will run
library(STRINGdb)
string_db <- STRINGdb$new(version="12.0", species=9606, score_threshold=400, input_directory="")
mapped <- string_db$map(pan_clean, "symbol", removeUnmappedRows=TRUE)

png("merged_results/PPI_Pan_final.png", width=3000, height=3000, res=300)
string_db$plot_network(mapped$STRING_id)
dev.off()
