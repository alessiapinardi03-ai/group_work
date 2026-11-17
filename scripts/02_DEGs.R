library(TCGAbiolinks)
library(SummarizedExperiment)

data_KIRC_matched <- readRDS("data/data_KIRC_matched.rds")

# 1. Take the ID for the two groups
normal_samples <- data_KIRC_matched$barcode[data_KIRC_matched$sample_type == "Solid Tissue Normal"]
tumor_samples <- data_KIRC_matched$barcode[data_KIRC_matched$sample_type == "Primary Tumor"]

# 2. DEA
# Normalization and DEGs
DEGs_results <- TCGAanalyze_DEA(
  mat1 = assay(data_KIRC_matched)[, normal_samples],
  mat2 = assay(data_KIRC_matched)[, tumor_samples],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  pipeline = "edgeR",
  fdr.cut = 0.05,               # P-value corretto per confronti multipli <= 0.05
  logFC.cut = log2(15)           # Soglia FC >= 4 convertita in log2FC
)

# 3. Verify the number of genes
n_genes <- nrow(DEGs_results)
message("Number of DEGs: ", n_genes)

# 4. Save the results (rds and csv)
saveRDS(DEGs_results, file = "data/DEGs_KIRC_results.rds")
write.csv(DEGs_results, file = "data/DEGs_KIRC_results.csv", row.names = TRUE)