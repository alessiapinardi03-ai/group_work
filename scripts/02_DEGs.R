library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)

data_KIRC_matched <- readRDS("data/data_KIRC_matched.rds")

# 1. Take the ID for the two groups
normal_samples <- data_KIRC_matched$barcode[data_KIRC_matched$sample_type == "Solid Tissue Normal"]
tumor_samples <- data_KIRC_matched$barcode[data_KIRC_matched$sample_type == "Primary Tumor"]

# 2. DEA
# Normalization 
all_genes <- TCGAanalyze_DEA(
  mat1 = assay(data_KIRC_matched)[, normal_samples],
  mat2 = assay(data_KIRC_matched)[, tumor_samples],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  pipeline = "edgeR",
  fdr.cut = 1,               
  logFC.cut = 0           
)

# thresholds
fdr_thresh <- 1e-100      
fc_thresh  <- log2(4) 

#DEGs
DEGs_results <- all_genes[all_genes$FDR <= fdr_thresh & abs(all_genes$logFC) >= fc_thresh, ]

# 3. Verify the number of genes
n_genes <- nrow(DEGs_results)
message("Number of DEGs: ", n_genes)

# 4. Save the results (rds and csv)
saveRDS(all_genes, file = "data/all_genes_KIRC.rds")
saveRDS(DEGs_results, file = "data/DEGs_KIRC_results.rds")
write.csv(DEGs_results, file = "data/DEGs_KIRC_results.csv", row.names = TRUE)

# 4. volcano plot
volcano_data <- all_genes

# Classification of genes
volcano_data$status <- "Not Significant"
volcano_data$status[volcano_data$FDR <= fdr_thresh & volcano_data$logFC >= fc_thresh]  <- "Up-regulated"
volcano_data$status[volcano_data$FDR <= fdr_thresh & volcano_data$logFC <= -fc_thresh] <- "Down-regulated"

# to avoid infinite values
volcano_data$neg_log10_FDR <- -log10(pmax(volcano_data$FDR, 1e-300))

p_volcano <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_FDR, color = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down-regulated" = "#2b83ba",
    "Not Significant" = "grey75",
    "Up-regulated" = "#d7191c"
  )) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot - TCGA-KIRC (Tumor vs Normal)",
    x = "log2(Fold Change)",
    y = "-log10(FDR)",
    color = "Expression"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_volcano)

# save the volcano plot in png 
if (!dir.exists("plots")) dir.create("plots")
ggsave("plots/volcano_plot_KIRC.png", plot = p_volcano, width = 8, height = 6, dpi = 300)