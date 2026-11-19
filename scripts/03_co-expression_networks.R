# task 3
# computation

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
library(SummarizedExperiment)
library(igraph)

data_KIRC_matched <- readRDS("data/data_KIRC_matched.rds")
DEGs_results      <- readRDS("data/DEGs_KIRC_results.rds")

normal_samples <- data_KIRC_matched$barcode[data_KIRC_matched$sample_type == "Solid Tissue Normal"]
tumor_samples  <- data_KIRC_matched$barcode[data_KIRC_matched$sample_type == "Primary Tumor"]

# 1. Take the names of 393 DEGs
deg_genes <- rownames(DEGs_results)

# 2. Take the expression matrix just for DEGs and use log2(x+1) to avoid the falsification of Pearson correlation 
expr_normal <- log2(assay(data_KIRC_matched)[deg_genes, normal_samples] + 1)
expr_tumor  <- log2(assay(data_KIRC_matched)[deg_genes, tumor_samples] + 1)

# 3. Compute the matrix of Pearson correlation between genes
cor_normal <- cor(t(expr_normal), method = "pearson")
cor_tumor  <- cor(t(expr_tumor), method = "pearson")

# 4. Matrix of binary adiacence 
rho_thresh <- 0.7

adj_normal <- ifelse(abs(cor_normal) >= rho_thresh, 1, 0)
adj_tumor  <- ifelse(abs(cor_tumor)  >= rho_thresh, 1, 0)

# Put the diagonal axes to 0 because a gene must not be linked to itself 
diag(adj_normal) <- 0
diag(adj_tumor)  <- 0

# 5. Create the nets
net_normal <- graph_from_adjacency_matrix(adj_normal, mode = "undirected", diag = FALSE)
net_tumor  <- graph_from_adjacency_matrix(adj_tumor,  mode = "undirected", diag = FALSE)

# 6. Check the nets
message("--- Normal net (threshold ", rho_thresh, ") ---")
message("N°links: ", ecount(net_normal))
message("Linked nodes (degree > 0): ", sum(degree(net_normal) > 0), " su ", vcount(net_normal))
message("Net density: ", round(edge_density(net_normal), 4))

message("---  TUMOR NET (threshold ", rho_thresh, ") ---")
message("N°links: ", ecount(net_tumor))
message("Linked nodes (degree > 0): ", sum(degree(net_tumor) > 0), " su ", vcount(net_tumor))
message("Net density: ", round(edge_density(net_tumor), 4))

# analysis
library(ggplot2)

# 7. Compute genes' degree
deg_normal <- degree(net_normal)
deg_tumor  <- degree(net_tumor)

# 8. Make the graph of the Degree Distribution 
df_deg <- data.frame(
  Degree = c(deg_normal, deg_tumor),
  Condition = rep(c("Normal", "Tumor"), each = length(deg_normal))
)

p_deg <- ggplot(df_deg, aes(x = Degree, fill = Condition)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.7) +
  facet_wrap(~ Condition, scales = "free") +
  scale_fill_manual(values = c("Normal" = "#2b83ba", "Tumor" = "#d7191c")) +
  labs(
    title = "Degree Distribution - Co-expression Networks (TCGA-KIRC)",
    x = "Degree (k)",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

print(p_deg)
ggsave("plots/degree_distribution_task3.png", plot = p_deg, width = 9, height = 4.5, dpi = 300)

# 9. Check the Scale-Free (Power-law fit on nodes with degree > 0)
check_scale_free <- function(deg_vector, name) {
  deg_pos <- deg_vector[deg_vector > 0]
  
  # P(k): probability of observing a node of grade k 
  deg_tab <- table(deg_pos)
  k  <- as.numeric(names(deg_tab))
  pk <- as.numeric(deg_tab) / sum(deg_tab)
  
  # Linear regression log10(P(k)) ~ log10(k)
  lm_fit <- lm(log10(pk) ~ log10(k))
  r2    <- summary(lm_fit)$r.squared
  gamma <- -coef(lm_fit)[2]
  
  message("--- Scale-Free check: ", name, " ---")
  message("R^2 (fit log-log): ", round(r2, 4))
  message("Gamma exponent: ", round(gamma, 3))
}

check_scale_free(deg_normal, "Normal")
check_scale_free(deg_tumor,  "Tumor")

# 10. Find the HUB (5% of the nodes with the highest degree)
top_n <- ceiling(0.05 * vcount(net_normal)) # 5% di 393 = 20 genes

hubs_normal <- names(sort(deg_normal, decreasing = TRUE))[1:top_n]
hubs_tumor  <- names(sort(deg_tumor,  decreasing = TRUE))[1:top_n]

# 11. Compare HUB between normal and tumor
common_hubs        <- intersect(hubs_normal, hubs_tumor)
only_tumor_hubs    <- setdiff(hubs_tumor, hubs_normal)
only_normal_hubs   <- setdiff(hubs_normal, hubs_tumor)

message("\n--- COMPARISON HUB (Top ", top_n, " genes) ---")
message("Common Hub: ", length(common_hubs))
message("Tumor exclusive HUB: ", length(only_tumor_hubs))
message("Normal exclusive HUB: ", length(only_normal_hubs))

# Print the names of Tumor exclusive HUB
print("Tumor exclusive HUB:")
print(only_tumor_hubs)

# traslate the ID Ensembl of Tumor HUB in Gene Symbols 
gene_info <- as.data.frame(rowData(data_KIRC_matched))
tumor_hubs_symbols <- gene_info$gene_name[match(only_tumor_hubs, rownames(gene_info))]
print(tumor_hubs_symbols)

# 12. Save nets and HUB
saveRDS(list(
  net_normal = net_normal,
  net_tumor = net_tumor,
  hubs_normal = hubs_normal,
  hubs_tumor = hubs_tumor,
  only_tumor_hubs = only_tumor_hubs,
  only_normal_hubs = only_normal_hubs,
  tumor_hubs_symbols = tumor_hubs_symbols
), file = "data/task3_results.rds")

