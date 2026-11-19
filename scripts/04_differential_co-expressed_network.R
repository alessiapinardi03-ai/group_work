# task 4
# computation

library(SummarizedExperiment)
library(igraph)

data_KIRC_matched <- readRDS("data/data_KIRC_matched.rds")
DEGs_results      <- readRDS("data/DEGs_KIRC_results.rds")
task3_results     <- readRDS("data/task3_results.rds")

normal_samples <- data_KIRC_matched$barcode[data_KIRC_matched$sample_type == "Solid Tissue Normal"]
tumor_samples  <- data_KIRC_matched$barcode[data_KIRC_matched$sample_type == "Primary Tumor"]
deg_genes      <- rownames(DEGs_results)

# 1. Expression matrix log2(x + 1) e Pearson's correlation
expr_normal <- log2(assay(data_KIRC_matched)[deg_genes, normal_samples] + 1)
expr_tumor  <- log2(assay(data_KIRC_matched)[deg_genes, tumor_samples] + 1)

cor_normal <- cor(t(expr_normal), method = "pearson")
cor_tumor  <- cor(t(expr_tumor), method = "pearson")

diag(cor_normal) <- 0
diag(cor_tumor)  <- 0

# 2. Fisher's Z transformation and computation of Z-score matrix (Tumor vs Normal)
#### atanh(r) = 0.5 * log((1 + r) / (1 - r))
z_normal <- atanh(cor_normal)
z_tumor  <- atanh(cor_tumor)

n_normal <- length(normal_samples)
n_tumor  <- length(tumor_samples)

Z_mat <- (z_tumor - z_normal) / sqrt(1 / (n_tumor - 3) + 1 / (n_normal - 3))
diag(Z_mat) <- 0

# 3. Binary Adiacence matrix with threshold |Z| >= 3 (a_ij = 0 if |Z| < 3)
z_thresh <- 6
adj_diff <- ifelse(abs(Z_mat) >= z_thresh, 1, 0)
diag(adj_diff) <- 0

# 4. Create the complete differential net
net_diff <- graph_from_adjacency_matrix(adj_diff, mode = "undirected", diag = FALSE)

# 5. Check the statistics
n_pos_edges <- sum(Z_mat[upper.tri(Z_mat)] >= z_thresh)
n_neg_edges <- sum(Z_mat[upper.tri(Z_mat)] <= -z_thresh)

message("--- Differential Net (soglia |Z| >= ", z_thresh, ") ---")
message("Tot N°links: ", ecount(net_diff))
message("  -> POSITIVE links (Z >= ", z_thresh, ", increase in Tumor): ", n_pos_edges)
message("  -> NEGATIVE links (Z <= -", z_thresh, ", decrease in Tumor): ", n_neg_edges)
message("Linked nodes (degree > 0): ", sum(degree(net_diff) > 0), " su ", vcount(net_diff))
message("Net density: ", round(edge_density(net_diff), 4))

#analysis
library(ggplot2)

# Gene Symbol Table 
gene_info <- as.data.frame(rowData(data_KIRC_matched))
get_symbols <- function(ens_ids) {
  syms <- gene_info$gene_name[match(ens_ids, rownames(gene_info))]
  ifelse(is.na(syms) | syms == "", ens_ids, syms)
}

# 6. Degree of the differential net and distribution graph
deg_diff <- degree(net_diff)

df_deg_diff <- data.frame(Degree = deg_diff)
p_deg_diff <- ggplot(df_deg_diff, aes(x = Degree)) +
  geom_histogram(bins = 30, fill = "#7570b3", color = "black", alpha = 0.8) +
  labs(
    title = paste0("Degree Distribution - Differential Co-expression Network (|Z| >= ", z_thresh, ")"),
    x = "Degree (k)",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 12)

print(p_deg_diff)
ggsave("plots/degree_distribution_task4.png", plot = p_deg_diff, width = 7, height = 4.5, dpi = 300)

# 7. Hub of differential net (Top 5% = 20 genes) + comparison with task 3
top_n <- ceiling(0.05 * vcount(net_diff))
hubs_diff <- names(sort(deg_diff, decreasing = TRUE))[1:top_n]
hubs_diff_symbols <- get_symbols(hubs_diff)

overlap_tumor_hubs  <- intersect(hubs_diff, task3_results$hubs_tumor)
overlap_normal_hubs <- intersect(hubs_diff, task3_results$hubs_normal)

message("\n--- HUB Differential Net (Top ", top_n, ") ---")
print(hubs_diff_symbols)
message("Differential Hub common with Hub Tumor (Task 3): ", length(overlap_tumor_hubs))
if (length(overlap_tumor_hubs) > 0) print(get_symbols(overlap_tumor_hubs))
message("Differential Hub common with Hub Normal (Task 3): ", length(overlap_normal_hubs))
if (length(overlap_normal_hubs) > 0) print(get_symbols(overlap_normal_hubs))

# 8. Take the subnets: Positive (Z >= 6) and Negative (Z <= -6)
adj_pos <- ifelse(Z_mat >= z_thresh, 1, 0)
adj_neg <- ifelse(Z_mat <= -z_thresh, 1, 0)
diag(adj_pos) <- 0
diag(adj_neg) <- 0

net_diff_pos <- graph_from_adjacency_matrix(adj_pos, mode = "undirected", diag = FALSE)
net_diff_neg <- graph_from_adjacency_matrix(adj_neg, mode = "undirected", diag = FALSE)

deg_pos_net <- degree(net_diff_pos)
deg_neg_net <- degree(net_diff_neg)

hubs_pos <- names(sort(deg_pos_net, decreasing = TRUE))[1:top_n]
hubs_neg <- names(sort(deg_neg_net, decreasing = TRUE))[1:top_n]

message("\n--- POSITIVE subnet (New correlation in Tumor, Z >= ", z_thresh, ") ---")
message("Positive Links: ", ecount(net_diff_pos), " | Linked nodes: ", sum(deg_pos_net > 0))
print("Top Positive Hub :")
print(get_symbols(hubs_pos))

message("\n---  NEGATIVE subnet (Lost of correlation in Tumor, Z <= -", z_thresh, ") ---")
message("Negative Links: ", ecount(net_diff_neg), " | Linked nodes: ", sum(deg_neg_net > 0))
print("Top Negative Hub :")
print(get_symbols(hubs_neg))

message("Common Hub between subnets: ", length(intersect(hubs_pos, hubs_neg)))

# 9. Plot of the subnet around the nodes with the highest degree 
top_hub_id  <- hubs_diff[1]
top_hub_sym <- get_symbols(top_hub_id)
message("\n Absolute Top Hub of differential net: ", top_hub_sym, " (Degree = ", deg_diff[top_hub_id], ")")


##rivedere definizioni dalle slide ############################################


# Estrai il vicinato diretto (ego network di ordine 1) del Top Hub
ego_nodes   <- ego(net_diff, order = 1, nodes = top_hub_id)[[1]]
sub_top_hub <- induced_subgraph(net_diff, ego_nodes)

# Color the links: Red = Z positive (gain), Blue = Z negativo (loss)
edge_ends <- ends(sub_top_hub, es = E(sub_top_hub))
edge_z    <- mapply(function(u, v) Z_mat[u, v], edge_ends[, 1], edge_ends[, 2])
E(sub_top_hub)$color <- ifelse(edge_z > 0, "#d7191c", "#2b83ba")

# Names and colors of nodes (with top hub in the center)
V(sub_top_hub)$label <- get_symbols(V(sub_top_hub)$name)
V(sub_top_hub)$color <- ifelse(V(sub_top_hub)$name == top_hub_id, "#ffae19", "#e0e0e0")
V(sub_top_hub)$size  <- ifelse(V(sub_top_hub)$name == top_hub_id, 18, 9)

# Save the figure
####################### verificare grafico se leggibile #######################
png("plots/top_hub_subnetwork_task4.png", width = 3000, height = 2800, res = 300)

center_idx  <- which(V(sub_top_hub)$name == top_hub_id)
star_layout <- layout_as_star(sub_top_hub, center = center_idx)

V(sub_top_hub)$size <- ifelse(V(sub_top_hub)$name == top_hub_id, 15, 4.5)

# Disegna la rete senza etichette (rescale = FALSE mantiene le coordinate esatte [-1, 1])
plot(
  sub_top_hub,
  layout = star_layout,
  rescale = FALSE,
  xlim = c(-1.45, 1.45),
  ylim = c(-1.45, 1.45),
  vertex.label = NA,
  vertex.frame.color = "grey40",
  edge.width = 1.2,
  main = paste0("Differential Co-expression Subnetwork of Top Hub: ", top_hub_sym)
)

# Aggiungi le etichette ruotate a raggiera verso l'esterno (zero sovrapposizioni)
for (i in seq_len(vcount(sub_top_hub))) {
  if (i == center_idx) {
    text(0, 0, labels = V(sub_top_hub)$label[i], cex = 0.7, font = 2, col = "black")
  } else {
    x <- star_layout[i, 1]
    y <- star_layout[i, 2]
    deg <- atan2(y, x) * 180 / pi
    
    if (x >= 0) {
      text(x * 1.04, y * 1.04, labels = V(sub_top_hub)$label[i],
           srt = deg, adj = c(0, 0.5), cex = 0.55, col = "black")
    } else {
      text(x * 1.04, y * 1.04, labels = V(sub_top_hub)$label[i],
           srt = deg + 180, adj = c(1, 0.5), cex = 0.55, col = "black")
    }
  }
}

legend(
  "bottomleft",
  legend = c("Gain of correlation (Z >= 6)", "Loss of correlation (Z <= -6)", paste0("Top Hub (", top_hub_sym, ")")),
  col = c("#d7191c", "#2b83ba", "#ffae19"),
  lty = c(1, 1, NA),
  pch = c(NA, NA, 19),
  lwd = 2,
  pt.cex = 1.6,
  bty = "n",
  cex = 0.85
)

dev.off()
# 10. Save the results
saveRDS(list(
  z_thresh = z_thresh,
  Z_mat = Z_mat,
  net_diff = net_diff,
  net_diff_pos = net_diff_pos,
  net_diff_neg = net_diff_neg,
  hubs_diff = hubs_diff,
  hubs_diff_symbols = hubs_diff_symbols,
  hubs_pos = hubs_pos,
  hubs_neg = hubs_neg
), file = "data/task4_results.rds")