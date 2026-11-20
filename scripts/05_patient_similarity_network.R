# task 5
# a
if (!requireNamespace("SNFtool", quietly = TRUE)) install.packages("SNFtool")
library(SummarizedExperiment)
library(SNFtool)
library(igraph)

data_KIRC_matched <- readRDS("data/data_KIRC_matched.rds")
DEGs_results      <- readRDS("data/DEGs_KIRC_results.rds")

# 1. Expression matrix Tumor just for 393 DEGs
tumor_idx     <- which(data_KIRC_matched$sample_type == "Primary Tumor")
tumor_samples <- data_KIRC_matched$barcode[tumor_idx]
deg_genes     <- rownames(DEGs_results)

# Transformation in log2(x + 1) + row = patient | columns = genes
expr_tumor_genes <- log2(assay(data_KIRC_matched)[deg_genes, tumor_samples] + 1)
expr_tumor_pat   <- t(expr_tumor_genes)

# Save patients' ID
patient_ids <- substr(rownames(expr_tumor_pat), 1, 12)

# Rimuovi duplicates for the same patient
keep_unique    <- !duplicated(patient_ids)
expr_tumor_pat <- expr_tumor_pat[keep_unique, ]
rownames(expr_tumor_pat) <- patient_ids[keep_unique]

message("N° of tumor patients in the PSN: ", nrow(expr_tumor_pat))

# 2. Patient Similarity Network 
# Standardization per gene (mean 0, SD 1) + #####distanza euclidea quadratica#####
expr_norm <- standardNormalization(expr_tumor_pat)
dist_expr <- dist2(as.matrix(expr_norm), as.matrix(expr_norm))

# Parameters SNF: K = local near (around N/5), alpha/sigma = variance of kernel
K     <- 15
alpha <- 0.5
W_expr <- affinityMatrix(dist_expr, K = K, sigma = alpha)
rownames(W_expr) <- colnames(W_expr) <- rownames(expr_tumor_pat)

# Create the weighted net
W_expr_nod <- W_expr
diag(W_expr_nod) <- 0
psn_expr <- graph_from_adjacency_matrix(W_expr_nod, mode = "undirected", weighted = TRUE, diag = FALSE)

# b
# 3. Community detection with Louvain algorithm
set.seed(42)
comm_expr <- cluster_louvain(psn_expr, weights = E(psn_expr)$weight)

message("\n--- COMMUNITY DETECTION Results (PSN Gene Expression)---")
message("Modularity: ", round(modularity(comm_expr), 4))
message("Community dimension (number of patient per cluster):")
print(sizes(comm_expr))

# (we show just the 15% of the strongest link)
edge_thresh <- quantile(E(psn_expr)$weight, 0.85)
psn_expr_vis <- delete_edges(psn_expr, E(psn_expr)[weight < edge_thresh])

# color distinction
comm_vec  <- membership(comm_expr)
n_comm    <- length(comm_expr)
comm_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")[1:n_comm]

V(psn_expr_vis)$color <- comm_cols[comm_vec]
# last 4 digits of the ID for a clear view
V(psn_expr_vis)$label <- substr(V(psn_expr_vis)$name, 9, 12)

# groups to highlight community areas in the graph 
comm_groups <- split(V(psn_expr_vis)$name, comm_vec)
edge_ends <- ends(psn_expr_vis, E(psn_expr_vis), names = TRUE)
is_cross  <- comm_vec[edge_ends[, 1]] != comm_vec[edge_ends[, 2]]
layout_weights <- ifelse(is_cross, 0.05, 2.5)

png("plots/psn_expr_communities_task5b.png", width = 2600, height = 2200, res = 300)
set.seed(42)
layout_expr <- layout_with_fr(psn_expr_vis, weights = E(psn_expr_vis)$weight)

plot(
  psn_expr_vis,
  layout = layout_expr,
  mark.groups = comm_groups,
  mark.col = adjustcolor(comm_cols, alpha.f = 0.15),
  mark.border = comm_cols,
  vertex.size = 11,
  vertex.frame.color = "white",
  vertex.label.cex = 0.50,
  vertex.label.color = "white",
  vertex.label.font = 2,
  edge.color = adjustcolor("grey45", alpha.f = 0.45),
  edge.width = E(psn_expr_vis)$weight / max(E(psn_expr_vis)$weight) * 3.5,
  main = paste0("Patient Similarity Network - Gene Expression (Louvain, Q = ", round(modularity(comm_expr), 3), ")")
)

legend(
  "bottomleft",
  legend = paste0("Community ", seq_len(n_comm), " (n = ", as.numeric(sizes(comm_expr)), ")"),
  pch = 21,
  pt.bg = comm_cols,
  col = "white",
  pt.cex = 2,
  bty = "n",
  cex = 0.95
)
dev.off()

# c
library(TCGAbiolinks)

# 4. Download data of somatic mutation (MAF) just for the 72 patients
if (!file.exists("data/maf_KIRC_matched.rds")) {
  query_maf <- GDCquery(
    project = "TCGA-KIRC",
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
    barcode = rownames(expr_tumor_pat)
  )
  GDCdownload(query_maf)
  maf_data <- GDCprepare(query_maf)
  saveRDS(maf_data, file = "data/maf_KIRC_matched.rds")
} else {
  maf_data <- readRDS("data/maf_KIRC_matched.rds")
}

# 5. Create the mutatinal binary matrix (patients x genes: 1 = mutate, 0 = wild-type)
# Exclusion of synonyme mutations ("Silent") which do not alter the protein 
maf_nonsilent <- maf_data[maf_data$Variant_Classification != "Silent", ]
maf_patients  <- substr(maf_nonsilent$Tumor_Sample_Barcode, 1, 12)

# Cross the patients with gene experession and somatic mutation 
common_patients <- intersect(rownames(expr_tumor_pat), unique(maf_patients))
message("Patients with both Gene-Expression and Somatic Mutations: ", length(common_patients), " su ", nrow(expr_tumor_pat))

# Filter the MAF table with just common patients for the 0/1 matrix
keep_maf <- maf_patients %in% common_patients
mut_tab  <- table(maf_patients[keep_maf], maf_nonsilent$Hugo_Symbol[keep_maf])
mut_mat  <- ifelse(mut_tab > 0, 1, 0)
mut_mat  <- mut_mat[common_patients, , drop = FALSE]

message("Tot N° of mutated genes (non-silent) in the cohort: ", ncol(mut_mat))

# 6. Build 2 ##matrici di affinità## on the same set of patients (common_patients)
# Layer 1: Gene expression 
expr_common    <- expr_tumor_pat[common_patients, , drop = FALSE]
expr_norm_comm <- standardNormalization(expr_common)
dist_expr_comm <- dist2(as.matrix(expr_norm_comm), as.matrix(expr_norm_comm))
W_expr_comm    <- affinityMatrix(dist_expr_comm, K = K, sigma = alpha)

# Layer 2: Mutational profile (binary distance / Jaccard)
dist_mut <- as.matrix(dist(mut_mat, method = "binary"))
W_mut    <- affinityMatrix(dist_mut, K = K, sigma = alpha)

# 7. Similarity Network Fusion (SNF) on 2 layer
T_iter <- 20  # n° of iterations of diffusion (standard value SNF)
W_snf  <- SNF(list(W_expr_comm, W_mut), K = K, t = T_iter)
rownames(W_snf) <- colnames(W_snf) <- common_patients

# 8. Community detection (Louvain) on the fused net SNF
W_snf_nod <- W_snf
diag(W_snf_nod) <- 0
psn_snf <- graph_from_adjacency_matrix(W_snf_nod, mode = "undirected", weighted = TRUE, diag = FALSE)

set.seed(42)
comm_snf <- cluster_louvain(psn_snf, weights = E(psn_snf)$weight)

message("\n---  COMMUNITY DETECTION Results (Fused Net SNF) ---")
message("Community number SNF: ", length(comm_snf))
message("Modularity SNF: ", round(modularity(comm_snf), 4))
message("Community dimension SNF:")
print(sizes(comm_snf))

# 9. Comparison between communities (task a/b vs. c) 
comm_expr_common <- membership(comm_expr)[common_patients]
comm_snf_vec     <- membership(comm_snf)[common_patients]

message("\n--- Community comparison' (Expression vs SNF) ---")
print(table(Expression = comm_expr_common, SNF = comm_snf_vec))

# Quantitative concordance metrics between clustering 
nmi_val <- compare(comm_expr_common, comm_snf_vec, method = "nmi")
ari_val <- compare(comm_expr_common, comm_snf_vec, method = "adjusted.rand")
message("Normalized Mutual Information (NMI): ", round(nmi_val, 4))
message("Adjusted Rand Index (ARI): ", round(ari_val, 4))

# 10. Graph of fusion net SNF 
edge_thresh_snf <- quantile(E(psn_snf)$weight, 0.85)
psn_snf_vis     <- delete_edges(psn_snf, E(psn_snf)[weight < edge_thresh_snf])

comm_snf_vec2 <- membership(comm_snf)
n_comm_snf    <- length(comm_snf)
snf_cols      <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")[1:n_comm_snf]

V(psn_snf_vis)$color <- snf_cols[comm_snf_vec2]
V(psn_snf_vis)$label <- substr(V(psn_snf_vis)$name, 9, 12)
snf_groups           <- split(V(psn_snf_vis)$name, comm_snf_vec2)

# weights to separe groups
edge_ends_snf <- ends(psn_snf_vis, E(psn_snf_vis), names = TRUE)
is_cross_snf  <- comm_snf_vec2[edge_ends_snf[, 1]] != comm_snf_vec2[edge_ends_snf[, 2]]
layout_w_snf  <- ifelse(is_cross_snf, 0.05, 2.5)

png("plots/psn_snf_communities_task5c.png", width = 2600, height = 2200, res = 300)
set.seed(42)
layout_snf <- layout_with_fr(psn_snf_vis, weights = layout_w_snf)

plot(
  psn_snf_vis,
  layout = layout_snf,
  mark.groups = snf_groups,
  mark.col = adjustcolor(snf_cols, alpha.f = 0.15),
  mark.border = snf_cols,
  vertex.size = 8,
  vertex.frame.color = "white",
  vertex.label.cex = 0.42,
  vertex.label.color = "white",
  vertex.label.font = 2,
  edge.color = ifelse(is_cross_snf, adjustcolor("grey65", alpha.f = 0.20), adjustcolor("grey25", alpha.f = 0.65)),
  edge.width = E(psn_snf_vis)$weight / max(E(psn_snf_vis)$weight) * 3.0,
  main = paste0("Fused PSN (SNF: Expression + Mutations, Louvain Q = ", round(modularity(comm_snf), 3), ")")
)

legend(
  "topleft",
  legend = paste0("SNF Comm. ", seq_len(n_comm_snf), " (n = ", as.numeric(sizes(comm_snf)), ")"),
  pch = 21,
  pt.bg = snf_cols,
  col = "white",
  pt.cex = 2,
  bty = "n",
  cex = 0.90
)
dev.off()

# Save the results
task5_results <- list(
  W_expr = W_expr,
  comm_expr = comm_expr,
  W_expr_comm = W_expr_comm,
  W_mut = W_mut,
  W_snf = W_snf,
  comm_snf = comm_snf,
  common_patients = common_patients,
  comparison_table = table(Expression = comm_expr_common, SNF = comm_snf_vec),
  nmi = nmi_val,
  ari = ari_val
)
saveRDS(task5_results, file = "data/task5_results.rds")
