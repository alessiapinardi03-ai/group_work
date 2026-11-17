# 1. Install packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

# 2. Set query TCGA-KIRC (Tumor and normal tissues)
query_KIRC <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# 3. Download 
GDCdownload(query_KIRC, files.per.chunk = 50)
data_KIRC <- GDCprepare(query_KIRC)

# 4. Select patients with both type of tissue
# Take ID patient from 12 fist number of the barcode
patients <- substr(data_KIRC$barcode, 1, 12)

# Create two lists
tumor_patients <- patients[data_KIRC$sample_type == "Primary Tumor"]
normal_patients <- patients[data_KIRC$sample_type == "Solid Tissue Normal"]

# Find the intersection
patients_with_both <- intersect(tumor_patients, normal_patients)

# Keep just the requested patients
data_KIRC_matched <- data_KIRC[, patients %in% patients_with_both]

# 5. Save the data in "data" folder
saveRDS(data_KIRC_matched, file = "data/data_KIRC_matched.rds")