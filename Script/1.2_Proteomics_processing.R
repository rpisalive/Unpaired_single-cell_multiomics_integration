# =========================================================
# Package versions
# =========================================================
# R version: 4.4.3 (2025-02-28 ucrt)
# rstudioapi 0.17.1
# tidyverse 2.0.0 
# dplyr 1.1.4 
# tidyr 1.3.1 
# tibble 3.3.0 
# readr 2.1.5 
# org.Mm.eg.db 3.19.1 
# AnnotationDbi 1.66.0 
# limma 3.60.6 

library(tidyverse)
library(readr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(limma)

#1 Raw data loading
raw <- read_tsv("PATH/TO/C10SVEC_singlecells_Protein_intensities.tsv",
                show_col_types = FALSE)

# Remove contaminants
raw <- raw[!grepl("contam_sp", raw$PROTID), ]

sample_cols <- setdiff(colnames(raw), "PROTID")

# Replace zeros with NA
raw[sample_cols] <- lapply(raw[sample_cols], function(x) {
  x <- as.numeric(x)
  x[x == 0] <- NA
  x
})

expr_mat <- as.matrix(raw[, sample_cols])
rownames(expr_mat) <- raw$PROTID
storage.mode(expr_mat) <- "numeric"

cat("Initial protein rows:", nrow(expr_mat), "\n")
cat("Initial samples:", ncol(expr_mat), "\n")

#2 Log transformation & Sample-wise median centering
log_norm <- log2(expr_mat)

# Median-centre each sample
sample_medians <- apply(log_norm, 2, median, na.rm = TRUE)
log_norm <- sweep(log_norm, 2, sample_medians, FUN = "-")

#3 Protein <-> gene mapping
log_norm_df <- as.data.frame(log_norm) %>%
  rownames_to_column(var = "PROTID")

log_norm_long <- log_norm_df %>%
  pivot_longer(cols = -PROTID, names_to = "SampleID", values_to = "Intensity") %>%
  dplyr::filter(!is.na(Intensity))

uniprot_ids <- unique(log_norm_long$PROTID)

geneSymbols_db <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = uniprot_ids,
  columns = c("SYMBOL", "UNIPROT"),
  keytype = "UNIPROT") %>%
  as_tibble() %>%
  dplyr::rename(Gene_db = SYMBOL)

# Supplementary mapping table from Nanosplit
convert <- read_tsv("PATH/TO/convert_uniprot.tsv", show_col_types = FALSE)

if (!all(c("UNIPROT", "Gene") %in% colnames(convert))) {
  stop("convert_uniprot.tsv must contain columns named 'UNIPROT' and 'Gene'")
}

convert <- convert %>%
  dplyr::rename(Gene_conv = Gene)

# Prioritise conversion-table gene names if available
gene_map <- geneSymbols_db %>%
  left_join(convert, by = "UNIPROT") %>%
  mutate(Gene = coalesce(Gene_conv, Gene_db)) %>%
  dplyr::select(UNIPROT, Gene) %>%
  dplyr::filter(!is.na(Gene)) %>%
  distinct()

# Attach gene names to measured proteins
log_norm_long <- log_norm_long %>%
  mutate(UNIPROT = PROTID) %>%
  left_join(gene_map, by = "UNIPROT") %>%
  dplyr::filter(!is.na(Gene))

# If multiple proteins map to same gene in same sample, keep the most intense one
log_norm_long_gene <- log_norm_long %>%
  group_by(SampleID, Gene) %>%
  slice_max(order_by = Intensity, n = 1, with_ties = FALSE) %>%
  ungroup()

# Convert back to gene x sample matrix
expr_gene <- log_norm_long_gene %>%
  dplyr::select(SampleID, Gene, Intensity) %>%
  pivot_wider(names_from = SampleID, values_from = Intensity) %>%
  column_to_rownames(var = "Gene") %>%
  as.matrix()

storage.mode(expr_gene) <- "numeric"

cat("Mapped gene rows:", nrow(expr_gene), "\n")
cat("Samples after mapping:", ncol(expr_gene), "\n")

if (any(duplicated(rownames(expr_gene)))) {
  rownames(expr_gene) <- make.unique(rownames(expr_gene))
}

#4 Prepare metadata based on cell names
sample_ids <- colnames(expr_gene)

coldata <- data.frame(SampleID = sample_ids, Chip = factor(substr(sample_ids, 1, 3)),
                      CellLine = factor(ifelse(grepl("SVEC", sample_ids), "SVEC", "C10"),
                                        levels = c("C10", "SVEC")), stringsAsFactors = FALSE)

rownames(coldata) <- coldata$SampleID

cat("Cell-line composition:\n")
print(table(coldata$CellLine))

cat("Chip composition:\n")
print(table(coldata$Chip))

cat("Chip x CellLine table:\n")
tab <- table(coldata$Chip, coldata$CellLine)
print(tab)

#5 Cell QC

detected_per_sample <- colSums(!is.na(expr_gene))
coldata$DetectedProteins <- detected_per_sample[coldata$SampleID]

keep_samples <- detected_per_sample >= 800

expr_gene <- expr_gene[, keep_samples, drop = FALSE]
coldata <- coldata[colnames(expr_gene), , drop = FALSE]

cat("Samples retained after proteomics QC:", ncol(expr_gene), "\n")
print(table(coldata$CellLine))

#6 Feature filtering
# Within-group presence thresholds
expr_c10 <- expr_gene[, coldata$CellLine == "C10", drop = FALSE]
expr_svec <- expr_gene[, coldata$CellLine == "SVEC", drop = FALSE]

present_c10 <- rowMeans(!is.na(expr_c10))
present_svec <- rowMeans(!is.na(expr_svec))

# Keep features in at least 50% of either cell line
keep_genes <- (present_c10 >= 0.5) | (present_svec >= 0.5)

expr_gene <- expr_gene[keep_genes, , drop = FALSE]

cat("Proteins/genes retained after missingness filter:", nrow(expr_gene), "\n")

#7 Batch effect correction (modify when necessary, not applicable in this study)
design_test <- model.matrix(~ CellLine + Chip, data = coldata)

if (qr(design_test)$rank < ncol(design_test)) {
  message("Chip and CellLine are confounded. Skipping batch correction.")
  expr_final <- expr_gene
} else {
  expr_final <- removeBatchEffect(
    x = expr_gene,
    batch = coldata$Chip,
    design = model.matrix(~ CellLine, data = coldata)
  )
}

#8 Vatiable feature selection
# Chip and CellLine are fully confounded, batch-corrected variance modelling not possible.
# Rank features directly by variance in the final processed matrix
feature_var <- apply(expr_final, 1, var, na.rm = TRUE)
feature_var <- feature_var[!is.na(feature_var) & feature_var > 0]

n_proteins <- min(1500, length(feature_var))
top_features <- names(sort(feature_var, decreasing = TRUE))[seq_len(n_proteins)]

prot_matrix <- expr_final[top_features, , drop = FALSE]

cat("Number of variable proteomic features selected:", nrow(prot_matrix), "\n")

if (any(duplicated(rownames(prot_matrix)))) {
  rownames(prot_matrix) <- make.unique(rownames(prot_matrix))
}

#9 Export
prot_df <- data.frame(Gene = rownames(prot_matrix), prot_matrix, check.names = FALSE)
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
processed_dir <- file.path(script_dir, "processed_data")
write.csv(prot_df, file.path(processed_dir, "processed_proteomics.csv"), row.names = FALSE,
          quote = FALSE)

write.csv(data.frame(coldata, check.names = FALSE),
          file.path(processed_dir, "processed_proteomics_metadata.csv"), row.names = FALSE,
          quote = FALSE)

#10 PCA diagnostics

complete_features <- rowSums(is.na(prot_matrix)) == 0
cat("Complete-case features used for PCA:", sum(complete_features), "\n")

if (sum(complete_features) >= 2) {
  pca <- prcomp(t(prot_matrix[complete_features, , drop = FALSE]), scale. = FALSE)
  
  plot(pca$x[, 1], pca$x[, 2], col = ifelse(coldata$CellLine == "C10", "red", "blue"),
       pch = 16, xlab = "PC1", ylab = "PC2", main = "Proteomics PCA after preprocessing")
  legend("topright", legend = c("C10", "SVEC"), col = c("red", "blue"), pch = 16)
  
  missing_per_cell <- colSums(is.na(prot_matrix))
  cat("Correlation PC1 vs missingness:", round(cor(pca$x[, 1], missing_per_cell), 3), "\n")
  if (ncol(pca$x) >= 2) {
    cat("Correlation PC2 vs missingness:", round(cor(pca$x[, 2], missing_per_cell), 3), "\n")
  }
  
  plot(pca$x[, 1], pca$x[, 2], col = coldata$Chip, pch = 16, xlab = "PC1", ylab = "PC2",
       main = "Proteomics PCA coloured by Chip (confounded with CellLine)")
  legend("topright", legend = levels(coldata$Chip), col = seq_along(levels(coldata$Chip)), pch = 16)
} else {
  cat("Too few complete-case proteomic features for PCA.\n")
}