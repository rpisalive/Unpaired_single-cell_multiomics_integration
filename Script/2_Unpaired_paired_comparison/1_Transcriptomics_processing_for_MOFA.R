library(SingleCellExperiment)
library(scran)
library(scater)
library(limma)
library(dplyr)
library(readr)

setwd("C:/Users/49152/Downloads/Multi-omics/")

#1 Raw data loading
raw_counts <- read_tsv("Github_data/Transcriptomics/C10SVEC_singlecells_Counts.txt", show_col_types = FALSE)
colnames(raw_counts)[1] <- "Gene"
count_matrix <- as.matrix(raw_counts[, -1])
rownames(count_matrix) <- raw_counts$Gene
sce <- SingleCellExperiment(assays = list(counts = count_matrix))

#2 Prepare metadata based on cell names
cell_names <- colnames(sce)
cell_info <- strsplit(cell_names, "_")

chip <- sapply(cell_info, function(x) x[1])
cell_line <- sapply(cell_info, function(x) x[2])

# Special case
cell_line[cell_names == "07J_SVEC_C10"] <- "SVEC"

colData(sce)$SampleID <- cell_names
colData(sce)$Chip <- factor(chip)
colData(sce)$CellLine <- factor(cell_line, levels = c("C10", "SVEC"))

# Sanity check
table(colData(sce)$CellLine, useNA = "ifany")
table(colData(sce)$Chip, useNA = "ifany")

#3 Cell QC

detected_genes_5 <- colSums(counts(sce) >= 5)
library_size <- colSums(counts(sce))
colData(sce)$DetectedGenes5 <- detected_genes_5
colData(sce)$LibrarySize <- library_size

# Cell line–specific thresholds
keep_cells <- ((colData(sce)$CellLine == "C10"  & detected_genes_5 >= 3000) |
                 (colData(sce)$CellLine == "SVEC" & detected_genes_5 >= 2000))

sce <- sce[, keep_cells]

cat("Cells retained after QC:", ncol(sce), "\n")
print(table(colData(sce)$CellLine))

#4 Gene filtering

keep_genes <- rowSums(counts(sce) > 1) >= 10
sce <- sce[keep_genes, ]

cat("Genes retained after filtering:", nrow(sce), "\n")

#5 Normalisation

sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)

# Sanity check
summary(sizeFactors(sce))

#6 Batch effect correction (modify when necessary, not applicable in this study)

tab <- table(colData(sce)$Chip, colData(sce)$CellLine)
print(tab)

design_test <- model.matrix(~ CellLine + Chip, data = as.data.frame(colData(sce)))

if (qr(design_test)$rank < ncol(design_test)) {
  message("Chip and CellLine are confounded after filtering. Skipping batch correction.")
  assay(sce, "logcounts") <- logcounts(sce)
} else {
  assay(sce, "logcounts") <- removeBatchEffect(
    x = logcounts(sce),
    batch = colData(sce)$Chip,
    design = model.matrix(~ CellLine, data = as.data.frame(colData(sce)))
  )
}

#7 HVG Selection
# Chip and CellLine are fully confounded, batch-corrected variance modelling not possible.
# Rank genes by variance in the processed log-expression matrix, results are reflecting combined CellLine/Chip structure.

log_mat <- assay(sce, "logcounts")

gene_var <- apply(log_mat, 1, var, na.rm = TRUE)
gene_var <- gene_var[!is.na(gene_var) & gene_var > 0]

n_hvgs <- min(1500, length(gene_var))
top_hvgs <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_hvgs)]

sce_hvg <- sce[top_hvgs, ]

cat("Number of HVGs selected:", nrow(sce_hvg), "\n")

#8 Matrix preparation for MOFA

rna_matrix <- as.matrix(assay(sce_hvg, "logcounts"))

if (any(duplicated(rownames(rna_matrix)))) {
  rownames(rna_matrix) <- make.unique(rownames(rna_matrix))
}

#9 Export matrix and metadata
rna_df <- data.frame(Gene = rownames(rna_matrix), rna_matrix, check.names = FALSE)

write.csv(rna_df, file = "MOFA/input/4_Unpaired_sc/transcriptomics_for_MOFA.csv",
          row.names = FALSE, quote = FALSE)

metadata <- as.data.frame(colData(sce_hvg))
metadata$SampleID <- colnames(sce_hvg)
metadata <- metadata[, c("SampleID", setdiff(colnames(metadata), "SampleID"))]

write.csv(metadata, file = "MOFA/input/4_Unpaired_sc/transcriptomics_metadata.csv",
          row.names = FALSE, quote = FALSE)

#10 PCA Diagnostics
pca <- prcomp(t(rna_matrix), scale. = FALSE)

plot(pca$x[, 1], pca$x[, 2], col = ifelse(metadata$CellLine == "C10", "red", "blue"),
     pch = 16, xlab = "PC1", ylab = "PC2", main = "RNA PCA after preprocessing")
legend("topright", legend = c("C10", "SVEC"), col = c("red", "blue"), pch = 16)

lib_size_final <- metadata$LibrarySize
detected_final <- metadata$DetectedGenes5

cor_pc1_lib <- cor(pca$x[, 1], lib_size_final)
cor_pc1_detected <- cor(pca$x[, 1], detected_final)

cat("Correlation PC1 vs library size:", round(cor_pc1_lib, 3), "\n")
cat("Correlation PC1 vs detected genes:", round(cor_pc1_detected, 3), "\n")

plot(pca$x[, 1], pca$x[, 2], col = metadata$Chip, pch = 16, xlab = "PC1", ylab = "PC2",
     main = "RNA PCA coloured by Chip (confounded with CellLine)")
legend("topright", legend = levels(metadata$Chip), col = seq_along(levels(metadata$Chip)), pch = 16)