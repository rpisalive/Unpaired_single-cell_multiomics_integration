library(SingleCellExperiment)
library(scran)
library(scater)
library(dplyr)
library(readr)

setwd("PATH_TO_WORKING_DIRECTORY")

#Read counts
raw_counts <- read_tsv("PATH/TO/C10SVEC_singlecells_Counts.txt", show_col_types = FALSE)
colnames(raw_counts)[1] <- "Gene"

count_matrix <- as.matrix(raw_counts[, -1])
rownames(count_matrix) <- raw_counts$Gene

sce <- SingleCellExperiment(assays = list(counts = count_matrix))

#Metadata parsing
cell_names <- colnames(sce)
cell_info <- strsplit(cell_names, "_")

chip <- sapply(cell_info, function(x) x[1])
cell_line <- sapply(cell_info, function(x) x[2])

#Special case: "07J_SVEC_C10"
cell_line[cell_names == "07J_SVEC_C10"] <- "SVEC"

colData(sce)$Chip <- chip
colData(sce)$CellLine <- cell_line

#Filtering
keep_genes <- rowSums(counts(sce) > 0) > 3
sce <- sce[keep_genes, ]

#Count number of detected genes per cell
detected_genes <- apply(counts(sce), 2, function(x) sum(x >= 5))

#Add to colData
colData(sce)$DetectedGenes5 <- detected_genes

#Apply thresholds by CellLine
keep_cells <- (colData(sce)$CellLine == "C10"  & colData(sce)$DetectedGenes5 >= 3000) |
  (colData(sce)$CellLine == "SVEC" & colData(sce)$DetectedGenes5 >= 2000)

sce <- sce[, keep_cells]

#Library size normalization + log-transform
sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)

#HVG selection
dec <- modelGeneVar(sce)
top_hvgs <- getTopHVGs(dec, n = 2000)
sce <- sce[top_hvgs, ]

#Export
rna_matrix <- as.matrix(logcounts(sce))

rna_df <- data.frame(rna_matrix, check.names = FALSE)

write.csv(rna_df,
          file = "PATH/TO/transcriptomics_for_SCOT.csv",
          row.names = TRUE, quote = FALSE)

metadata <- as.data.frame(colData(sce))

metadata$SampleID <- rownames(metadata)
metadata <- metadata[, c("SampleID", setdiff(names(metadata), "SampleID"))]

write.csv(metadata,
          file = "PATH/TO/transcriptomics_metadata.csv",
          row.names = FALSE, quote = FALSE)

#PCA
pca <- prcomp(t(rna_matrix), scale. = FALSE)
plot(pca$x[,1], pca$x[,2],
     col = ifelse(metadata$CellLine == "C10", "red", "blue"),
     pch = 16,
     xlab = "PC1", ylab = "PC2",
     main = "RNA PCA: C10 (red) vs SVEC (blue)")
