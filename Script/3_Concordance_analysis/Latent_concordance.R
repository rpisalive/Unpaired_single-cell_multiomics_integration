library(MOFA2)
library(pheatmap)
library(RColorBrewer)
library(clue)

outdir <- "PATH_TO_OUTPUT_DIRECTORY"

#Load models
m_single <- load_model("PATH/TO/single_cell_trained_model.hdf5")
m_scot   <- load_model("PATH/TO/SCOTplus_trained_model.hdf5")

#Extract factor matrices
Z_single <- get_factors(m_single, factors = "all", as.data.frame = FALSE)[[1]]
Z_scot   <- get_factors(m_scot, factors = "all", as.data.frame = FALSE)[[1]]

#Align shared cells
common_cells <- intersect(rownames(Z_single), rownames(Z_scot))
Z_single_aligned <- Z_single[common_cells, , drop = FALSE]
Z_scot_aligned   <- Z_scot[common_cells, , drop = FALSE]
colnames(Z_single_aligned) <- paste0("Single_F", seq_len(ncol(Z_single_aligned)))
colnames(Z_scot_aligned)   <- paste0("SCOT_F", seq_len(ncol(Z_scot_aligned)))

#Compute correlation matrix
cor_mat <- cor(Z_single_aligned, Z_scot_aligned, use = "pairwise.complete.obs", method = "pearson")

#Reorder SCOT+ factors for best match using Hungarian algorithm
cost_mat <- 1 - abs(cor_mat)
assignment <- solve_LSAP(cost_mat)
cor_mat_ordered <- cor_mat[, assignment]
write.csv(cor_mat_ordered, file.path(outdir, "sample_level_concordance.csv"), row.names = TRUE)

#Heatmap colors
breaks <- seq(-1, 1, length.out = 101)
colors <- colorRampPalette(c("blue", "white", "red"))(100)

#Sample-level factor concordance SVG
svg(file.path(outdir, "Figure_SampleConcordance.svg"), width = 10, height = 8)
pheatmap(cor_mat_ordered, cluster_rows = TRUE, cluster_cols = FALSE,
         color = colors, breaks = breaks,
         main = "Sample-level concordance of between reference and SCOT+ models",
         fontsize = 14, fontsize_row = 14, fontsize_col = 14, fontsize_number = 18,
         display_numbers = TRUE, angle_col = 45)

dev.off()

#Summary statistics for sample-level concordance
mean_abs_corr <- mean(abs(cor_mat), na.rm = TRUE)
cat("Mean absolute factor correlation:", round(mean_abs_corr, 3), "\n")
max_corr <- apply(cor_mat, 1, function(x) max(abs(x), na.rm = TRUE))
cat("Mean of max correlations per reference factor:", round(mean(max_corr), 3), "\n")
strong_matches <- sum(max_corr > 0.5)
cat("Number of reference factors with strong match (>0.5):", strong_matches, "\n")

#Feature level correlation
W_single <- get_weights(m_single, views = "all", factors = "all")
W_scot   <- get_weights(m_scot, views = "all", factors = "all")

#Align shared features
common_genes <- intersect(rownames(W_single[["Transcriptomics"]]), rownames(W_scot[["Transcriptomics"]]))
common_proteins <- intersect(rownames(W_single[["Proteomics"]]), rownames(W_scot[["Proteomics"]]))
rna_corr  <- numeric(length(assignment))
prot_corr <- numeric(length(assignment))

for (i in seq_along(assignment)) {
  j <- assignment[i]
  
  #RNA
  w1_rna <- W_single[["Transcriptomics"]][common_genes, i]
  w2_rna <- W_scot[["Transcriptomics"]][common_genes, j]
  rna_corr[i] <- cor(abs(w1_rna), abs(w2_rna), use = "pairwise.complete.obs")
  
  #Protein
  w1_prot <- W_single[["Proteomics"]][common_proteins, i]
  w2_prot <- W_scot[["Proteomics"]][common_proteins, j]
  prot_corr[i] <- cor(abs(w1_prot), abs(w2_prot), use = "pairwise.complete.obs")
}

loading_mat <- cbind(Transcriptomics = rna_corr, Proteomics = prot_corr)
rownames(loading_mat) <- paste0("F", seq_along(assignment))
write.csv(loading_mat, file.path(outdir, "feature_level_concordance.csv"), row.names = TRUE)

#Feature-level loading concordance SVG
svg(file.path(outdir, "Figure_FeatureConcordance.svg"), width = 10, height = 8)

pheatmap(loading_mat, cluster_rows = TRUE, cluster_cols = FALSE, color = colors,
         breaks = breaks,
         main = "Feature-level concordance of between reference and SCOT+ models",
         fontsize = 14, fontsize_row = 14, fontsize_col = 14, fontsize_number = 18,
         display_numbers = TRUE, angle_col = 0)

dev.off()

#Summary statistics for loading correlations
cat("Mean RNA loading correlation:", round(mean(rna_corr), 3), "\n")
cat("Mean protein loading correlation:", round(mean(prot_corr), 3), "\n")
cat("RNA factors > 0.5:", sum(rna_corr > 0.5), "\n")
cat("Protein factors > 0.5:", sum(prot_corr > 0.5), "\n")