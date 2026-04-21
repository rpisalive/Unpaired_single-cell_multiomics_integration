library(MOFA2)
library(pheatmap)
library(RColorBrewer)
library(clue)
library(dplyr)

outdir <- "C:/Users/49152/Downloads/Multi-omics/sc_vs_scot+/MOFA_complied/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#1 Models loading
m_single <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/single_cell_trained_model_MOFAcomplied.hdf5")
m_scot   <- load_model("C:/Users/49152/Downloads/Multi-omics/SCOT_plus/MOFA_output/SCOTplus_trained_model_MOFAcomp.hdf5")

#2 Extract factor matrices
Z_single <- get_factors(m_single, factors = "all", as.data.frame = FALSE)[[1]]
Z_scot   <- get_factors(m_scot, factors = "all", as.data.frame = FALSE)[[1]]

if (is.null(rownames(Z_single)) || is.null(rownames(Z_scot))) {
  stop("Factor matrices must have sample names as rownames.")
}

#3 Aligning shared samples
common_cells <- intersect(rownames(Z_single), rownames(Z_scot))

if (length(common_cells) == 0) {
  stop("No shared sample IDs between the reference and SCOT+ MOFA models.")
}

Z_single_aligned <- Z_single[match(common_cells, rownames(Z_single)), , drop = FALSE]
Z_scot_aligned   <- Z_scot[match(common_cells, rownames(Z_scot)), , drop = FALSE]

stopifnot(identical(rownames(Z_single_aligned), rownames(Z_scot_aligned)))

colnames(Z_single_aligned) <- paste0("Reference_F", seq_len(ncol(Z_single_aligned)))
colnames(Z_scot_aligned)   <- paste0("SCOT_F", seq_len(ncol(Z_scot_aligned)))

cat("Shared samples used for concordance:", length(common_cells), "\n")
cat("Reference factors:", ncol(Z_single_aligned), "\n")
cat("SCOT+ factors:", ncol(Z_scot_aligned), "\n")

#4 Sample-level factor correlation
cor_mat <- cor(Z_single_aligned, Z_scot_aligned, use = "pairwise.complete.obs",
               method = "pearson")

# Hungarian matching
n_ref  <- nrow(cor_mat)
n_scot <- ncol(cor_mat)

if (n_ref <= n_scot) {
  cost_mat <- 1 - abs(cor_mat)
  assignment <- solve_LSAP(cost_mat)
  
  matched_pairs <- data.frame(ReferenceFactor = rownames(cor_mat),
                              SCOTFactor = colnames(cor_mat)[assignment],
                              Correlation = cor_mat[cbind(seq_len(n_ref), assignment)],
                              stringsAsFactors = FALSE)
  
  cor_mat_ordered <- cor_mat[, assignment, drop = FALSE]
  colnames(cor_mat_ordered) <- matched_pairs$SCOTFactor
} else {
  # solve on transposed matrix if SCOT has fewer factors
  cost_mat_t <- 1 - abs(t(cor_mat))
  assignment_t <- solve_LSAP(cost_mat_t)
  
  matched_pairs <- data.frame(ReferenceFactor = rownames(cor_mat)[assignment_t],
                              SCOTFactor = colnames(cor_mat),
                              Correlation = cor_mat[cbind(assignment_t, seq_len(n_scot))],
                              stringsAsFactors = FALSE)
  
  cor_mat_ordered <- cor_mat[assignment_t, , drop = FALSE]
  rownames(cor_mat_ordered) <- matched_pairs$ReferenceFactor
}

write.csv(cor_mat_ordered, file.path(outdir, "sample_level_concordance.csv"),
          row.names = TRUE)
write.csv(matched_pairs, file.path(outdir, "matched_factor_pairs.csv"), row.names = FALSE)

# Sample-level heatmap
breaks <- seq(-1, 1, length.out = 101)
colors <- colorRampPalette(c("blue", "white", "red"))(100)

svg(file.path(outdir, "Figure_SampleConcordance.svg"), width = 10, height = 8)
pheatmap(cor_mat_ordered, cluster_rows = TRUE, cluster_cols = FALSE, color = colors,
         breaks = breaks,
         main = "Sample-level concordance between reference and SCOT+ models",
         fontsize = 14, fontsize_row = 14, fontsize_col = 14, fontsize_number = 18,
         display_numbers = TRUE, angle_col = 45)
dev.off()

# Sample-level summary statistics
mean_abs_corr <- mean(abs(cor_mat), na.rm = TRUE)
cat("Mean absolute factor correlation:", round(mean_abs_corr, 3), "\n")

max_corr <- apply(abs(cor_mat), 1, max, na.rm = TRUE)
cat("Mean of max correlations per reference factor:", round(mean(max_corr), 3), "\n")

strong_matches <- sum(max_corr > 0.5)
cat("Number of reference factors with strong match (>0.5):", strong_matches, "\n")

#5 Feature-level concordance
W_single <- get_weights(m_single, views = "all", factors = "all")
W_scot   <- get_weights(m_scot, views = "all", factors = "all")

# Shared features by view
common_genes <- intersect(rownames(W_single[["Transcriptomics"]]), 
                          rownames(W_scot[["Transcriptomics"]]))

common_proteins <- intersect(rownames(W_single[["Proteomics"]]),
                             rownames(W_scot[["Proteomics"]]))

if (length(common_genes) == 0) {
  stop("No shared transcriptomic features between the two MOFA models.")
}
if (length(common_proteins) == 0) {
  stop("No shared proteomic features between the two MOFA models.")
}

# Matched factor pairs from sample-level concordance
rna_corr_signed  <- numeric(nrow(matched_pairs))
prot_corr_signed <- numeric(nrow(matched_pairs))

rna_corr_abs  <- numeric(nrow(matched_pairs))
prot_corr_abs <- numeric(nrow(matched_pairs))

for (idx in seq_len(nrow(matched_pairs))) {
  ref_name  <- matched_pairs$ReferenceFactor[idx]
  scot_name <- matched_pairs$SCOTFactor[idx]
  sign_flip <- sign(matched_pairs$Correlation[idx])
  if (sign_flip == 0) sign_flip <- 1
  
  i <- match(ref_name, colnames(Z_single_aligned))
  j <- match(scot_name, colnames(Z_scot_aligned))
  
  # Transcriptomics weights
  w1_rna <- W_single[["Transcriptomics"]][common_genes, i]
  w2_rna <- W_scot[["Transcriptomics"]][common_genes, j] * sign_flip
  
  # Proteomics weights
  w1_prot <- W_single[["Proteomics"]][common_proteins, i]
  w2_prot <- W_scot[["Proteomics"]][common_proteins, j] * sign_flip
  
  # Signed correlations after sign alignment
  rna_corr_signed[idx]  <- cor(w1_rna,  w2_rna,  use = "pairwise.complete.obs", method = "pearson")
  prot_corr_signed[idx] <- cor(w1_prot, w2_prot, use = "pairwise.complete.obs", method = "pearson")
  
  # Absolute-loading correlations (optional complementary summary)
  rna_corr_abs[idx]  <- cor(abs(w1_rna),  abs(w2_rna),  use = "pairwise.complete.obs", method = "pearson")
  prot_corr_abs[idx] <- cor(abs(w1_prot), abs(w2_prot), use = "pairwise.complete.obs", method = "pearson")
}

loading_mat_signed <- cbind(Transcriptomics = rna_corr_signed,
                            Proteomics = prot_corr_signed)
rownames(loading_mat_signed) <- matched_pairs$ReferenceFactor

loading_mat_abs <- cbind(Transcriptomics = rna_corr_abs,
                         Proteomics = prot_corr_abs)
rownames(loading_mat_abs) <- matched_pairs$ReferenceFactor

write.csv(loading_mat_signed, file.path(outdir, "feature_level_concordance_signed.csv"),
          row.names = TRUE)

write.csv(loading_mat_abs, file.path(outdir, "feature_level_concordance_absolute.csv"),
          row.names = TRUE)

# Feature-level signed loading concordance heatmap
svg(file.path(outdir, "Figure_FeatureConcordance_signed.svg"), width = 10, height = 8)
pheatmap(loading_mat_signed, cluster_rows = TRUE, cluster_cols = FALSE,
         color = colors, breaks = breaks, 
         main = "Feature-level concordance between reference and SCOT+ models",
         fontsize = 14, fontsize_row = 14, fontsize_col = 14, fontsize_number = 18,
         display_numbers = TRUE, angle_col = 0)
dev.off()

# Absolute-loading heatmap
svg(file.path(outdir, "Figure_FeatureConcordance_absolute.svg"), width = 10, height = 8)
pheatmap(loading_mat_abs, cluster_rows = TRUE, cluster_cols = FALSE,
  color = colors,
  breaks = breaks,
  main = "Feature-level absolute-loading concordance between reference and SCOT+ models",
  fontsize = 12, fontsize_row = 14, fontsize_col = 14, fontsize_number = 18,
  display_numbers = TRUE, angle_col = 0)
dev.off()

# Feature-level summary statistics
cat("Mean RNA signed loading correlation:", round(mean(rna_corr_signed, na.rm = TRUE), 3), "\n")
cat("Mean protein signed loading correlation:", round(mean(prot_corr_signed, na.rm = TRUE), 3), "\n")
cat("RNA signed correlations > 0.5:", sum(rna_corr_signed > 0.5, na.rm = TRUE), "\n")
cat("Protein signed correlations > 0.5:", sum(prot_corr_signed > 0.5, na.rm = TRUE), "\n")

cat("Mean RNA absolute-loading correlation:", round(mean(rna_corr_abs, na.rm = TRUE), 3), "\n")
cat("Mean protein absolute-loading correlation:", round(mean(prot_corr_abs, na.rm = TRUE), 3), "\n")
cat("RNA absolute-loading correlations > 0.5:", sum(rna_corr_abs > 0.5, na.rm = TRUE), "\n")
cat("Protein absolute-loading correlations > 0.5:", sum(prot_corr_abs > 0.5, na.rm = TRUE), "\n")