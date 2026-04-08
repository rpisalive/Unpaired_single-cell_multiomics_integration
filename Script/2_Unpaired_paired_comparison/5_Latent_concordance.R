library(MOFA2)
library(dplyr)
library(pheatmap)
library(clue)

outdir <- "C:/Users/49152/Downloads/Multi-omics/MOFA/output/graphs/Unpaired_single_cell/Minimal_baseline/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#1 Model loading
m_ref <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/single_cell_trained_model_MOFAcomplied.hdf5")
m_unp <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/Unpaired_single_cell_model.hdf5")

#2 Factor matrices
Z_ref <- get_factors(m_ref, factors = "all", as.data.frame = FALSE)[[1]]
Z_unp <- get_factors(m_unp, factors = "all", as.data.frame = FALSE)[[1]]

Z_unp_rna  <- Z_unp[grepl("_rna$", rownames(Z_unp)), , drop = FALSE]
Z_unp_prot <- Z_unp[grepl("_prot$", rownames(Z_unp)), , drop = FALSE]

rownames(Z_unp_rna)  <- sub("_rna$", "", rownames(Z_unp_rna))
rownames(Z_unp_prot) <- sub("_prot$", "", rownames(Z_unp_prot))

colnames(Z_ref)      <- paste0("Reference_F", seq_len(ncol(Z_ref)))
colnames(Z_unp_rna)  <- paste0("Unpaired_F", seq_len(ncol(Z_unp_rna)))
colnames(Z_unp_prot) <- paste0("Unpaired_F", seq_len(ncol(Z_unp_prot)))

#3 Align shared samples
common_rna  <- intersect(rownames(Z_ref), rownames(Z_unp_rna))
common_prot <- intersect(rownames(Z_ref), rownames(Z_unp_prot))

Z_ref_rna <- Z_ref[match(common_rna, rownames(Z_ref)), , drop = FALSE]
Z_unp_rna <- Z_unp_rna[match(common_rna, rownames(Z_unp_rna)), , drop = FALSE]

Z_ref_prot <- Z_ref[match(common_prot, rownames(Z_ref)), , drop = FALSE]
Z_unp_prot <- Z_unp_prot[match(common_prot, rownames(Z_unp_prot)), , drop = FALSE]

#4 Modality-specific correlation matrices
cor_mat_rna <- cor(Z_ref_rna, Z_unp_rna, use = "pairwise.complete.obs", method = "pearson")
cor_mat_prot <- cor(Z_ref_prot, Z_unp_prot, use = "pairwise.complete.obs", method = "pearson")

combined_score <- (abs(cor_mat_rna) + abs(cor_mat_prot)) / 2

cost_mat <- 1 - combined_score
assignment <- solve_LSAP(cost_mat)

matched_pairs <- data.frame(ReferenceFactor = rownames(combined_score),
                            UnpairedFactor = colnames(combined_score)[assignment],
                            MeanAbsCorrelation = combined_score[cbind(seq_len(nrow(combined_score)), assignment)],
                            RNA_Correlation = cor_mat_rna[cbind(seq_len(nrow(cor_mat_rna)), assignment)],
                            Proteomics_Correlation = cor_mat_prot[cbind(seq_len(nrow(cor_mat_prot)), assignment)],
                            stringsAsFactors = FALSE)

write.csv(matched_pairs, file.path(outdir, "matched_factor_pairs_reference_vs_unpaired.csv"),
          row.names = FALSE)

#5 Heatmaps
svg(file.path(outdir, "LatentConcordance_RNA.svg"), width = 8, height = 6)
pheatmap(cor_mat_rna[, assignment, drop = FALSE], cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, main = "Reference vs unpaired latent concordance (RNA)")
dev.off()

svg(file.path(outdir, "LatentConcordance_Proteomics.svg"), width = 8, height = 6)
pheatmap(cor_mat_prot[, assignment, drop = FALSE], cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, main = "Reference vs unpaired latent concordance (Proteomics)")
dev.off()

svg(file.path(outdir, "LatentConcordance_meanAbs.svg"), width = 8, height = 6)
pheatmap(combined_score[, assignment, drop = FALSE], cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, main = "Reference vs unpaired latent concordance (mean absolute)")
dev.off()

#6 Summary
cat("Mean matched factor mean absolute correlation:",
    round(mean(matched_pairs$MeanAbsCorrelation, na.rm = TRUE), 3), "\n")
cat("Matched factors with mean abs correlation > 0.5:",
    sum(matched_pairs$MeanAbsCorrelation > 0.5, na.rm = TRUE), "\n")