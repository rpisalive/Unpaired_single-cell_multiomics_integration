library(MOFA2)
library(dplyr)
library(ggplot2)
library(clue)
library(readr)
library(tidyr)

# Output path
outdir <- "C:/Users/49152/Downloads/Multi-omics/sc_vs_scot+/MOFA_complied/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#1 Models loading
m_single <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/single_cell_trained_model_MOFAcomplied.hdf5")
m_scot   <- load_model("C:/Users/49152/Downloads/Multi-omics/SCOT_plus/MOFA_output/SCOTplus_trained_model_MOFAcomp.hdf5")

#2 Helper function for factor names standardisation
standardize_factor_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  # Already in correct format
  is_factor <- grepl("^Factor[0-9]+$", x)
  x[is_factor] <- x[is_factor]
  
  # Reference_F1 or SCOT_F1 -> Factor1
  idx <- grepl("^(Reference_F|SCOT_F)[0-9]+$", x)
  x[idx] <- sub("^(Reference_F|SCOT_F)([0-9]+)$", "Factor\\2", x[idx])
  
  # F1 -> Factor1
  idx <- grepl("^F[0-9]+$", x)
  x[idx] <- sub("^F([0-9]+)$", "Factor\\1", x[idx])
  
  # factor1 -> Factor1
  idx <- grepl("^factor[0-9]+$", x, ignore.case = TRUE)
  x[idx] <- sub("^factor([0-9]+)$", "Factor\\1", x[idx], ignore.case = TRUE)
  
  # Bare numeric -> FactorN
  idx <- grepl("^[0-9]+$", x)
  x[idx] <- paste0("Factor", x[idx])
  
  x
}

#3 Load or recompute matched factor pairs
matched_pairs_file <- file.path(outdir, "matched_factor_pairs.csv")

if (file.exists(matched_pairs_file)) {
  matched_pairs <- read.csv(matched_pairs_file, stringsAsFactors = FALSE)
  cat("Loaded matched factor pairs from existing file.\n")
} else {
  cat("matched_factor_pairs.csv not found. Recomputing factor matching from latent factors.\n")
  
  Z_single <- get_factors(m_single, factors = "all", as.data.frame = FALSE)[[1]]
  Z_scot   <- get_factors(m_scot, factors = "all", as.data.frame = FALSE)[[1]]
  
  common_cells <- intersect(rownames(Z_single), rownames(Z_scot))
  if (length(common_cells) == 0) {
    stop("No shared sample IDs between the reference and SCOT+ MOFA models.")
  }
  
  Z_single_aligned <- Z_single[match(common_cells, rownames(Z_single)), , drop = FALSE]
  Z_scot_aligned   <- Z_scot[match(common_cells, rownames(Z_scot)), , drop = FALSE]
  
  stopifnot(identical(rownames(Z_single_aligned), rownames(Z_scot_aligned)))
  
  colnames(Z_single_aligned) <- paste0("Reference_F", seq_len(ncol(Z_single_aligned)))
  colnames(Z_scot_aligned)   <- paste0("SCOT_F", seq_len(ncol(Z_scot_aligned)))
  
  cor_mat <- cor(Z_single_aligned, z_scot_aligned, use = "pairwise.complete.obs",
                 method = "pearson")
  cost_mat <- 1 - abs(cor_mat)
  assignment <- solve_LSAP(cost_mat)
  
  matched_pairs <- data.frame(ReferenceFactor = rownames(cor_mat),
                              SCOTFactor = colnames(cor_mat)[assignment],
                              Correlation = cor_mat[cbind(seq_len(nrow(cor_mat)), assignment)],
                              stringsAsFactors = FALSE)
  
  write.csv(matched_pairs, matched_pairs_file, row.names = FALSE)
}

# Standardise matched-pair column names
if (all(c("ReferenceFactor", "SCOTFactor", "Correlation") %in% colnames(matched_pairs))) {
  matched_pairs <- matched_pairs %>%
    rename(factor_ref = ReferenceFactor, factor_scot = SCOTFactor,
           factor_correlation = Correlation)
} else if (!all(c("factor_ref", "factor_scot", "factor_correlation") %in% colnames(matched_pairs))) {
  stop("matched_pairs file does not contain expected columns.")
}

matched_pairs <- matched_pairs %>%
  mutate(factor_ref = standardize_factor_name(factor_ref),
         factor_scot = standardize_factor_name(factor_scot))

print(matched_pairs)

#4 Extarct variance explained per factor
var_single_obj <- get_variance_explained(m_single, as.data.frame = TRUE)
var_scot_obj   <- get_variance_explained(m_scot, as.data.frame = TRUE)

if (is.null(var_single_obj$r2_per_factor) || is.null(var_scot_obj$r2_per_factor)) {
  stop("Could not extract r2_per_factor from one or both models.")
}

var_single_df <- var_single_obj$r2_per_factor
var_scot_df   <- var_scot_obj$r2_per_factor

cat("\nColumns in reference variance explained table:\n")
print(colnames(var_single_df))

cat("\nColumns in SCOT+ variance explained table:\n")
print(colnames(var_scot_df))

required_cols <- c("view", "factor", "value")

if (!all(required_cols %in% colnames(var_single_df))) {
  stop("Reference variance explained table does not contain expected columns: view, factor, value.")
}
if (!all(required_cols %in% colnames(var_scot_df))) {
  stop("SCOT+ variance explained table does not contain expected columns: view, factor, value.")
}

var_single_df <- var_single_df %>%
  rename(R2_ref = value, factor_ref = factor) %>%
  mutate(factor_ref = standardize_factor_name(factor_ref)) %>%
  select(view, factor_ref, R2_ref)

var_scot_df <- var_scot_df %>%
  rename(R2_scot = value, factor_scot = factor) %>%
  mutate(factor_scot = standardize_factor_name(factor_scot)) %>%
  select(view, factor_scot, R2_scot)

cat("\nUnique factor names in matched_pairs$factor_scot:\n")
print(unique(matched_pairs$factor_scot))

cat("\nUnique factor names in var_scot_df$factor_scot:\n")
print(unique(var_scot_df$factor_scot))

#5 Map SCOT factors to reference factors using matched pairs
var_scot_mapped <- var_scot_df %>%
  left_join(matched_pairs, by = "factor_scot")

if (all(is.na(var_scot_mapped$factor_ref))) {
  stop("Join failed: factor_ref is all NA after mapping SCOT factors to matched pairs.")
}

# Merge reference and SCOT+ variance explained by reference factor + view
var_merged <- var_single_df %>%
  left_join(var_scot_mapped %>%
      select(view, factor_ref, factor_scot, factor_correlation, R2_scot),
      by = c("view", "factor_ref")) %>%
  mutate(R2_difference = R2_scot - R2_ref, abs_R2_difference = abs(R2_difference))

write.csv(var_merged, file.path(outdir, "variance_consistency.csv"), row.names = FALSE)
print(var_merged)

#6 Correlation of variance explained per view
views <- unique(var_merged$view)

for (v in views) {
  df <- var_merged %>% filter(view == v)
  
  if (nrow(df) >= 2) {
    corr <- cor(df$R2_ref, df$R2_scot, use = "pairwise.complete.obs", method = "pearson")
    cat("Variance explained correlation for", v, ":", round(corr, 3), "\n")
  } else {
    cat("Variance explained correlation for", v, ": not enough matched factors\n")
  }
  
  cat("Mean absolute R2 difference for", v, ":", round(mean(df$abs_R2_difference, na.rm = TRUE), 3), "\n")
}

#7 Scatter plot
p <- ggplot(var_merged, aes(x = R2_ref, y = R2_scot, color = view, label = factor_ref)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text(nudge_y = 0.5, size = 3, show.legend = FALSE) +
  facet_wrap(~view, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Reference variance explained (%)", y = "SCOT+ variance explained (%)",
       title = "Variance explained concordance between MOFA+ models")

print(p)

ggsave(filename = file.path(outdir, "Variance_consistency.svg"), plot = p, device = "svg",
       width = 7, height = 5)

#8 Bar plot for direct comparison
plot_df <- var_merged %>%
  select(view, factor_ref, factor_scot, R2_ref, R2_scot) %>%
  pivot_longer(cols = c(R2_ref, R2_scot), names_to = "Model", values_to = "R2") %>%
  mutate(Model = recode(Model, R2_ref = "Reference", R2_scot = "SCOT+"),
         FactorPair = paste0(factor_ref, "\n", factor_scot))

p_bar <- ggplot(plot_df, aes(x = FactorPair, y = R2, fill = Model)) +
  geom_col(position = "dodge") +
  facet_wrap(~view, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(x = "Matched factor pair", y = "Variance explained (%)",
       title = "Per-factor variance explained comparison")

print(p_bar)

ggsave(filename = file.path(outdir, "Variance_consistency_bar.svg"), plot = p_bar,
       device = "svg", width = 8, height = 5)