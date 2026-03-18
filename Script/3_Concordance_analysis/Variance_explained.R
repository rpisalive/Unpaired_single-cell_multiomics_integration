library(MOFA2)
library(dplyr)
library(ggplot2)
library(clue)

outdir <- "PATH_TO_OUTPUT_DIRECTORY"
m_single <- load_model("PATH/TO/single_cell_trained_model.hdf5")
m_scot   <- load_model("PATH/TO/SCOTplus_trained_model.hdf5")

#Extract factor matrices
Z_single <- get_factors(m_single, factors = "all", as.data.frame = FALSE)[[1]]
Z_scot   <- get_factors(m_scot, factors = "all", as.data.frame = FALSE)[[1]]

#Align cells
common_cells <- intersect(rownames(Z_single), rownames(Z_scot))
Z_single_aligned <- Z_single[common_cells, , drop = FALSE]
Z_scot_aligned   <- Z_scot[common_cells, , drop = FALSE]

#Compute factor correlation matrix
cor_mat <- cor(Z_single_aligned, Z_scot_aligned, use = "pairwise.complete.obs", method = "pearson")
cost_mat <- 1 - abs(cor_mat)
assignment <- solve_LSAP(cost_mat)  # SCOT+ factor matched to reference factor

#Extract variance explained per factor
var_single_df <- get_variance_explained(m_single, as.data.frame = TRUE)$r2_per_factor
var_scot_df   <- get_variance_explained(m_scot, as.data.frame = TRUE)$r2_per_factor

#Rename columns for clarity
var_single_df <- var_single_df %>% rename(R2_ref = value)
var_scot_df   <- var_scot_df %>% rename(R2_scot = value)

#Map SCOT+ factors to reference using Hungarian assignment
factor_map <- setNames(paste0("Factor", assignment), paste0("Factor", seq_along(assignment)))
var_scot_df$factor <- factor_map[var_scot_df$factor]

#Merge reference and SCOT+ variance explained
var_merged <- var_single_df %>%  left_join(var_scot_df, by = c("view", "factor"))

#Compute correlation per view
views <- unique(var_merged$view)
for (v in views) {
  df <- var_merged %>% filter(view == v)
  corr <- cor(df$R2_ref, df$R2_scot)
  cat("Variance explained correlation for", v, ":", round(corr, 3), "\n")
}

p <- ggplot(var_merged, aes(x = R2_ref, y = R2_scot, color = view)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~view) +
  theme_bw() +
  labs(x = "Reference variance explained", y = "SCOT+ variance explained",
       title = "Variance explained consistency between MOFA+ models")
print(p)

write.csv(var_merged, file.path(outdir, "variance_consistency.csv"), row.names = TRUE)

ggsave(filename = file.path(outdir, "Variance_consistency.svg"), plot = p, device = "svg", width = 7, height = 5)