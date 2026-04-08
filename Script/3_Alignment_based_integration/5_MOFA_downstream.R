library(ggplot2)
library(dplyr)
library(MOFA2)
library(svglite)
library(ggokabeito)
library(extrafont)
library(readr)

setwd("C:/Users/49152/Downloads/Multi-omics/")

#1 Model loading
model <- load_model("SCOT_plus/MOFA_output/SCOTplus_trained_model_MOFAcomp.hdf5")

outdir <- "C:/Users/49152/Downloads/Multi-omics/SCOT_plus/MOFA_output/graphs/MOFA+complied_SCOTplus/Downstream_general/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

loadfonts(device = "win")
theme_set(theme_bw(base_size = 16, base_family = "Arial"))

#2 Metadata loading and alignment
# Use the aligned metadata exported after SCOT+,
aligned_metadata <- read.csv("SCOT_plus/AGW_results/aligned_metadata.csv",
                             stringsAsFactors = FALSE,check.names = FALSE)

# If the first column is an exported row-name column, rename it
if (colnames(aligned_metadata)[1] %in% c("", "X", "Unnamed: 0")) {
  colnames(aligned_metadata)[1] <- "SampleID"
}

# Model sample names
model_samples <- unname(unlist(samples_names(model)))

# Keep only metadata rows present in the model
metadata <- aligned_metadata %>%
  filter(SampleID %in% model_samples) %>%
  distinct(SampleID, .keep_all = TRUE)

# Align metadata row order to the model sample order
metadata <- metadata[match(model_samples, metadata$SampleID), , drop = FALSE]

# Safety checks
if (nrow(metadata) != length(model_samples)) {
  stop("Mismatch between model samples and aligned metadata rows after alignment.")
}
if (!all(metadata$SampleID == model_samples)) {
  stop("Sample order mismatch between model and aligned metadata.")
}

# Rename join key for MOFA samples_metadata
metadata_for_join <- metadata %>%
  rename(sample = SampleID)

model@samples_metadata <- model@samples_metadata %>%
  left_join(metadata_for_join, by = "sample")

print(head(model@samples_metadata))
print(colnames(model@samples_metadata))

#3 Definition
views <- c("Transcriptomics", "Proteomics")
K <- get_dimensions(model)$K
factors_to_plot <- seq_len(min(5, K))

#4 Data overview
p_overview <- plot_data_overview(model)

ggsave(filename = file.path(outdir, "Data_overview.svg"), plot = p_overview,
       device = svglite, width = 8, height = 6)

#5 Variance explained
if (!is.null(model@cache$variance_explained)) {
  print(model@cache$variance_explained$r2_total)
  print(model@cache$variance_explained$r2_per_factor)
} else {
  message("Variance explained cache is empty; plots will compute from model.")
}

p_var <- plot_variance_explained(model, x = "view", y = "factor", plot_total = FALSE)

ggsave(filename = file.path(outdir, "Figure2A.svg"), plot = p_var, device = svglite,
       width = 8, height = 6)

p_var_total <- plot_variance_explained(model, x = "view", y = "factor",
                                       plot_total = TRUE)

ggsave(filename = file.path(outdir, "Variance_explained_total.svg"),
       plot = p_var_total, device = svglite, width = 8, height = 6)

#6 Factor plots
p_factor <- plot_factor(model, factors = factors_to_plot, color_by = "CellLine",
                        dot_size = 1) +
  theme_bw(base_size = 12) +
  scale_colour_okabe_ito() +
  theme(axis.title.x = element_blank(), axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(), axis.title   = element_text(size = 12),
        axis.text.y  = element_text(size = 12), legend.title = element_text(size = 12),
        legend.text  = element_text(size = 12))

ggsave(filename = file.path(outdir, "Figure2B.svg"), plot = p_factor, device = svglite,
       width = 6, height = 4)

p_violin <- plot_factor(model, factors = factors_to_plot, color_by = "CellLine", 
                        dot_size = 3, dodge = TRUE, legend = FALSE, add_violin = TRUE,
                        violin_alpha = 0.25) +
  scale_color_manual(values = c("C10" = "black", "SVEC" = "red")) +
  scale_fill_manual(values = c("C10" = "black", "SVEC" = "red"))

print(p_violin)

ggsave(filename = file.path(outdir, "Factor_violin.svg"), plot = p_violin,
       device = svglite, width = 8, height = 6)

p_factors_grid <- plot_factors(model, factors = factors_to_plot,
                               color_by = "CellLine")

ggsave(filename = file.path(outdir, "Factors_grid.svg"), plot = p_factors_grid,
       device = svglite, width = 10, height = 8)

#7 Feature weights, top features, scatters and heatmaps
for (v in views) {
  for (f in seq_len(K)) {
    
    message("Processing: ", v, " | Factor ", f)
    
    # Weights plot
    p_weights <- plot_weights(model, view = v, factor = f, nfeatures = 10,
                              scale = TRUE, abs = FALSE, text_size = 2)
    
    ggsave(file.path(outdir, paste0(v, "_F", f, "_weights.svg")), p_weights,
           width = 6, height = 4, device = svglite)
    
    # Top weights plot
    p_top <- plot_top_weights(model, view = v, factor = f, nfeatures = 20)
    
    ggsave(file.path(outdir, paste0(v, "_F", f, "_topweights.svg")), p_top,
           width = 6, height = 4, device = svglite)
    
    # Extract and save weight table
    weights_df <- get_weights(model, views = v, factors = f, as.data.frame = TRUE)
    
    # Standardise column names
    if (!"feature" %in% colnames(weights_df)) {
      stop("Expected column 'feature' not found in weights table.")
    }
    if (!"value" %in% colnames(weights_df)) {
      stop("Expected column 'value' not found in weights table.")
    }
    
    weights_df <- weights_df %>%
      mutate(feature_clean = sub("_(rna|prot)$", "", feature),
             abs_value = abs(value)) %>%
      arrange(desc(abs_value))
    
    # Save full table
    write.csv(weights_df, file = file.path(outdir, paste0(v, "_F", f, "_weights_table.csv")),
              row.names = FALSE, quote = FALSE)
    
    # Save top 20 cleaned table
    top_weights_df <- weights_df %>%
      slice_head(n = 20)
    
    write.csv(top_weights_df, file = file.path(outdir, paste0(v, "_F", f, "_top20_weights_table.csv")),
              row.names = FALSE, quote = FALSE)
    
    # Scatter plot
    p_scatter <- plot_data_scatter(model, view = v, factor = f, features = 12,
                                   add_lm = TRUE, color_by = "CellLine")
    
    ggsave(file.path(outdir, paste0(v, "_F", f, "_scatter.svg")), p_scatter,
           width = 10, height = 8, device = svglite)
    
    # Heatmap
    svglite(file.path(outdir, paste0(v, "_F", f, "_heatmap.svg")), width = 10, height = 8)
    
    plot_data_heatmap(model, view = v, factor = f, features = 30, cluster_rows = TRUE,
                      cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
                      fontsize_col = 6)
    
    dev.off()
  }
}

#8 UMAP
set.seed(42)

# Use leading factors
factors_umap <- seq_len(min(3, K))
model <- run_umap(model, factors = factors_umap)

p_umap <- plot_dimred(model, method = "UMAP", color_by = "CellLine",
                      color_name = "Cell Line", shape_by = NULL, label = FALSE,
                      dot_size = 2, legend = TRUE)

p_umap <- p_umap +
  scale_colour_okabe_ito() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16),
        legend.title = element_text(size = 16), legend.text = element_text(size = 16))

ggsave(filename = file.path(outdir, "UMAP_cellline.svg"), plot = p_umap, device = svglite,
       width = 8, height = 6)