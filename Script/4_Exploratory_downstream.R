library(ggplot2)
library(dplyr)
library(MOFA2)
library(svglite)
library(ggokabeito)
library(extrafont)
library(readr)

# Formatting and path settings
loadfonts(device = "win")
theme_set(theme_bw(base_size = 16, base_family = "Arial"))

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
model_dir <- file.path(script_dir, "trained_model")
meta_dir <- file.path(script_dir, "processed_data")
align_dir <- file.path(script_dir, "aligned_data")
downstream_dir <- file.path(script_dir, "exploratory_downstream")
dir.create(downstream_dir, recursive = TRUE, showWarnings = FALSE)

views <- c("Transcriptomics", "Proteomics")

# User setting
model_name <- "SCOT+"  # Choose from "Paired", "Unpaired", "SCOT+"

allowed_models <- c("Paired", "Unpaired", "SCOT+")
if (!(model_name %in% allowed_models)) {
  stop(sprintf("Warning: invalid model '%s'. Choose from: Paired, Unpaired, or SCOT+.", model_name))
}
cat("Model input correct!\n")

# Functions
attach_metadata_to_model <- function(mofa_model, metadata, sample_col = "SampleID") {
  model_samples <- unname(unlist(samples_names(mofa_model)))
  
  metadata <- metadata %>%
    filter(.data[[sample_col]] %in% model_samples) %>%
    distinct(.data[[sample_col]], .keep_all = TRUE)
  
  metadata <- metadata[match(model_samples, metadata[[sample_col]]), , drop = FALSE]
  
  if (nrow(metadata) != length(model_samples)) {
    stop("Mismatch between model samples and metadata rows after alignment.")
  }
  if (!all(metadata[[sample_col]] == model_samples)) {
    stop("Sample order mismatch between model and metadata.")
  }
  
  metadata_for_join <- metadata %>%
    rename(sample = all_of(sample_col))
  
  mofa_model@samples_metadata <- mofa_model@samples_metadata %>%
    left_join(metadata_for_join, by = "sample")
  
  print(head(mofa_model@samples_metadata))
  print(colnames(mofa_model@samples_metadata))
  
  return(mofa_model)
}

save_data_overview <- function(mofa_model, downstream_dir) {
  p_overview <- plot_data_overview(mofa_model)
  ggsave(filename = file.path(downstream_dir, "Data_overview.svg"),
         plot = p_overview, device = svglite, width = 8, height = 6)
}

save_variance_plots <- function(mofa_model, downstream_dir) {
  if (!is.null(mofa_model@cache$variance_explained)) {
    print(mofa_model@cache$variance_explained$r2_total)
    print(mofa_model@cache$variance_explained$r2_per_factor)
  } else {
    message("Variance explained cache is empty; plots will compute from model.")
  }
  
  p_var <- plot_variance_explained(mofa_model, x = "view", y = "factor", plot_total = FALSE)
  ggsave(filename = file.path(downstream_dir, "Variance_plot.svg"),
         plot = p_var, device = svglite, width = 8, height = 6)
  
  p_var_total <- plot_variance_explained(mofa_model, x = "view", y = "factor", plot_total = TRUE)
  ggsave(filename = file.path(downstream_dir, "Variance_explained_total.svg"),
         plot = p_var_total, device = svglite, width = 8, height = 6)
}

save_factor_plots_shared <- function(mofa_model, downstream_dir, factors_to_plot) {
  p_factor <- plot_factor(mofa_model, factors = factors_to_plot, color_by = "CellLine", dot_size = 1) +
    theme_bw(base_size = 12) +
    scale_colour_okabe_ito() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12), legend.text = element_text(size = 12))
  
  ggsave(filename = file.path(downstream_dir, "Factor_by_CellLine.svg"),
         plot = p_factor, device = svglite, width = 6, height = 4)
  
  p_violin <- plot_factor(mofa_model, factors = factors_to_plot, color_by = "CellLine",
                          dot_size = 3, dodge = TRUE, legend = FALSE, add_violin = TRUE,
                          violin_alpha = 0.25) +
    scale_color_manual(values = c("C10" = "black", "SVEC" = "red")) +
    scale_fill_manual(values = c("C10" = "black", "SVEC" = "red"))
  
  print(p_violin)
  
  ggsave(filename = file.path(downstream_dir, "Factor_violin.svg"),
         plot = p_violin, device = svglite, width = 8, height = 6)
  
  p_factors_grid <- eval(bquote(plot_factors(.(mofa_model), factors = .(as.integer(factors_to_plot)),
                                             color_by = "CellLine")))
  
  print(p_factors_grid)
  
  ggsave(filename = file.path(downstream_dir, "Factors_grid.svg"),
         plot = p_factors_grid, device = svglite, width = 10, height = 8)
}

save_feature_level_outputs <- function(mofa_model, downstream_dir, views, K) {
  for (v in views) {
    for (f in seq_len(K)) {
      message("Processing: ", v, " | Factor ", f)
      
      p_weights <- plot_weights(mofa_model, view = v, factor = f,
                                nfeatures = 10, scale = TRUE, abs = FALSE, text_size = 2)
      ggsave(file.path(downstream_dir, paste0(v, "_F", f, "_weights.svg")),
             p_weights, width = 6, height = 4, device = svglite)
      
      p_top <- plot_top_weights(mofa_model, view = v, factor = f, nfeatures = 20)
      ggsave(file.path(downstream_dir, paste0(v, "_F", f, "_topweights.svg")),
             p_top, width = 6, height = 4, device = svglite)
      
      weights_df <- get_weights(mofa_model, views = v, factors = f, as.data.frame = TRUE)
      
      if (!"feature" %in% colnames(weights_df)) {
        stop("Expected column 'feature' not found in weights table.")
      }
      if (!"value" %in% colnames(weights_df)) {
        stop("Expected column 'value' not found in weights table.")
      }
      
      weights_df <- weights_df %>%
        mutate(feature_clean = sub("_(rna|prot)$", "", feature), abs_value = abs(value)) %>%
        arrange(desc(abs_value))
      
      write.csv(weights_df, file = file.path(downstream_dir, paste0(v, "_F", f, "_weights_table.csv")),
                row.names = FALSE, quote = FALSE)
      
      top_weights_df <- weights_df %>%
        slice_head(n = 20)
      
      write.csv(top_weights_df,
                file = file.path(downstream_dir, paste0(v, "_F", f, "_top20_weights_table.csv")),
                row.names = FALSE, quote = FALSE)
      
      p_scatter <- plot_data_scatter(mofa_model, view = v, factor = f, features = 12,
                                     add_lm = TRUE, color_by = "CellLine")
      ggsave(file.path(downstream_dir, paste0(v, "_F", f, "_scatter.svg")),
             p_scatter, width = 10, height = 8, device = svglite)
      
      svglite(file.path(downstream_dir, paste0(v, "_F", f, "_heatmap.svg")),
              width = 10, height = 8)
      
      plot_data_heatmap(mofa_model, view = v, factor = f, features = 20, cluster_rows = TRUE,
                        cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
                        fontsize_col = 6)
      dev.off()
    }
  }
}

save_umap_shared <- function(mofa_model, downstream_dir, K) {
  set.seed(42)
  factors_umap <- seq_len(min(3, K))
  mofa_model <- run_umap(mofa_model, factors = factors_umap)
  
  p_umap <- plot_dimred(mofa_model, method = "UMAP", color_by = "CellLine", color_name = "Cell Line",
    shape_by = NULL, label = FALSE, dot_size = 2, legend = TRUE) +
    scale_colour_okabe_ito() +
    theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16),
          legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  ggsave(filename = file.path(downstream_dir, "UMAP_cellline.svg"),
         plot = p_umap, device = svglite, width = 8, height = 6)
  
  return(mofa_model)
}

run_full_downstream <- function(mofa_model, downstream_dir, views) {
  K <- get_dimensions(mofa_model)$K
  factors_to_plot <- seq_len(min(5, K))
  
  save_data_overview(mofa_model, downstream_dir)
  save_variance_plots(mofa_model, downstream_dir)
  save_factor_plots_shared(mofa_model, downstream_dir, factors_to_plot)
  save_feature_level_outputs(mofa_model, downstream_dir, views, K)
  mofa_model <- save_umap_shared(mofa_model, downstream_dir, K)
  
  return(mofa_model)
}

run_unpaired_downstream <- function(mofa_model, downstream_dir) {
  K <- get_dimensions(mofa_model)$K
  factors_to_plot <- seq_len(min(5, K))
  
  p_factor_cellline <- plot_factor(mofa_model, factors = factors_to_plot,
                                   color_by = "CellLine", dot_size = 1) +
    scale_colour_okabe_ito() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  ggsave(filename = file.path(downstream_dir, "Factor_by_CellLine.svg"),
         plot = p_factor_cellline, device = svglite, width = 7, height = 5)
  
  p_factor_modality <- plot_factor(mofa_model, factors = factors_to_plot,
                                   color_by = "Modality", dot_size = 1) +
    
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  ggsave(filename = file.path(downstream_dir, "Factor_by_Modality.svg"),
         plot = p_factor_modality, device = svglite, width = 7, height = 5)
  
  set.seed(42)
  factors_umap <- seq_len(min(3, K))
  mofa_model <- run_umap(mofa_model, factors = factors_umap)
  
  p_umap_cellline <- plot_dimred(mofa_model, method = "UMAP", color_by = "CellLine",
                                 shape_by = "Modality", label = FALSE, dot_size = 2, legend = TRUE) +
    scale_colour_okabe_ito()
  
  ggsave(filename = file.path(downstream_dir, "UMAP_CellLine_shape_Modality.svg"),
         plot = p_umap_cellline, device = svglite, width = 8, height = 6)
  
  p_umap_modality <- plot_dimred(mofa_model, method = "UMAP", color_by = "Modality", shape_by = "CellLine",
    label = FALSE, dot_size = 2, legend = TRUE)
  
  ggsave(filename = file.path(downstream_dir, "UMAP_Modality_shape_CellLine.svg"),
         plot = p_umap_modality, device = svglite, width = 8, height = 6)
  
  return(mofa_model)
}

# Main
if (model_name %in% c("Paired", "Unpaired")) {
  rna_metadata <- read.csv(file.path(meta_dir, "Processed_transcriptomics_metadata.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)
  
  prot_metadata <- read.csv(file.path(meta_dir, "Processed_proteomics_metadata.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
}

if (model_name == "Paired") {
  mofa_model <- load_model(file.path(model_dir, "Paired_model.hdf5"))
  
  metadata <- rna_metadata %>%
    distinct(SampleID, .keep_all = TRUE)
  
  mofa_model <- attach_metadata_to_model(mofa_model, metadata = metadata, sample_col = "SampleID")
  
  mofa_model <- run_full_downstream(mofa_model = mofa_model, downstream_dir = downstream_dir,
                                    views = views)
  
} else if (model_name == "SCOT+") {
  aligned_metadata <- read.csv(file.path(align_dir, "aligned_metadata.csv"),
                               stringsAsFactors = FALSE, check.names = FALSE)
  
  mofa_model <- load_model(file.path(model_dir, "SCOT+-aligned_model.hdf5"))
  
  mofa_model <- attach_metadata_to_model(mofa_model, metadata = aligned_metadata,
                                         sample_col = "SampleID")
  
  mofa_model <- run_full_downstream(mofa_model = mofa_model, downstream_dir = downstream_dir,
                                    views = views)
  
} else if (model_name == "Unpaired") {
  mofa_model <- load_model(file.path(model_dir, "Unpaired_model.hdf5"))
  
  rna_metadata$SampleID <- paste0(as.character(rna_metadata$SampleID), "_rna")
  prot_metadata$SampleID <- paste0(as.character(prot_metadata$SampleID), "_prot")
  rna_metadata$Modality <- "Transcriptomics"
  prot_metadata$Modality <- "Proteomics"
  
  metadata <- bind_rows(prot_metadata, rna_metadata) %>%
    distinct(SampleID, .keep_all = TRUE)
  
  mofa_model <- attach_metadata_to_model(mofa_model, metadata = metadata, sample_col = "SampleID")
  
  mofa_model <- run_unpaired_downstream(mofa_model = mofa_model, downstream_dir = downstream_dir)
}
