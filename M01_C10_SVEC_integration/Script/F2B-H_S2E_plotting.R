# =========================================================
# Package versions
# =========================================================
# ggplot2 4.0.0
# dplyr 1.1.4
# tidyr 1.3.1
# MOFA2 1.14.0
# svglite 2.2.1
# extrafont 0.20
# cowplot 1.2.0
# uwot 0.2.3
# pheatmap 1.0.13
# clue 0.3.66
# grid 4.4.3
# ggplotify 0.1.3
# RColorBrewer 1.1.3
# ggVennDiagram 1.5.7
# patchwork 1.3.2
# grDevices 4.4.3

library(ggplot2)
library(dplyr)
library(tidyr)
library(MOFA2)
library(svglite)
library(extrafont)
library(cowplot)
library(uwot)
library(pheatmap)
library(clue)
library(grid)
library(ggplotify)
library(RColorBrewer)
library(ggVennDiagram)
library(patchwork)
library(grDevices)

FONT_PT <- 8.5
FONT_MM <- FONT_PT / ggplot2::.pt
FONT_FAMILY <- "Arial"

theme_pub <- function() {
  theme_bw(base_size = FONT_PT, base_family = FONT_FAMILY) +
    theme(text = element_text(family = FONT_FAMILY, size = FONT_PT),
          axis.text = element_text(family = FONT_FAMILY, size = FONT_PT),
          axis.title = element_text(family = FONT_FAMILY, size = FONT_PT),
          legend.text = element_text(family = FONT_FAMILY, size = FONT_PT),
          legend.title = element_text(family = FONT_FAMILY, size = FONT_PT),
          strip.text = element_text(family = FONT_FAMILY, size = FONT_PT),
          plot.title = element_text(family = FONT_FAMILY, size = FONT_PT),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

theme_set(theme_pub())

# Paths
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
ref_model_path  <- file.path(script_dir, "trained_model/Paired_model.hdf5")
scot_model_path <- file.path(script_dir, "trained_model/SCOT+-aligned_model.hdf5")
ref_rna_meta_path  <- file.path(script_dir, "processed_data/processed_transcriptomics_metadata.csv")
ref_prot_meta_path <- file.path(script_dir, "processed_data/processed_proteomics_metadata.csv")
scot_meta_path     <- file.path(script_dir, "aligned_data/aligned_metadata.csv")

outdir <- file.path(script_dir, "figure2B-H_S2E")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load models
ref_model  <- load_model(ref_model_path)
scot_model <- load_model(scot_model_path)

# Attach metadata to reference model
ref_rna_metadata <- read.csv(ref_rna_meta_path, stringsAsFactors = FALSE)
ref_prot_metadata <- read.csv(ref_prot_meta_path, stringsAsFactors = FALSE)

ref_model_samples <- unname(unlist(samples_names(ref_model)))

ref_metadata <- ref_rna_metadata %>%
  filter(SampleID %in% ref_model_samples) %>%
  distinct(SampleID, .keep_all = TRUE)

ref_metadata <- ref_metadata[match(ref_model_samples, ref_metadata$SampleID), , drop = FALSE]

if (nrow(ref_metadata) != length(ref_model_samples)) {
  stop("Mismatch between reference model samples and metadata rows.")
}
if (!all(ref_metadata$SampleID == ref_model_samples)) {
  stop("Reference metadata sample order mismatch.")
}

ref_metadata_for_join <- ref_metadata %>%
  rename(sample = SampleID)

ref_model@samples_metadata <- ref_model@samples_metadata %>%
  left_join(ref_metadata_for_join, by = "sample")

# Attach metadata to SCOT+ model
scot_metadata <- read.csv(scot_meta_path, stringsAsFactors = FALSE, check.names = FALSE)

if (colnames(scot_metadata)[1] %in% c("", "X", "Unnamed: 0")) {
  colnames(scot_metadata)[1] <- "SampleID"
}

scot_model_samples <- unname(unlist(samples_names(scot_model)))

scot_metadata <- scot_metadata %>%
  filter(SampleID %in% scot_model_samples) %>%
  distinct(SampleID, .keep_all = TRUE)

scot_metadata <- scot_metadata[match(scot_model_samples, scot_metadata$SampleID), , drop = FALSE]

if (nrow(scot_metadata) != length(scot_model_samples)) {
  stop("Mismatch between SCOT+ model samples and metadata rows.")
}
if (!all(scot_metadata$SampleID == scot_model_samples)) {
  stop("SCOT+ metadata sample order mismatch.")
}

scot_metadata_for_join <- scot_metadata %>%
  rename(sample = SampleID)

scot_model@samples_metadata <- scot_model@samples_metadata %>%
  left_join(scot_metadata_for_join, by = "sample")

# Helper: SVG export
save_svg <- function(filename, plot, width, height) {
  ggsave(filename = filename, plot = plot, device = svglite::svglite, width = width,
         height = height,  units = "in", system_fonts = list(sans = "Arial", Arial = "Arial"))
}

# Helper: variance explained extraction
extract_var_explained_df <- function(model, model_label) {
  var_exp <- get_variance_explained(model)
  r2_list <- var_exp$r2_per_factor
  
  df_list <- lapply(names(r2_list), function(grp) {
    x <- as.data.frame(r2_list[[grp]])
    x$Factor <- rownames(x)
    
    x_long <- x %>%
      tidyr::pivot_longer(cols = -Factor, names_to = "View", values_to = "R2") %>%
      mutate(Group = grp, Model = model_label)
    
    x_long
  })
  
  df <- bind_rows(df_list) %>%
    mutate(View = factor(View, levels = c("Transcriptomics", "Proteomics")))
  
  return(df)
}

# Helper: factor value extraction
extract_factor_df <- function(model, model_label) {
  factor_df <- get_factors(model, factors = "all", as.data.frame = TRUE)
  
  expected_cols <- c("sample", "factor", "value")
  if (!all(expected_cols %in% colnames(factor_df))) {
    stop("Factor table does not contain expected columns: sample, factor, value")
  }
  
  meta_df <- model@samples_metadata %>%
    dplyr::select(sample, CellLine)
  
  factor_df <- factor_df %>%
    left_join(meta_df, by = "sample") %>%
    mutate(Model = model_label,
           Factor = factor(factor, levels = paste0("Factor",
                                                   seq_len(max(as.integer(gsub("Factor", "", unique(factor))))))))
  
  return(factor_df)
}

# Helper: save pheatmap gtable to SVG
save_pheatmap_svg <- function(ph, filename, width = 5, height = 4) {
  svglite::svglite(filename, width = width, height = height)
  grid::grid.draw(ph$gtable)
  dev.off()
}

# Helper: build combined matched-factor heatmap
build_matched_factor_heatmap <- function(ref_model, scot_model, matched_pairs,
                                         ref_factor_name = c("Factor1", "Factor2", "Factor3"),
                                         view_name = c("Transcriptomics", "Proteomics"),
                                         nfeatures = 20,
                                         force_feature = NULL) {
  ref_factor_name <- match.arg(ref_factor_name)
  view_name <- match.arg(view_name)
  
  W_ref  <- get_weights(ref_model,  views = view_name, factors = "all")[[1]]
  W_scot <- get_weights(scot_model, views = view_name, factors = "all")[[1]]
  
  D_ref  <- get_data(ref_model,  views = view_name)[[1]]
  D_scot <- get_data(scot_model, views = view_name)[[1]]
  
  if (is.list(D_ref))  D_ref  <- D_ref[[1]]
  if (is.list(D_scot)) D_scot <- D_scot[[1]]
  
  D_ref  <- as.matrix(D_ref)
  D_scot <- as.matrix(D_scot)
  
  if (!ref_factor_name %in% matched_pairs$ReferenceFactor) {
    stop(paste("Reference", ref_factor_name, "not found in matched_pairs."))
  }
  
  row_match <- matched_pairs %>%
    filter(ReferenceFactor == ref_factor_name) %>%
    slice(1)
  
  scot_factor_name <- row_match$SCOTFactor[[1]]
  sign_flip <- sign(row_match$Correlation[[1]])
  if (sign_flip == 0) sign_flip <- 1
  
  ref_idx  <- match(ref_factor_name, colnames(W_ref))
  scot_idx <- match(scot_factor_name, colnames(W_scot))
  
  if (is.na(ref_idx) || is.na(scot_idx)) {
    stop("Could not match factor indices for heatmap construction.")
  }
  
  common_features <- intersect(rownames(W_ref), rownames(W_scot))
  common_features <- intersect(common_features, intersect(rownames(D_ref), rownames(D_scot)))
  
  w_ref_factor <- W_ref[common_features, ref_idx]
  top_ref <- names(sort(abs(w_ref_factor), decreasing = TRUE))
  
  if (!is.null(force_feature) && force_feature %in% common_features) {
    top_ref <- c(force_feature, setdiff(top_ref, force_feature))
  }
  
  features_use <- head(top_ref, nfeatures)
  
  Z_ref  <- get_factors(ref_model,  factors = "all", as.data.frame = FALSE)[[1]]
  Z_scot <- get_factors(scot_model, factors = "all", as.data.frame = FALSE)[[1]]
  
  ref_factor_values  <- Z_ref[, ref_factor_name]
  scot_factor_values <- Z_scot[, scot_factor_name] * sign_flip
  
  ref_order  <- names(sort(ref_factor_values, decreasing = TRUE))
  scot_order <- names(sort(scot_factor_values, decreasing = TRUE))
  
  ref_order  <- intersect(ref_order, colnames(D_ref))
  scot_order <- intersect(scot_order, colnames(D_scot))
  
  mat_ref  <- D_ref[features_use, ref_order, drop = FALSE]
  mat_scot <- D_scot[features_use, scot_order, drop = FALSE]
  
  colnames(mat_ref)  <- paste0("Reference__", colnames(mat_ref))
  colnames(mat_scot) <- paste0("SCOTplus__", colnames(mat_scot))
  
  mat_ref_scaled  <- t(scale(t(mat_ref)))
  mat_scot_scaled <- t(scale(t(mat_scot)))
  
  mat_ref_scaled[is.na(mat_ref_scaled)] <- 0
  mat_scot_scaled[is.na(mat_scot_scaled)] <- 0
  
  mat_scaled <- cbind(mat_ref_scaled, mat_scot_scaled)
  
  ann_ref <- ref_model@samples_metadata %>%
    dplyr::select(sample, CellLine) %>%
    filter(sample %in% ref_order) %>%
    distinct(sample, .keep_all = TRUE)
  
  ann_scot <- scot_model@samples_metadata %>%
    dplyr::select(sample, CellLine) %>%
    filter(sample %in% scot_order) %>%
    distinct(sample, .keep_all = TRUE)
  
  ann_ref <- ann_ref[match(ref_order, ann_ref$sample), , drop = FALSE]
  ann_scot <- ann_scot[match(scot_order, ann_scot$sample), , drop = FALSE]
  
  annotation_col <- rbind(
    data.frame(sample = paste0("Reference__", ann_ref$sample), Model = "Reference",
               CellLine = ann_ref$CellLine),
    data.frame(sample = paste0("SCOTplus__", ann_scot$sample), Model = "SCOT+-aligned",
               CellLine = ann_scot$CellLine))
  
  rownames(annotation_col) <- annotation_col$sample
  annotation_col$sample <- NULL
  
  annotation_colors <- list(Model = c("Reference" = "#7F7F7F", "SCOT+-aligned" = "#377EB8"),
                            CellLine = c("C10" = "#000000", "SVEC" = "#D55E00"))
  
  heat_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  heat_breaks <- seq(-2.5, 2.5, length.out = 101)
  
  ph <- pheatmap::pheatmap(mat_scaled, cluster_rows = TRUE, cluster_cols = FALSE,
                           show_rownames = TRUE, show_colnames = FALSE, annotation_col = annotation_col,
                           annotation_colors = annotation_colors, gaps_col = ncol(mat_ref),
                           fontsize = FONT_PT, fontsize_row = FONT_PT, fontsize_col = FONT_PT,
                           fontfamily = FONT_FAMILY, silent = TRUE, border_color = NA, 
                           color = heat_colors, breaks = heat_breaks)
  
  return(ph)
}

extract_pheatmap_grob <- function(ph, grob_name) {
  idx <- which(ph$gtable$layout$name == grob_name)
  if (length(idx) == 0) return(grid::nullGrob())
  ph$gtable$grobs[[idx[1]]]
}

remove_pheatmap_all_legends <- function(ph) {
  legend_names <- c("legend", "annotation_legend")
  keep <- !ph$gtable$layout$name %in% legend_names
  
  ph$gtable$grobs <- ph$gtable$grobs[keep]
  ph$gtable$layout <- ph$gtable$layout[keep, , drop = FALSE]
  
  ph
}

remove_pheatmap_legend <- function(ph) {
  keep <- ph$gtable$layout$name != "annotation_legend"
  ph$gtable$grobs <- ph$gtable$grobs[keep]
  ph$gtable$layout <- ph$gtable$layout[keep, , drop = FALSE]
  ph
}

# Extract plotting data
var_ref  <- extract_var_explained_df(ref_model, "Reference")
var_scot <- extract_var_explained_df(scot_model, "SCOT+-aligned")

var_df <- bind_rows(var_ref, var_scot)

factor_ref  <- extract_factor_df(ref_model, "Reference")
factor_scot <- extract_factor_df(scot_model, "SCOT+-aligned")

factor_df <- bind_rows(factor_ref, factor_scot)

all_factor_levels <- paste0("Factor",
                            seq_len(max(as.integer(gsub("Factor", "", as.character(unique(factor_df$Factor)))))))
factor_df$Factor <- factor(factor_df$Factor, levels = all_factor_levels)
var_df$Factor <- factor(var_df$Factor, levels = all_factor_levels)

factor_df <- factor_df %>%
  mutate(Model_CellLine = paste(Model, CellLine, sep = "_"))

factor_df$Model_CellLine <- factor(factor_df$Model_CellLine,
                                   levels = c("Reference_C10", "Reference_SVEC", "SCOT+-aligned_C10", "SCOT+-aligned_SVEC"))

# Figure 2B: variance explained
var_df <- var_df %>%
  mutate(Model = factor(Model, levels = c("Reference", "SCOT+-aligned")),
         Factor = factor(Factor, levels = c("Factor1", "Factor2", "Factor3", "Factor4")),
         View = factor(View, levels = c("Transcriptomics", "Proteomics")))

factor_pos <- data.frame(Model = c(rep("Reference", 4), rep("SCOT+-aligned", 3)),
                         Factor = c("Factor1", "Factor2", "Factor3", "Factor4",
                                    "Factor1", "Factor2", "Factor3"),
                         x = c(1, 2, 3, 4, 5.6, 6.6, 7.6), stringsAsFactors = FALSE)

factor_pos <- factor_pos %>%
  mutate(Model = factor(Model, levels = c("Reference", "SCOT+-aligned")),
         Factor = factor(Factor, levels = c("Factor1", "Factor2", "Factor3", "Factor4")))

var_df <- var_df %>%
  dplyr::select(-dplyr::any_of(c("x_ref", "x_scot", "x", "Factor_Model"))) %>%
  left_join(factor_pos, by = c("Model", "Factor")) %>%
  filter(!is.na(x))

p2B <- ggplot(var_df, aes(x = x, y = R2, fill = View)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Transcriptomics" = "#E69F00", "Proteomics" = "#56B4E9")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5.6, 6.6, 7.6),
                     labels = c("Factor1", "Factor2", "Factor3", "Factor4",
                                "Factor1", "Factor2", "Factor3"), limits = c(0.5, 8.5),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Variance explained (%)", fill = NULL, colour = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(family = FONT_FAMILY, size = FONT_PT, margin = margin(t = 4)),
    axis.title.y = element_text(family = FONT_FAMILY, size = FONT_PT, margin = margin(r = 8)),
    legend.position = c(0.837, 0.88),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.key = element_blank(),
    legend.margin = margin(3, 3, 3, 3),
    legend.box.margin = margin(1, 1, 1, 1),
    legend.spacing.y = unit(0.05, "cm"),
    legend.key.size = unit(0.32, "cm"),
    plot.margin = margin(t = 8, r = 8, b = 38, l = 8)) +
  annotate("text", x = mean(c(1, 4)), y = -Inf, label = "Reference paired model", vjust = 4.0,
           size = FONT_MM, family = FONT_FAMILY) +
  annotate("text", x = mean(c(5.6, 7.6)), y = -Inf, label = "SCOT+-aligned model", vjust = 4.0,
           size = FONT_MM, family = FONT_FAMILY)

save_svg(filename = file.path(outdir, "Figure2B_variance_explained.svg"), plot = p2B, width = 4.5,
         height = 3.2)

write.csv(var_df, file.path(outdir, "variance_explained.csv"))

# Figure 2C: combined factor dot plot
cellline_colors <- c("C10" = "#000000", "SVEC" = "#D55E00")

factor_df <- factor_df %>%
  mutate(Model = factor(Model, levels = c("Reference", "SCOT+-aligned")),
         Factor = factor(Factor, levels = c("Factor1", "Factor2", "Factor3", "Factor4")))

factor_pos_B <- data.frame(Model = c(rep("Reference", 4), rep("SCOT+-aligned", 3)),
                           Factor = c("Factor1", "Factor2", "Factor3", "Factor4",
                                      "Factor1", "Factor2", "Factor3"),
                           x = c(1, 1.7, 2.4, 3.1, 4.1, 4.8, 5.5), stringsAsFactors = FALSE)

factor_pos_B <- factor_pos_B %>%
  mutate(Model = factor(Model, levels = c("Reference", "SCOT+-aligned")),
         Factor = factor(Factor, levels = c("Factor1", "Factor2", "Factor3", "Factor4")))

factor_df <- factor_df %>%
  dplyr::select(-dplyr::any_of(c("x_ref", "x_scot", "x", "Factor_Model"))) %>%
  left_join(factor_pos_B, by = c("Model", "Factor")) %>%
  filter(!is.na(x))

p2C <- ggplot(factor_df, aes(x = x, y = value, colour = CellLine)) +
  geom_jitter(width = 0.12, height = 0, size = 1.2, alpha = 0.7) +
  scale_colour_manual(values = cellline_colors) +
  scale_x_continuous(breaks = c(1, 1.7, 2.4, 3.1, 4.1, 4.8, 5.5),
                     labels = c("Factor1", "Factor2", "Factor3", "Factor4",
                                "Factor1", "Factor2", "Factor3"), limits = c(0.5, 6.7),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Factor value", colour = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(family = FONT_FAMILY, size = FONT_PT, margin = margin(t = 4)),
        axis.title.y = element_text(family = FONT_FAMILY, size = FONT_PT, margin = margin(r = 8)),
    # Legend inside the dot plot box, top-right
    legend.position = c(0.91, 0.88),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.key = element_blank(),
    legend.margin = margin(3, 3, 3, 3),
    legend.box.margin = margin(1, 1, 1, 1),
    legend.spacing.y = unit(0.08, "cm"),
    legend.key.size = unit(0.32, "cm"),
    
    # Extra bottom space for model labels
    plot.margin = margin(t = 8, r = 8, b = 38, l = 8)) +
  annotate("text", x = mean(c(1, 3.1)), y = -Inf, label = "Reference paired model", vjust = 4.0,
           size = FONT_MM, family = FONT_FAMILY) +
  annotate("text", x = mean(c(4.1, 5.5)), y = -Inf, label = "SCOT+-aligned model", vjust = 4.0,
           size = FONT_MM, family = FONT_FAMILY)

save_svg(filename = file.path(outdir, "Figure2C_factor_values.svg"), plot = p2C, width = 5.0,
         height = 3.2)

write.csv(factor_df, file.path(outdir, "factor_values.csv"))

# Figure 2D: sample-level concordance heatmap
Z_ref  <- get_factors(ref_model, factors = "all", as.data.frame = FALSE)[[1]]
Z_scot <- get_factors(scot_model, factors = "all", as.data.frame = FALSE)[[1]]

if (is.null(rownames(Z_ref)) || is.null(rownames(Z_scot))) {
  stop("Factor matrices must have sample names as rownames.")
}

common_cells <- intersect(rownames(Z_ref), rownames(Z_scot))
if (length(common_cells) == 0) {
  stop("No shared sample IDs between the reference and SCOT+ MOFA models.")
}

Z_ref_aligned  <- Z_ref[match(common_cells, rownames(Z_ref)), , drop = FALSE]
Z_scot_aligned <- Z_scot[match(common_cells, rownames(Z_scot)), , drop = FALSE]

stopifnot(identical(rownames(Z_ref_aligned), rownames(Z_scot_aligned)))

cor_mat <- cor(Z_ref_aligned, Z_scot_aligned, use = "pairwise.complete.obs", method = "pearson")

n_ref  <- nrow(cor_mat)
n_scot <- ncol(cor_mat)

if (n_ref <= n_scot) {
  cost_mat <- 1 - abs(cor_mat)
  assignment <- clue::solve_LSAP(cost_mat)
  
  matched_pairs <- data.frame(ReferenceFactor = rownames(cor_mat), 
                              SCOTFactor = colnames(cor_mat)[assignment],
                              Correlation = cor_mat[cbind(seq_len(n_ref), assignment)],
                              stringsAsFactors = FALSE)
  
  cor_mat_ordered <- cor_mat[, assignment, drop = FALSE]
  colnames(cor_mat_ordered) <- matched_pairs$SCOTFactor
  
} else {
  cost_mat_t <- 1 - abs(t(cor_mat))
  assignment_t <- clue::solve_LSAP(cost_mat_t)
  
  matched_pairs <- data.frame(ReferenceFactor = rownames(cor_mat)[assignment_t],
                              SCOTFactor = colnames(cor_mat),
                              Correlation = cor_mat[cbind(assignment_t, seq_len(n_scot))],
                              stringsAsFactors = FALSE)
  
  cor_mat_ordered <- cor_mat[assignment_t, , drop = FALSE]
  rownames(cor_mat_ordered) <- matched_pairs$ReferenceFactor
}

write.csv(cor_mat_ordered, file.path(outdir, "Figure2D_sample_level_concordance.csv"))

write.csv(matched_pairs, file.path(outdir, "Figure2_matched_factor_pairs.csv"), row.names = FALSE)

breaks <- seq(-1, 1, length.out = 101)
colors <- colorRampPalette(c("blue", "white", "red"))(100)

make_corr_scale <- function(colors, font_size = FONT_PT, font_family = FONT_FAMILY) {
  scale_df <- data.frame(ymin = seq(-1, 1, length.out = 200),
                         ymax = c(seq(-1, 1, length.out = 200)[-1], 1),
                         value = seq(-1, 1, length.out = 200))
  
  ggplot(scale_df) +
    geom_rect(aes(xmin = 0.00, xmax = 0.3, ymin = ymin, ymax = ymax, fill = value), colour = NA) +
    scale_fill_gradientn(colours = colors, limits = c(-1, 1), guide = "none") +
    scale_y_continuous(position = "right", breaks = c(-1, -0.5, 0, 0.5, 1),
                       labels = c("-1", "-0.5", "0", "0.5", "1"), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0.00, 0.48), expand = c(0, 0)) +
    coord_cartesian(ylim = c(-1, 1), clip = "off") +
    theme_void(base_family = font_family) +
    theme(axis.text.y = element_text(family = font_family, size = font_size, colour = "black",
                                     margin = margin(l = 1)),
          axis.ticks.y = element_blank(),
          axis.ticks.length = unit(2, "pt"), plot.margin = margin(0, 0, 0, 0))
}

add_right_corr_scale <- function(heatmap_plot, scale_plot, scale_x = 0.82, scale_y = 0.31,
                                 scale_width = 0.05, scale_height = 0.38) {
  cowplot::ggdraw() +
    cowplot::draw_plot(heatmap_plot, x = 0.00, y = 0.00, width = 0.78, height = 1.00) +
    cowplot::draw_plot(scale_plot, x = scale_x, y = scale_y, width = scale_width,
                       height = scale_height)
}

corr_scale <- make_corr_scale(colors)


ph2D <- pheatmap::pheatmap(cor_mat_ordered, cluster_rows = FALSE, cluster_cols = FALSE, color = colors,
                           breaks = breaks, fontsize = FONT_PT, fontsize_row = FONT_PT,
                           fontsize_col = FONT_PT, fontsize_number = FONT_PT, display_numbers = TRUE,
                           angle_col = 0, silent = TRUE, border_color = NA, legend = FALSE, 
                           legend_breaks = c(-1, -0.5, 0, 0.5, 1),
                           legend_labels = c("-1", "-0.5", "0", "0.5", "1"), fontfamily = FONT_FAMILY)

p2D_heatmap <- ggplotify::as.ggplot(ph2D$gtable) +
  theme(plot.margin = margin(t = 8, r = 8, b = 8, l = 8))

p2D <- add_right_corr_scale(heatmap_plot = p2D_heatmap, scale_plot = corr_scale, scale_x = 0.82,
                            scale_y = 0.31, scale_width = 0.085, scale_height = 0.38)

save_svg(filename = file.path(outdir, "Figure2D_sample_level_concordance.svg"), plot = p2D,
         width = 4.2, height = 3.2)

# Figure 2E: signed feature-level concordance heatmap
W_ref  <- get_weights(ref_model, views = "all", factors = "all")
W_scot <- get_weights(scot_model, views = "all", factors = "all")

common_genes <- intersect(rownames(W_ref[["Transcriptomics"]]), rownames(W_scot[["Transcriptomics"]]))

common_proteins <- intersect(rownames(W_ref[["Proteomics"]]), rownames(W_scot[["Proteomics"]]))

if (length(common_genes) == 0) stop("No shared transcriptomic features.")
if (length(common_proteins) == 0) stop("No shared proteomic features.")

rna_corr_signed  <- numeric(nrow(matched_pairs))
prot_corr_signed <- numeric(nrow(matched_pairs))

for (idx in seq_len(nrow(matched_pairs))) {
  ref_name  <- matched_pairs$ReferenceFactor[idx]
  scot_name <- matched_pairs$SCOTFactor[idx]
  
  sign_flip <- sign(matched_pairs$Correlation[idx])
  if (sign_flip == 0) sign_flip <- 1
  
  i <- match(ref_name, colnames(Z_ref_aligned))
  j <- match(scot_name, colnames(Z_scot_aligned))
  
  w1_rna  <- W_ref[["Transcriptomics"]][common_genes, i]
  w2_rna  <- W_scot[["Transcriptomics"]][common_genes, j] * sign_flip
  
  w1_prot <- W_ref[["Proteomics"]][common_proteins, i]
  w2_prot <- W_scot[["Proteomics"]][common_proteins, j] * sign_flip
  
  rna_corr_signed[idx]  <- cor(w1_rna,  w2_rna,  use = "pairwise.complete.obs", method = "pearson")
  prot_corr_signed[idx] <- cor(w1_prot, w2_prot, use = "pairwise.complete.obs", method = "pearson")
}

loading_mat_signed <- cbind(Transcriptomics = rna_corr_signed, Proteomics = prot_corr_signed)

rownames(loading_mat_signed) <- matched_pairs$ReferenceFactor

write.csv(loading_mat_signed, file.path(outdir, "Figure2E_feature_level_concordance_signed.csv"))

ph2E <- pheatmap::pheatmap(loading_mat_signed, cluster_rows = FALSE, cluster_cols = FALSE,
                           color = colors, breaks = breaks, fontsize = FONT_PT,
                           fontsize_row = FONT_PT, fontsize_col = FONT_PT, fontsize_number = FONT_PT,
                           display_numbers = TRUE, angle_col = 0, silent = TRUE, border_color = NA,
                           legend = FALSE, fontfamily = FONT_FAMILY)

p2E_heatmap <- ggplotify::as.ggplot(ph2E$gtable) +
  theme(plot.margin = margin(t = 8, r = 8, b = 8, l = 8))

p2E <- add_right_corr_scale(heatmap_plot = p2E_heatmap, scale_plot = corr_scale, scale_x = 0.82,
                            scale_y = 0.31, scale_width = 0.085, scale_height = 0.38)

save_svg(filename = file.path(outdir, "Figure2E_feature_level_concordance_signed.svg"), plot = p2E,
         width = 4.2, height = 3.2)

# Figure 2H: Factor 1 heatmaps
ph2H_rna <- build_matched_factor_heatmap(ref_model = ref_model, scot_model = scot_model,
                                         matched_pairs = matched_pairs, ref_factor_name = "Factor1",
                                         view_name = "Transcriptomics", nfeatures = 20,
                                         force_feature = "H2-K1_rna")

ph2H_prot <- build_matched_factor_heatmap(ref_model = ref_model, scot_model = scot_model,
                                          matched_pairs = matched_pairs, ref_factor_name = "Factor1",
                                          view_name = "Proteomics", nfeatures = 20,
                                          force_feature = "H2-K1_prot")

# Extract legends from the proteomics heatmap before removing them
heat_legend <- extract_pheatmap_grob(ph2H_prot, "legend")
annotation_legend <- extract_pheatmap_grob(ph2H_prot, "annotation_legend")

# Remove legends from both heatmaps
ph2H_rna_noleg  <- remove_pheatmap_all_legends(ph2H_rna)
ph2H_prot_noleg <- remove_pheatmap_all_legends(ph2H_prot)

g2H_rna  <- ggplotify::as.ggplot(ph2H_rna_noleg$gtable)
g2H_prot <- ggplotify::as.ggplot(ph2H_prot_noleg$gtable)

g2H_rna_tight <- g2H_rna +
  theme(plot.margin = margin(t = 5, r = 0, b = 5, l = 30, unit = "pt"))

# Keep proteomics heatmap same height; only shift horizontally as before
g2H_prot_shift <- cowplot::ggdraw() +
  cowplot::draw_plot(g2H_prot, x = -0.40, y = 0, width = 1.32, height = 1)

# Separate centred legend column, closer to proteomics heatmap
g2H_legend <- cowplot::ggdraw() +
  cowplot::draw_grob(heat_legend, x = -2.5, y = 0.54, width = 0.95, height = 0.28) +
  cowplot::draw_grob(annotation_legend, x = -2.5, y = 0.18, width = 0.95, height = 0.30)

fig2H_combined <- cowplot::plot_grid(g2H_rna_tight, g2H_prot_shift, g2H_legend, ncol = 3,
                                     rel_widths = c(1.5, 1, 0.2), align = "h", axis = "tb")

save_svg(filename = file.path(outdir, "Figure2H_Factor1_combined_heatmap.svg"), plot = fig2H_combined,
         width = 12, height = 5)

# Supplementary Figure S2E: Factor 2 feature heatmaps
phS2E_rna <- build_matched_factor_heatmap(ref_model = ref_model, scot_model = scot_model,
                                         matched_pairs = matched_pairs, ref_factor_name = "Factor2",
                                         view_name = "Transcriptomics", nfeatures = 20,
                                         force_feature = NULL)

phS2E_prot <- build_matched_factor_heatmap(ref_model = ref_model, scot_model = scot_model,
                                          matched_pairs = matched_pairs, ref_factor_name = "Factor2",
                                          view_name = "Proteomics", nfeatures = 20,
                                          force_feature = NULL)

# Extract legends from the proteomics heatmap before removing them
heat_legend_S2E <- extract_pheatmap_grob(phS2E_prot, "legend")
annotation_legend_S2E <- extract_pheatmap_grob(phS2E_prot, "annotation_legend")

# Remove legends from both heatmaps
phS2E_rna_noleg  <- remove_pheatmap_all_legends(phS2E_rna)
phS2E_prot_noleg <- remove_pheatmap_all_legends(phS2E_prot)

gS2E_rna  <- ggplotify::as.ggplot(phS2E_rna_noleg$gtable)
gS2E_prot <- ggplotify::as.ggplot(phS2E_prot_noleg$gtable)

gS2E_rna_tight <- gS2E_rna +
  theme(plot.margin = margin(t = 5, r = 0, b = 5, l = 30, unit = "pt"))

gS2E_prot_shift <- cowplot::ggdraw() +
  cowplot::draw_plot(gS2E_prot, x = -0.40, y = 0, width = 1.32, height = 1)

gS2E_legend <- cowplot::ggdraw() +
  cowplot::draw_grob(heat_legend_S2E, x = -2.5, y = 0.54, width = 0.95, height = 0.28) +
  cowplot::draw_grob(annotation_legend_S2E, x = -2.5, y = 0.18, width = 0.95, height = 0.30)

figS2E_combined <- cowplot::plot_grid(gS2E_rna_tight, gS2E_prot_shift, gS2E_legend, ncol = 3,
                                     rel_widths = c(1.5, 1, 0.2), align = "h", axis = "tb")

save_svg(filename = file.path(outdir, "S2E_Factor2_combined_heatmap.svg"),
         plot = figS2E_combined, width = 12, height = 5)

# Supplementary Figure 2G: joint UMAP on shared latent factors
ref_factor_long  <- get_factors(ref_model, factors = 1:3, as.data.frame = TRUE)
scot_factor_long <- get_factors(scot_model, factors = 1:3, as.data.frame = TRUE)

ref_factor_long <- ref_factor_long %>% dplyr::select(sample, factor, value) %>% mutate(Model = "Reference paired model")

scot_factor_long <- scot_factor_long %>% dplyr::select(sample, factor, value) %>% mutate(Model = "SCOT+-aligned model")

joint_factor_long <- bind_rows(ref_factor_long, scot_factor_long)

joint_factor_wide <- joint_factor_long %>% tidyr::pivot_wider(names_from = factor, values_from = value)

ref_meta_plot <- ref_model@samples_metadata %>% dplyr::select(sample, CellLine) %>% mutate(Model = "Reference paired model")

scot_meta_plot <- scot_model@samples_metadata %>% dplyr::select(sample, CellLine) %>% mutate(Model = "SCOT+-aligned model")

joint_meta <- bind_rows(ref_meta_plot, scot_meta_plot)

joint_plot_df <- joint_factor_wide %>% left_join(joint_meta, by = c("sample", "Model"))

if (any(is.na(joint_plot_df$CellLine))) {
  stop("Missing CellLine annotations in joint UMAP metadata.")
}

# Joint UMAP
umap_input <- joint_plot_df %>% dplyr::select(Factor1, Factor2, Factor3) %>% as.matrix()

set.seed(1)

umap_coords <- uwot::umap(umap_input, n_neighbors = 15, min_dist = 0.3, metric = "euclidean", verbose = TRUE)

joint_plot_df$UMAP1 <- umap_coords[, 1]
joint_plot_df$UMAP2 <- umap_coords[, 2]

cellline_colors_2G <- c("C10"  = "#0072B2", "SVEC" = "#D55E00")

model_shapes_2G <- c("Reference paired model" = 16, "SCOT+-aligned model" = 17)

p2G <- ggplot(joint_plot_df, aes(x = UMAP1, y = UMAP2, colour = CellLine, shape = Model)) +
  geom_point(size = 2.1, alpha = 0.65, stroke = 0) +
  scale_colour_manual(values = cellline_colors_2G) +
  scale_shape_manual(values = model_shapes_2G) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = NULL, shape = NULL) +
  theme_bw(base_size = FONT_PT, base_family = FONT_FAMILY) +
  theme(text = element_text(family = FONT_FAMILY, size = FONT_PT),
        axis.text = element_text(family = FONT_FAMILY, size = FONT_PT, colour = "black"),
        axis.title = element_text(family = FONT_FAMILY, size = FONT_PT, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        axis.ticks.length = unit(2, "pt"), legend.position = "right",
        legend.background = element_blank(), legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = FONT_FAMILY, size = FONT_PT),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 4),
        legend.spacing.y = unit(0.08, "cm"), plot.margin = margin(t = 8, r = 8, b = 8, l = 8)) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 3, alpha = 1, shape = 16)),
         shape = guide_legend(order = 2, override.aes = list(size = 3, alpha = 1, colour = "black")))

save_svg(filename = file.path(outdir, "Figure2G_joint_UMAP.svg"), plot = p2G, width = 6.0, height = 4.0)

# Supplementary 2F: combined signed Venn diagrams
F2F_outdir <- file.path(outdir, "Figure2F_Venn")
dir.create(supp5_outdir, recursive = TRUE, showWarnings = FALSE)

matched_pairs_2F <- matched_pairs[seq_len(min(3, nrow(matched_pairs))), , drop = FALSE]

views <- c("Transcriptomics", "Proteomics")
signs <- c("positive", "negative")
n_top <- 20
sign_colors <- c("positive" = "#D55E00", "negative" = "#0072B2")

# Helper: extract top features by sign
get_top_features_signed <- function(model, view, factor_name, sign = "positive", n_top = 20) {
  factor_idx <- as.integer(gsub("Factor", "", factor_name))
  
  w <- get_weights(model, views = view, factors = factor_idx, as.data.frame = TRUE)
  
  if (!all(c("feature", "value") %in% colnames(w))) {
    stop("Expected columns 'feature' and 'value' not found in weights table.")
  }
  
  w <- w %>%
    mutate(feature_clean = sub("_(rna|prot)$", "", feature, ignore.case = TRUE))
  
  if (sign == "positive") {
    w <- w %>% filter(value > 0) %>% arrange(desc(value))
  } else if (sign == "negative") {
    w <- w %>% filter(value < 0) %>% arrange(value)
  } else {
    stop("sign must be 'positive' or 'negative'")
  }
  
  if (nrow(w) == 0) return(character(0))
  unique(w$feature_clean[seq_len(min(n_top, nrow(w)))])
}

# Helper: clean Venn plot
make_small_venn_plot <- function(ref_features, scot_features, sign = "positive") {
  venn_list <- list(Reference = ref_features, `SCOT+` = scot_features)
  
  high_col <- sign_colors[[sign]]
  
  p <- ggVennDiagram(venn_list, label = "count", set_size = 0, edge_size = 0.35) +
    scale_fill_gradient(low = "white", high = high_col) + theme_void(base_family = FONT_FAMILY) +
    theme(text = element_text(family = FONT_FAMILY, size = FONT_PT), legend.position = "none",
          plot.margin = margin(0, 0, 0, 0))
  
  gb <- ggplotGrob(p)
  
  idx <- grep("setLabel", gb$layout$name)
  if (length(idx) > 0) {
    for (i in idx) {
      gb$grobs[[i]] <- grid::nullGrob()
    }
  }
  
  patchwork::wrap_elements(gb)
}

# Bottom factor labels
make_factor_label_block <- function(label_text) {
  ggplot() +
    annotate("text", x = 1, y = 1, label = label_text, hjust = 0.5, vjust = 0.5, size = 3.4, fontface = "bold") +
    xlim(0, 2) + ylim(0, 2) + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}

# Build one 2-row block for a given view
build_view_block <- function(view_name, matched_pairs_df, ref_model, scot_model, n_top = 20) {
  panel_store <- vector("list", length = length(signs) * nrow(matched_pairs_df))
  summary_list_local <- list()
  
  for (row_idx in seq_along(signs)) {
    sign_i <- signs[row_idx]
    
    for (col_idx in seq_len(nrow(matched_pairs_df))) {
      ref_factor  <- matched_pairs_df$ReferenceFactor[col_idx]
      scot_factor <- matched_pairs_df$SCOTFactor[col_idx]
      
      ref_feats <- get_top_features_signed(model = ref_model, view = view_name,
                                           factor_name = ref_factor, sign = sign_i, n_top = n_top)
      
      scot_feats <- get_top_features_signed(model = scot_model, view = view_name,
                                            factor_name = scot_factor, sign = sign_i, n_top = n_top)
      
      overlap_feats <- intersect(ref_feats, scot_feats)
      union_feats   <- union(ref_feats, scot_feats)
      jaccard <- if (length(union_feats) == 0) NA_real_ else length(overlap_feats) / length(union_feats)
      
      plot_idx <- (row_idx - 1) * nrow(matched_pairs_df) + col_idx
      
      panel_store[[plot_idx]] <- make_small_venn_plot(ref_features = ref_feats,
                                                      scot_features = scot_feats, sign = sign_i)
      
      summary_list_local[[length(summary_list_local) + 1]] <- data.frame(View = view_name,
                                                                         Sign = sign_i,
                                                                         ReferenceFactor = ref_factor,
                                                                         SCOTFactor = scot_factor,
                                                                         N_reference = length(ref_feats),
                                                                         N_scot = length(scot_feats),
                                                                         N_overlap = length(overlap_feats),
                                                                         Jaccard = jaccard,
                                                                         OverlapFeatures = paste(overlap_feats, collapse = "; "),
                                                                         stringsAsFactors = FALSE)
    }
  }
  
  row_pos <- panel_store[[1]] + panel_store[[2]] + panel_store[[3]] +
    plot_layout(ncol = 3)
  
  row_neg <- panel_store[[4]] + panel_store[[5]] + panel_store[[6]] +
    plot_layout(ncol = 3)
  
  block <- row_pos / row_neg +
    plot_layout(heights = c(1, 1))
  
  list(plot = block, summary = bind_rows(summary_list_local))
}

# Build transcriptomics and proteomics blocks
tx_block <- build_view_block(view_name = "Transcriptomics", matched_pairs_df = matched_pairs_2F,
                             ref_model = ref_model, scot_model = scot_model, n_top = n_top)

prot_block <- build_view_block(view_name = "Proteomics", matched_pairs_df = matched_pairs_2F,
                               ref_model = ref_model, scot_model = scot_model, n_top = n_top)

summary_df <- bind_rows(tx_block$summary, prot_block$summary)

# Bottom labels shown once only
bottom_labels <- make_factor_label_block("Matched factor 1") + make_factor_label_block("Matched factor 2") +
  make_factor_label_block("Matched factor 3") + plot_layout(ncol = 3)

# Final combined Venn plot
add_block_box <- function(p, line_size = 0.4) {
  cowplot::ggdraw() + cowplot::draw_plot(p, x = 0.015, y = 0.02, width = 0.97, height = 0.96) +
    theme(plot.background = element_rect(fill = NA, colour = "black", linewidth = line_size))
}

tx_boxed   <- add_block_box(tx_block$plot)
prot_boxed <- add_block_box(prot_block$plot)

make_F2F_global_legend <- function() {
  ggplot() +
    annotate("point", x = 0.05, y = 0.65, shape = 22, size = 4, fill = sign_colors[["positive"]],
             colour = "black", stroke = 0.3) +
    annotate("text", x = 0.10, y = 0.65, label = "Positively-weighted features", hjust = 0,
             vjust = 0.5, size = FONT_MM, family = FONT_FAMILY) +
    annotate("point", x = 0.05, y = 0.35, shape = 22, size = 4, fill = sign_colors[["negative"]],
             colour = "black", stroke = 0.3) +
    annotate("text", x = 0.10, y = 0.35, label = "Negatively-weighted features", hjust = 0,
             vjust = 0.5, size = FONT_MM, family = FONT_FAMILY) +
    annotate("text", x = 0.65, y = 0.65, label = "Upper circle: Reference paired model", hjust = 0,
             vjust = 0.5, size = FONT_MM, family = FONT_FAMILY) +
    annotate("text", x = 0.65, y = 0.35, label = "Lower circle: SCOT+-aligned model", hjust = 0,
             vjust = 0.5, size = FONT_MM, family = FONT_FAMILY) +
    coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1), clip = "off") +
    theme_void(base_family = FONT_FAMILY) +
    theme(plot.margin = margin(t = 2, r = 12, b = 2, l = 6))
}

F2F_legend <- make_F2F_global_legend()

F2F_plot <- (patchwork::wrap_elements(F2F_legend) / patchwork::wrap_elements(tx_boxed) /
                 patchwork::plot_spacer() / patchwork::wrap_elements(prot_boxed) /
                 patchwork::wrap_elements(bottom_labels)) +
  plot_layout(heights = c(0.24, 1, 0.03, 1, 0.12))

save_svg(filename = file.path(F2F_outdir, "Figure2F_combined_venn_signed.svg"), plot = F2F_plot, width = 5, height = 7.1)

write.csv(summary_df, file.path(F2F_outdir, "Figure2F_feature_overlap_venn_signed_summary.csv"), row.names = FALSE)

print(summary_df)