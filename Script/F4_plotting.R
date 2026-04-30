library(ggplot2)
library(dplyr)
library(tidyr)
library(MOFA2)
library(MOFAdata)
library(svglite)
library(cowplot)
library(patchwork)
library(clue)
library(scales)

# Paths
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
ref_model_path <- file.path(script_dir, "trained_model/Paired_model.hdf5")
unp_model_path <- file.path(script_dir, "trained_model/Unpaired_model.hdf5")
ref_rna_meta_path  <- file.path(script_dir, "processed_data/processed_transcriptomics_metadata.csv")
ref_prot_meta_path <- file.path(script_dir, "processed_data/processed_proteomics_metadata.csv")

outdir <- file.path(script_dir, "figure4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Global style
theme_set(theme_bw(base_size = 13, base_family = "sans"))

cellline_colors <- c("C10"  = "#000000", "SVEC" = "#D55E00")

modality_colors <- c("Transcriptomics" = "#E69F00", "Proteomics" = "#56B4E9")

sign_colors_5C <- c("Positive" = "#B40426", "Negative" = "#3B4CC0")

# Load models and gene sets
ref_model <- load_model(ref_model_path)
unp_model <- load_model(unp_model_path)

data("MSigDB_v6.0_C5_mouse")

# Helper: extract factor values for unpaired model
extract_unpaired_factor_df <- function(model, rna_meta_path, prot_meta_path) {
  
  factor_df <- get_factors(model, factors = "all", as.data.frame = TRUE)
  
  expected_cols <- c("sample", "factor", "value")
  if (!all(expected_cols %in% colnames(factor_df))) {
    stop("Factor table does not contain expected columns: sample, factor, value")
  }
  
  rna_meta <- read.csv(rna_meta_path, stringsAsFactors = FALSE)
  prot_meta <- read.csv(prot_meta_path, stringsAsFactors = FALSE)
  
  rna_meta <- rna_meta %>%
    dplyr::select(SampleID, CellLine) %>%
    distinct(SampleID, .keep_all = TRUE) %>%
    mutate(sample = paste0(SampleID, "_rna"), Modality = "Transcriptomics")
  
  prot_meta <- prot_meta %>%
    dplyr::select(SampleID, CellLine) %>%
    distinct(SampleID, .keep_all = TRUE) %>%
    mutate(sample = paste0(SampleID, "_prot"), Modality = "Proteomics")
  
  meta_df <- bind_rows(rna_meta, prot_meta)
  
  factor_df <- factor_df %>%
    left_join(meta_df, by = "sample") %>%
    mutate(Factor = factor(factor,
                           levels = paste0("Factor",
                                           seq_len(max(as.integer(gsub("Factor", "", unique(factor))))))))
  
  if (any(is.na(factor_df$CellLine))) {
    warning("Some unpaired samples were not matched to metadata.")
  }
  
  factor_df
}

# Helper: compute reference vs unpaired latent concordance
compute_unpaired_concordance <- function(m_ref, m_unp) {
  
  Z_ref <- get_factors(m_ref, factors = "all", as.data.frame = FALSE)[[1]]
  Z_unp <- get_factors(m_unp, factors = "all", as.data.frame = FALSE)[[1]]
  
  Z_unp_rna  <- Z_unp[grepl("_rna$", rownames(Z_unp)), , drop = FALSE]
  Z_unp_prot <- Z_unp[grepl("_prot$", rownames(Z_unp)), , drop = FALSE]
  
  rownames(Z_unp_rna)  <- sub("_rna$", "", rownames(Z_unp_rna))
  rownames(Z_unp_prot) <- sub("_prot$", "", rownames(Z_unp_prot))
  
  colnames(Z_ref)      <- paste0("Reference_F", seq_len(ncol(Z_ref)))
  colnames(Z_unp_rna)  <- paste0("Unpaired_F", seq_len(ncol(Z_unp_rna)))
  colnames(Z_unp_prot) <- paste0("Unpaired_F", seq_len(ncol(Z_unp_prot)))
  
  common_rna  <- intersect(rownames(Z_ref), rownames(Z_unp_rna))
  common_prot <- intersect(rownames(Z_ref), rownames(Z_unp_prot))
  
  Z_ref_rna <- Z_ref[match(common_rna, rownames(Z_ref)), , drop = FALSE]
  Z_unp_rna <- Z_unp_rna[match(common_rna, rownames(Z_unp_rna)), , drop = FALSE]
  
  Z_ref_prot <- Z_ref[match(common_prot, rownames(Z_ref)), , drop = FALSE]
  Z_unp_prot <- Z_unp_prot[match(common_prot, rownames(Z_unp_prot)), , drop = FALSE]
  
  cor_mat_rna <- cor(Z_ref_rna, Z_unp_rna, use = "pairwise.complete.obs", method = "pearson")
  cor_mat_prot <- cor(Z_ref_prot, Z_unp_prot, use = "pairwise.complete.obs", method = "pearson")
  
  combined_score <- (abs(cor_mat_rna) + abs(cor_mat_prot)) / 2
  
  cost_mat <- 1 - combined_score
  assignment <- solve_LSAP(cost_mat)
  
  matched_pairs <- data.frame(
    ReferenceFactor = rownames(combined_score),
    UnpairedFactor = colnames(combined_score)[assignment],
    MeanAbsCorrelation = combined_score[cbind(seq_len(nrow(combined_score)), assignment)],
    RNA_Correlation = cor_mat_rna[cbind(seq_len(nrow(cor_mat_rna)), assignment)],
    Proteomics_Correlation = cor_mat_prot[cbind(seq_len(nrow(cor_mat_prot)), assignment)],
    stringsAsFactors = FALSE)
  
  list(cor_mat_rna = cor_mat_rna, cor_mat_prot = cor_mat_prot, combined_score = combined_score,
       assignment = assignment, matched_pairs = matched_pairs)
}

# Helper: normalize feature names for enrichment
normalize_features <- function(model) {
  for (v in names(features_names(model))) {
    feats <- features_names(model)[[v]]
    feats <- toupper(feats)
    feats <- sub("_(RNA|PROT)$", "", feats)
    features_names(model)[[v]] <- feats
  }
  model
}

# Helper: enrichment utilities
get_sig_paths <- function(model, view, factor, sign_name) {
  enrich <- tryCatch(run_enrichment(model, view = view, factors = factor,
                                    feature.sets = MSigDB_v6.0_C5_mouse, sign = sign_name,
                                    statistical.test = "parametric"), error = function(e) NULL)
  
  if (is.null(enrich) || is.null(enrich$sigPathways)) return(character(0))
  unique(unlist(enrich$sigPathways))
}

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

# Compute concordance
conc <- compute_unpaired_concordance(ref_model, unp_model)

cor_mat_rna <- conc$cor_mat_rna
cor_mat_prot <- conc$cor_mat_prot
combined_score <- conc$combined_score
matched_pairs <- conc$matched_pairs

write.csv(matched_pairs, file.path(outdir, "matched_factor_pairs_reference_vs_unpaired.csv"),
          row.names = FALSE)

# Figure 4A data
factor_df_4A <- extract_unpaired_factor_df(model = unp_model, rna_meta_path = ref_rna_meta_path,
                                           prot_meta_path = ref_prot_meta_path)

factor_df_4A <- factor_df_4A %>%
  mutate(Modality = factor(Modality, levels = c("Transcriptomics", "Proteomics")),
         Factor   = factor(Factor, levels = c("Factor1", "Factor2", "Factor3", "Factor4")))

factor_pos_4A <- data.frame(Factor = factor(c("Factor1", "Factor2", "Factor3", "Factor4"),
                                            levels = c("Factor1", "Factor2", "Factor3", "Factor4")),
                            x_rna  = c(0.92, 1.42, 1.92, 2.42), x_prot = c(1.08, 1.58, 2.08, 2.58))

factor_df_4A <- factor_df_4A %>%
  dplyr::select(-dplyr::any_of(c("x_rna", "x_prot", "x"))) %>%
  left_join(factor_pos_4A, by = "Factor") %>%
  mutate(x = dplyr::case_when(Modality == "Transcriptomics" ~ .data$x_rna,
                              Modality == "Proteomics" ~ .data$x_prot)) %>%
  filter(!is.na(x))

yrange_4A <- range(factor_df_4A$value, na.rm = TRUE)
y_annot_4A <- yrange_4A[1] - 0.10 * diff(yrange_4A)

modality_colors <- c("Transcriptomics" = "#E69F00", "Proteomics" = "#56B4E9")

cellline_shapes <- c("C10" = 16, "SVEC" = 17)

p4A <- ggplot(factor_df_4A, aes(x = x, y = value, colour = Modality, shape = CellLine)) +
  geom_jitter(width = 0.03, height = 0, size = 1.4, alpha = 0.8) +
  scale_colour_manual(values = modality_colors) +
  scale_shape_manual(values = cellline_shapes) +
  scale_x_continuous(breaks = c(1.00, 1.50, 2.00, 2.50),
                     labels = c("Factor1", "Factor2", "Factor3", "Factor4")) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Factor value", colour = NULL, shape = NULL) +
  theme_bw(base_size = 13, base_family = "sans") +
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12, margin = margin(r = 10)), legend.title = element_blank(),
        legend.position = "right", legend.justification = "center",
        legend.background = element_blank(), legend.key = element_blank(),
        legend.text = element_text(size = 7), legend.margin = margin(6, 6, 6, 6),
        legend.box.margin = margin(2, 2, 2, 2), legend.spacing.y = unit(0.2, "cm"),
        legend.key.size = unit(0.45, "cm"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = margin(t = 12, r = 8, b = 18, l = 8))

ggsave(filename = file.path(outdir, "Figure4A_unpaired_factor_values.svg"),
       plot = p4A, device = svglite, width = 5.4, height = 4.2)

write.csv(factor_df_4A, file.path(outdir, "Figure5A_unpaired_factor_values.csv"), row.names = FALSE)

# Figure 4B data: heatmaps in ggplot style
ordered_unpaired <- matched_pairs$UnpairedFactor
cor_df_rna <- as.data.frame(as.table(cor_mat_rna[, ordered_unpaired, drop = FALSE]))
colnames(cor_df_rna) <- c("ReferenceFactor", "UnpairedFactor", "Correlation")
cor_df_rna$Modality <- "Transcriptomics"

cor_df_prot <- as.data.frame(as.table(cor_mat_prot[, ordered_unpaired, drop = FALSE]))
colnames(cor_df_prot) <- c("ReferenceFactor", "UnpairedFactor", "Correlation")
cor_df_prot$Modality <- "Proteomics"

# clean labels for plotting only
ref_levels <- rownames(cor_mat_rna)
unp_levels <- ordered_unpaired

ref_labels <- sub("^Reference_", "", ref_levels)
unp_labels <- sub("^Unpaired_", "", unp_levels)

cor_df_4B <- bind_rows(cor_df_rna, cor_df_prot) %>%
  mutate(Modality = factor(Modality, levels = c("Transcriptomics", "Proteomics")),
         ReferenceFactor = factor(ReferenceFactor, levels = rev(ref_levels),
                                  labels = rev(sub("^Reference_", "", ref_levels))),
         UnpairedFactor = factor(UnpairedFactor, levels = unp_levels,
                                 labels = sub("^Unpaired_", "", unp_levels)),
         label = sprintf("%.2f", Correlation))

p4B <- ggplot(cor_df_4B, aes(x = UnpairedFactor, y = ReferenceFactor, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  facet_wrap(~Modality, nrow = 1) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
                       limits = c(-1, 1), name = NULL) +
  labs(x = "Unpaired factor", y = "Reference factor") +
  theme_bw(base_size = 13, base_family = "sans") +
  theme(axis.text.x = element_text(size = 10,  hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, margin = margin(t = 10)),
        axis.title.y = element_text(size = 12, margin = margin(r = 10)),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 11), legend.position = "right",
        legend.background = element_blank(), legend.box.background = element_blank(),
        legend.key = element_blank(), legend.text = element_text(size = 8),
        legend.title = element_text(size = 9), panel.grid = element_blank(),
        plot.margin = margin(t = 12, r = 8, b = 8, l = 8))

ggsave(filename = file.path(outdir, "Figure4B_latent_concordance_heatmaps.svg"),
       plot = p4B, device = svglite, width = 7.6, height = 3.9)

write.csv(cor_df_4B, file.path(outdir, "Figure4B_latent_concordance_heatmaps.csv"),
          row.names = FALSE)

# Figure 4C data: pathway overlap for best-matched factor
m_ref_norm <- normalize_features(ref_model)
m_unp_norm <- normalize_features(unp_model)

top_pair <- matched_pairs %>%
  arrange(desc(MeanAbsCorrelation)) %>%
  slice(1)

ref_factor <- as.integer(sub("Reference_F", "", top_pair$ReferenceFactor))
unp_factor <- as.integer(sub("Unpaired_F", "", top_pair$UnpairedFactor))

# sign alignment
overall_sign <- sign(top_pair$RNA_Correlation + top_pair$Proteomics_Correlation)
if (overall_sign == 0) {
  overall_sign <- sign(top_pair$RNA_Correlation)
}
if (overall_sign == 0) {
  overall_sign <- 1
}

views <- c("Transcriptomics", "Proteomics")
signs <- c("positive", "negative")

res_list <- list()

for (v in views) {
  for (s in signs) {
    
    ref_paths <- get_sig_paths(m_ref_norm, v, ref_factor, s)
    
    # sign-align unpaired model if needed
    if (overall_sign > 0) {
      unp_sign <- s
    } else {
      unp_sign <- ifelse(s == "positive", "negative", "positive")
    }
    
    unp_paths <- get_sig_paths(m_unp_norm, v, unp_factor, unp_sign)
    
    res_list[[paste(v, s, sep = "_")]] <- data.frame(
      View = v,
      Sign = s,
      ReferenceFactor = top_pair$ReferenceFactor, UnpairedFactor = top_pair$UnpairedFactor,
      MeanAbsCorrelation = top_pair$MeanAbsCorrelation, RNA_Correlation = top_pair$RNA_Correlation,
      Proteomics_Correlation = top_pair$Proteomics_Correlation, n_ref = length(ref_paths),
      n_unpaired = length(unp_paths), n_intersection = length(intersect(ref_paths, unp_paths)),
      n_union = length(union(ref_paths, unp_paths)), Jaccard = jaccard(ref_paths, unp_paths),
      stringsAsFactors = FALSE)
  }
}

pathway_summary <- bind_rows(res_list)

write.csv(pathway_summary, file.path(outdir, "Figure4C_top_factor_pathway_comparison.csv"),
          row.names = FALSE)

pathway_df_4C <- pathway_summary %>%
  mutate(View = factor(View, levels = c("Transcriptomics", "Proteomics")),
         Sign = factor(Sign, levels = c("positive", "negative"), labels = c("Positive", "Negative")),
         Jaccard_plot = ifelse(is.na(Jaccard), 0, Jaccard),
         label = ifelse(is.na(Jaccard), "NA", sprintf("%.2f", Jaccard)), 
         alpha_group = ifelse(is.na(Jaccard), "NA", "Observed"))

p4C <- ggplot(pathway_df_4C, aes(x = Sign, y = Jaccard_plot, fill = Sign, alpha = alpha_group)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = label), vjust = -0.35, size = 3.2) +
  facet_wrap(~View, nrow = 1) +
  scale_fill_manual(values = sign_colors_5C) +
  scale_alpha_manual(values = c("Observed" = 1, "NA" = 0.35), guide = "none") +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, by = 0.2),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Jaccard similarity", fill = NULL) +
  theme_bw(base_size = 13, base_family = "sans") +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, margin = margin(r = 10)),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 11), legend.position = "none",
        legend.background = element_blank(), legend.box.background = element_blank(),
        legend.key = element_blank(), legend.text = element_text(size = 8),
        legend.justification = c(1, 1), legend.margin = margin(6, 6, 6, 6),
        legend.box.margin = margin(2, 2, 2, 2), panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), plot.margin = margin(t = 12, r = 8, b = 8, l = 8))

ggsave(filename = file.path(outdir, "Figure4C_top_factor_pathway_overlap.svg"),
       plot = p4C, device = svglite, width = 6.8, height = 3.8)

# Combined Figure 4
row2 <- p4B + p4C +
  plot_layout(widths = c(1.05, 0.85))

figure4 <- p4A / row2 +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "A")

top_row <- p4A + plot_spacer() +
  plot_layout(widths = c(6.4, 2.0))

row2 <- p4B + p4C +
  plot_layout(widths = c(1.15, 0.95))

figure4 <- top_row / row2 +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(filename = file.path(outdir, "Figure4_combined.svg"), plot = figure4, device = svglite,
       width = 8.4, height = 7.2)
