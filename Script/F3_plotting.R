library(MOFA2)
library(MOFAdata)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(readr)
library(svglite)
library(cowplot)
library(stringr)

set.seed(1)

# Paths
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
ref_model_path  <- file.path(script_dir, "trained_model/Paired_model.hdf5")
scot_model_path <- file.path(script_dir, "trained_model/SCOT+-aligned_model.hdf5")
fig2_dir <- file.path(script_dir, "figure2_S345")
outdir   <- file.path(script_dir, "figure3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

matched_pairs_file <- file.path(fig2_dir, "Figure2_matched_factor_pairs.csv")

# Global style
theme_set(theme_bw(base_size = 12, base_family = "sans"))

# Load models
m_ref  <- load_model(ref_model_path)
m_scot <- load_model(scot_model_path)

data("MSigDB_v6.0_C5_mouse")
gene_sets <- MSigDB_v6.0_C5_mouse

# Helper functions
standardize_factor_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  keep <- grepl("^Factor[0-9]+$", x)
  x[keep] <- x[keep]
  
  idx <- grepl("^(Reference_F|SCOT_F)[0-9]+$", x)
  x[idx] <- sub("^(Reference_F|SCOT_F)([0-9]+)$", "Factor\\2", x[idx])
  
  idx <- grepl("^F[0-9]+$", x)
  x[idx] <- sub("^F([0-9]+)$", "Factor\\1", x[idx])
  
  idx <- grepl("^factor[0-9]+$", x, ignore.case = TRUE)
  x[idx] <- sub("^factor([0-9]+)$", "Factor\\1", x[idx], ignore.case = TRUE)
  
  idx <- grepl("^[0-9]+$", x)
  x[idx] <- paste0("Factor", x[idx])
  
  x
}

normalize_features <- function(model) {
  for (v in names(features_names(model))) {
    feats <- features_names(model)[[v]]
    feats <- toupper(feats)
    feats <- sub("_(RNA|PROT)$", "", feats)
    features_names(model)[[v]] <- feats
  }
  model
}

get_sig_pathways <- function(enrich_obj) {
  if (is.null(enrich_obj)) return(character(0))
  sig_list <- enrich_obj$sigPathways
  if (is.null(sig_list)) return(character(0))
  unique(unlist(sig_list))
}

jaccard_similarity <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x | y)
  if (union == 0) return(NA_real_)
  intersection / union
}

build_dotplot_df <- function(enrichment, sign_label, factors_vec, model_label, top_n = 10) {
  
  df_list <- list()
  
  for (fac in factors_vec) {
    
    sig_paths <- enrichment$sigPathways[[fac]]
    if (is.null(sig_paths) || length(sig_paths) == 0) next
    
    scores_all <- enrichment$set.statistics[, fac]
    sig_paths_valid <- intersect(sig_paths, names(scores_all))
    if (length(sig_paths_valid) == 0) next
    
    scores <- scores_all[sig_paths_valid]
    top_paths <- names(sort(abs(scores), decreasing = TRUE))[seq_len(min(top_n, length(scores)))]
    if (length(top_paths) == 0) next
    
    valid_gene_paths <- intersect(top_paths, rownames(enrichment$feature.sets))
    if (length(valid_gene_paths) == 0) next
    
    gene_counts <- rowSums(enrichment$feature.sets[valid_gene_paths, , drop = FALSE])
    scores_subset <- scores[valid_gene_paths]
    
    score_signed <- if (sign_label == "negative") {
      -as.numeric(scores_subset)
    } else {
      as.numeric(scores_subset)
    }
    
    df_tmp <- data.frame(Model = model_label, Factor = paste0("Factor", fac), FactorNumber = fac,
                         Pathway = valid_gene_paths, Score = score_signed,
                         Genes = as.numeric(gene_counts), Sign = sign_label, stringsAsFactors = FALSE)
    
    df_list[[length(df_list) + 1]] <- df_tmp
  }
  
  bind_rows(df_list)
}

run_model_enrichment_all <- function(model, model_label, views = c("Transcriptomics", "Proteomics"),
                                     signs = c("positive", "negative"), max_pathways = 10) {
  
  K <- get_dimensions(model)$K
  factors_all <- seq_len(K)
  
  out_list <- list()
  
  for (view in views) {
    message("Running enrichment: ", model_label, " / ", view)
    
    enrich_pos <- run_enrichment(model, view = view, factors = factors_all, feature.sets = gene_sets,
                                 sign = "positive", statistical.test = "parametric")
    
    enrich_neg <- run_enrichment(model, view = view, factors = factors_all, feature.sets = gene_sets,
                                 sign = "negative", statistical.test = "parametric")
    
    df_pos <- build_dotplot_df(enrich_pos, "positive", factors_all, model_label, top_n = max_pathways)
    df_neg <- build_dotplot_df(enrich_neg, "negative", factors_all, model_label, top_n = max_pathways)
    
    out_list[[view]] <- bind_rows(df_pos, df_neg)
  }
  
  out_list
}

make_combined_dotplot <- function(df_view) {
  
  if (nrow(df_view) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "No significant pathways", size = 5) +
             theme_void())
  }
  
  df_view <- df_view %>%
    mutate(Model = factor(Model, levels = c("Reference", "SCOT+-aligned")),
           Factor = factor(Factor, levels = paste0("Factor", sort(unique(FactorNumber)))),
           ViewModel = ifelse(Model == "Reference", "Ref", "SCOT+"))
  
  # Order pathways by overall max absolute score
  pathway_order <- df_view %>%
    group_by(Pathway) %>%
    summarise(order_stat = max(abs(Score), na.rm = TRUE), .groups = "drop") %>%
    arrange(order_stat) %>%
    pull(Pathway)
  
  df_view <- df_view %>%
    mutate(Pathway = factor(Pathway, levels = pathway_order))
  
  # Manual x positions: Ref and SCOT+ side-by-side within each factor
  factor_pos <- data.frame(Factor = factor(paste0("Factor", sort(unique(df_view$FactorNumber))),
                                           levels = paste0("Factor", sort(unique(df_view$FactorNumber)))))
  
  factor_ids <- levels(factor_pos$Factor)
  n_factors <- length(factor_ids)
  
  factor_pos$x_ref  <- seq(0.9, by = 0.6, length.out = n_factors)
  factor_pos$x_scot <- seq(1.1, by = 0.6, length.out = n_factors)
  
  df_view <- df_view %>%
    left_join(factor_pos, by = "Factor") %>%
    mutate(x = case_when(Model == "Reference" ~ x_ref, Model == "SCOT+-aligned" ~ x_scot)) %>%
    filter(!is.na(x))
  
  # central factor tick positions
  x_breaks <- as.vector(rbind(factor_pos$x_ref, factor_pos$x_scot))
  x_breaks <- x_breaks[!is.na(x_breaks)]
  
  x_labels <- unlist(lapply(seq_len(n_factors), function(i) {
    labs <- c(paste0("Ref_F", i), paste0("SCOT+_F", i))
    if (is.na(factor_pos$x_scot[i])) labs <- labs[1]
    labs
  }))
  
  p <- ggplot(df_view, aes(x = x, y = Pathway, size = Genes, color = Score)) +
    geom_point(alpha = 0.9) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          name = "Signed\nenrichment\nscore", guide = guide_colorbar(order = 1)) +
    scale_size(range = c(2, 6), name = "# Genes", guide = guide_legend(order = 2)) +
    labs(x = NULL, y = "Enriched pathway") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 8), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), legend.title = element_text(size = 9, face = "bold"),
          legend.text = element_text(size = 8), plot.margin = margin(t = 8, r = 8, b = 18, l = 8))
  
  return(p)
}

# Normalise feature names in copies
m_ref_gsea  <- normalize_features(m_ref)
m_scot_gsea <- normalize_features(m_scot)

# Run enrichment
max_pathways <- 10

ref_enrich_list <- run_model_enrichment_all(model = m_ref_gsea, model_label = "Reference",
                                            max_pathways = max_pathways)

scot_enrich_list <- run_model_enrichment_all(model = m_scot_gsea, model_label = "SCOT+-aligned",
                                             max_pathways = max_pathways)

# Combine by view
df_rna <- bind_rows(ref_enrich_list[["Transcriptomics"]], scot_enrich_list[["Transcriptomics"]])

df_prot <- bind_rows(ref_enrich_list[["Proteomics"]], scot_enrich_list[["Proteomics"]])

# Save underlying data
write.csv(df_rna,  file.path(outdir, "Figure3A_Transcriptomics_combined_dotplot_data.csv"), row.names = FALSE)
write.csv(df_prot, file.path(outdir, "Figure3B_Proteomics_combined_dotplot_data.csv"), row.names = FALSE)

# Figure 3A and 3B dot plots
p3A <- make_combined_dotplot(df_rna)
p3B <- make_combined_dotplot(df_prot)

# Pathway concordance heatmap data
if (!file.exists(matched_pairs_file)) {
  stop("Matched factor pair file not found: ", matched_pairs_file)
}

matched_pairs <- read.csv(matched_pairs_file, stringsAsFactors = FALSE)

matched_pairs <- matched_pairs %>%
  rename(factor_ref = ReferenceFactor, factor_scot = SCOTFactor, factor_correlation = Correlation) %>%
  mutate(factor_ref = standardize_factor_name(factor_ref),
         factor_scot = standardize_factor_name(factor_scot))

views <- c("Transcriptomics", "Proteomics")
signs <- c("positive", "negative")

# Run enrichment factor-by-factor for concordance
run_model_enrichment_factorwise <- function(model, model_name) {
  enrichment_list <- list()
  K <- get_dimensions(model)$K
  
  for (view in views) {
    for (sign in signs) {
      for (f in seq_len(K)) {
        enrich_res <- tryCatch(
          {
            run_enrichment(model, view = view, factors = f, feature.sets = gene_sets, sign = sign,
                           statistical.test = "parametric")
          }, error = function(e) NULL)
        
        enrichment_list[[paste(model_name, view, sign, paste0("Factor", f), sep = "_")]] <- enrich_res
      }
    }
  }
  
  enrichment_list
}

flatten_enrichment <- function(enrich_list) {
  df_list <- lapply(names(enrich_list), function(nm) {
    obj <- enrich_list[[nm]]
    
    parts <- strsplit(nm, "_")[[1]]
    model_name <- parts[1]
    view <- parts[2]
    sign <- parts[3]
    factor <- parts[4]
    
    paths <- get_sig_pathways(obj)
    if (length(paths) == 0) return(NULL)
    
    data.frame(model = model_name, view = view, sign = sign, factor = standardize_factor_name(factor),
               pathway = paths, stringsAsFactors = FALSE)
  })
  
  bind_rows(df_list)
}

message("Running factor-wise enrichment for pathway concordance heatmap...")
enrich_ref_factorwise  <- run_model_enrichment_factorwise(m_ref_gsea, "Reference")
enrich_scot_factorwise <- run_model_enrichment_factorwise(m_scot_gsea, "SCOTplus")

df_ref_path  <- flatten_enrichment(enrich_ref_factorwise)
df_scot_path <- flatten_enrichment(enrich_scot_factorwise)

matched_pathway_concordance <- expand.grid(view = views, sign = signs, stringsAsFactors = FALSE) %>%
  tidyr::crossing(matched_pairs) %>%
  rowwise() %>%
  mutate(scot_sign_aligned = case_when(factor_correlation < 0 & sign == "positive" ~ "negative",
                                       factor_correlation < 0 & sign == "negative" ~ "positive",
      TRUE ~ sign),
      ref_paths = list({
        v <- view
        s <- sign
        fref <- factor_ref
        x <- df_ref_path %>% 
          filter(.data$view == v, .data$sign == s, .data$factor == fref) %>%
          pull(pathway)
        unique(x)
    }), scot_paths = list({
      v <- view
      s_aligned <- scot_sign_aligned
      fscot <- factor_scot
      x <- df_scot_path %>%
        filter(.data$view == v, .data$sign == s_aligned, .data$factor == fscot) %>%
        pull(pathway)
      unique(x)
      }), n_ref = length(unlist(ref_paths)), n_scot = length(unlist(scot_paths)),
    n_intersection = length(intersect(unlist(ref_paths), unlist(scot_paths))),
    n_union = length(union(unlist(ref_paths), unlist(scot_paths))),
    jaccard = ifelse(n_union == 0, NA_real_, n_intersection / n_union)) %>%
  ungroup()

matched_pathway_concordance_export <- matched_pathway_concordance %>%
  mutate(ref_paths = vapply(ref_paths, function(x) paste(x, collapse = "; "), character(1)),
         scot_paths = vapply(scot_paths, function(x) paste(x, collapse = "; "), character(1)))

write.csv(matched_pathway_concordance_export,
          file.path(outdir, "Figure3C_pathway_concordance_matched_factors.csv"), row.names = FALSE)

# Figure 3C heatmap

df_heat <- matched_pathway_concordance %>%
  mutate(Factor = factor_ref, RowLabel = case_when(
    view == "Proteomics" & sign == "positive"      ~ "Proteomics (+)",
    view == "Proteomics" & sign == "negative"      ~ "Proteomics (-)",
    view == "Transcriptomics" & sign == "positive" ~ "Transcriptomics (+)",
    view == "Transcriptomics" & sign == "negative" ~ "Transcriptomics (-)",
    TRUE ~ paste(view, sign))) %>%
  dplyr::select(RowLabel, Factor, jaccard) %>%
  distinct()

df_heat$RowLabel <- factor(df_heat$RowLabel, levels = c("Proteomics (+)", "Proteomics (-)",
                                                        "Transcriptomics (+)", "Transcriptomics (-)"))

factor_levels <- matched_pairs$factor_ref
df_heat$Factor <- factor(df_heat$Factor, levels = factor_levels)

df_heat <- df_heat %>%
  mutate(label = ifelse(is.na(jaccard), "NA", sprintf("%.2f", jaccard)))

p3C <- ggplot(df_heat, aes(x = Factor, y = RowLabel, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), na.value = "grey85",
                      name = "Jaccard\nsimilarity") +
  labs(x = "Matched factor", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8))

p3C_shift <- cowplot::ggdraw() +
  cowplot::draw_plot(p3C, x = -0.04, y = 0, width = 1.04, height = 1)

# Combine into Figure 3
blank_spacer <- ggplot() + theme_void()

fig3_left <- cowplot::plot_grid(p3A, p3B, labels = c("A", "B"), label_size = 14,
                                label_fontface = "bold", label_x = 0.01, label_y = 0.99, ncol = 1,
                                rel_heights = c(1, 1))

fig3_right <- cowplot::plot_grid(p3C, blank_spacer, labels = c("C", ""), label_size = 14,
                                 label_fontface = "bold", label_x = c(0.01, 0.01), 
                                 label_y = c(0.99, 0.99), ncol = 1, rel_heights = c(0.45, 1.55))

fig3 <- cowplot::plot_grid(fig3_left, fig3_right, ncol = 2, rel_widths = c(1.9, 1.0), align = "h")

# Export
ggsave(filename = file.path(outdir, "Figure3_combined_pathway_recovery.svg"), plot = fig3,
       device = svglite, width = 15, height = 20)

# Optional: also save individual panels
ggsave(file.path(outdir, "Figure3A_transcriptomics_combined_dotplot.svg"),
       p3A, device = svglite, width = 9, height = 10)

ggsave(file.path(outdir, "Figure3B_proteomics_combined_dotplot.svg"),
       p3B, device = svglite, width = 11, height = 10)

ggsave(file.path(outdir, "Figure3C_pathway_concordance_heatmap.svg"),
       p3C, device = svglite, width = 7, height = 4.5)
