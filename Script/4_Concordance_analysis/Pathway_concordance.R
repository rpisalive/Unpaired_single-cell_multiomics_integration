library(MOFA2)
library(dplyr)
library(MOFAdata)
library(ggplot2)
library(tidyr)
library(purrr)
library(pheatmap)
library(tibble)
library(readr)

# Output ptah
outdir <- "C:/Users/49152/Downloads/Multi-omics/sc_vs_scot+/MOFA_complied/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

m_single <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/single_cell_trained_model_MOFAcomplied.hdf5")
m_scot   <- load_model("C:/Users/49152/Downloads/Multi-omics/SCOT_plus/MOFA_output/SCOTplus_trained_model_MOFAcomp.hdf5")

data("MSigDB_v6.0_C5_mouse")
gene_sets <- MSigDB_v6.0_C5_mouse

matched_pairs_file <- file.path(outdir, "matched_factor_pairs.csv")

#1 Helper functions
standardize_factor_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  # Already standardised
  keep <- grepl("^Factor[0-9]+$", x)
  x[keep] <- x[keep]
  
  # Reference_F1 / SCOT_F1 -> Factor1
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

shorten_factor_names <- function(fnames) {
  fnames <- gsub("^Single", "Ref", fnames)
  fnames <- gsub("^SCOTplus", "SCOT", fnames)
  fnames <- gsub("Transcriptomics", "RNA", fnames)
  fnames <- gsub("Proteomics", "PROT", fnames)
  fnames <- gsub("positive", "+", fnames)
  fnames <- gsub("negative", "-", fnames)
  fnames <- gsub("_([+-])_(Factor[0-9]+)", "-\\2(\\1)", fnames)
  fnames <- gsub("_", "-", fnames)
  fnames
}

#2 Matched factor pairs loading
if (!file.exists(matched_pairs_file)) {
  stop("matched_factor_pairs.csv not found. Run latent concordance analysis first.")
}

matched_pairs <- read.csv(matched_pairs_file, stringsAsFactors = FALSE)

if (all(c("ReferenceFactor", "SCOTFactor", "Correlation") %in% colnames(matched_pairs))) {
  matched_pairs <- matched_pairs %>%
    rename(factor_ref = ReferenceFactor, factor_scot = SCOTFactor,
           factor_correlation = Correlation)
} else if (!all(c("factor_ref", "factor_scot", "factor_correlation") %in% colnames(matched_pairs))) {
  stop("matched_factor_pairs.csv does not contain expected columns.")
}

matched_pairs <- matched_pairs %>%
  mutate(factor_ref = standardize_factor_name(factor_ref),
         factor_scot = standardize_factor_name(factor_scot))

print(matched_pairs)

#3 Normalise feature names in model copies
m_single_gsea <- normalize_features(m_single)
m_scot_gsea   <- normalize_features(m_scot)

#4 Settings
views <- c("Transcriptomics", "Proteomics")
signs <- c("positive", "negative")
K_single <- get_dimensions(m_single_gsea)$K
K_scot   <- get_dimensions(m_scot_gsea)$K

#5 Run enrichment for all factors, views and signs
run_model_enrichment <- function(model, model_name) {
  enrichment_list <- list()
  K <- get_dimensions(model)$K
  
  for (view in views) {
    for (sign in signs) {
      for (f in seq_len(K)) {
        enrich_res <- tryCatch(
          {
            run_enrichment(model, view = view, factors = f, feature.sets = gene_sets,
                           sign = sign, statistical.test = "parametric")
          },
          error = function(e) NULL
        )
        
        enrichment_list[[paste(model_name, view, sign, paste0("Factor", f), sep = "_")]] <- enrich_res
      }
    }
  }
  
  enrichment_list
}

enrich_single <- run_model_enrichment(m_single_gsea, "Single")
enrich_scot   <- run_model_enrichment(m_scot_gsea, "SCOTplus")

#6 Flatten enrichment results
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
    
    data.frame(model  = model_name, view   = view, sign   = sign,
               factor = standardize_factor_name(factor), pathway = paths,
               stringsAsFactors = FALSE)
  })
  
  df <- bind_rows(df_list)
  if (nrow(df) == 0) return(df)
  rownames(df) <- NULL
  df
}

df_single <- flatten_enrichment(enrich_single)
df_scot   <- flatten_enrichment(enrich_scot)

if (nrow(df_single) == 0 || nrow(df_scot) == 0) {
  stop("No significant pathways found in one or both models.")
}

df_combined <- bind_rows(df_single, df_scot) %>%
  mutate(factor_id = paste(model, view, sign, factor, sep = "_"))

#7 All-vs-all pathway presence matrix
pathway_mat <- df_combined %>%
  select(factor_id, pathway) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = pathway, values_from = value, values_fill = 0) %>%
  column_to_rownames("factor_id")

factor_ids <- rownames(pathway_mat)

jaccard_mat <- matrix(0, nrow = length(factor_ids), ncol = length(factor_ids),
                      dimnames = list(factor_ids, factor_ids))

for (i in seq_along(factor_ids)) {
  for (j in seq_along(factor_ids)) {
    jaccard_mat[i, j] <- jaccard_similarity(pathway_mat[i, ], pathway_mat[j, ])
  }
}

cross_mat <- jaccard_mat[
  grep("^Single_", rownames(jaccard_mat)),
  grep("^SCOTplus_", colnames(jaccard_mat)),
  drop = FALSE
]

rownames(cross_mat) <- shorten_factor_names(rownames(cross_mat))
colnames(cross_mat) <- shorten_factor_names(colnames(cross_mat))

write.csv(cross_mat, file.path(outdir, "pathway_concordance_all_vs_all.csv"),
          row.names = TRUE)

svg(file.path(outdir, "pathway_concordance_heatmap.svg"), width = 10, height = 10)
pheatmap(cross_mat, cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(c("white", "red"))(50),
         main = paste0("Pathway concordance (Jaccard similarity)\n",
                       "Ref = Reference, SCOT = SCOT+ aligned\n",
                       "RNA = Transcriptomic view, PROT = Proteomic view"), fontsize = 14,
         fontsize_row = 12, fontsize_col = 12, angle_col = 45, border_color = NA,
         legend = TRUE, legend_breaks = c(0, 0.5, 1), legend_labels = c("0", "0.5", "1"))
dev.off()

#8 Matched-factor pathway concordance
# Compare reference vs SCOT only for matched latent factors,
matched_pathway_concordance <- expand.grid(view = views, sign = signs,
                                           stringsAsFactors = FALSE) %>%
  tidyr::crossing(matched_pairs) %>%
  rowwise() %>%
  mutate(ref_paths = list({
      v <- view; s <- sign; fref <- factor_ref
      x <- df_single %>%
        filter(.data$view == v, .data$sign == s, .data$factor == fref) %>%
        pull(pathway)
      unique(x)
    }),
    scot_paths = list({
      v <- view; s <- sign; fscot <- factor_scot
      x <- df_scot %>%
        filter(.data$view == v, .data$sign == s, .data$factor == fscot) %>%
        pull(pathway)
      unique(x)
    }),
    n_ref = length(unlist(ref_paths)),
    n_scot = length(unlist(scot_paths)),
    n_intersection = length(intersect(unlist(ref_paths), unlist(scot_paths))),
    n_union = length(union(unlist(ref_paths), unlist(scot_paths))),
    jaccard = ifelse(n_union == 0, NA_real_, n_intersection / n_union)
  ) %>%
  ungroup() %>%
  select(view, sign, factor_ref, factor_scot, factor_correlation, n_ref, n_scot,
         n_intersection, n_union, jaccard)

write.csv(matched_pathway_concordance,
          file.path(outdir, "pathway_concordance_matched_factors.csv"),
          row.names = FALSE)

print(matched_pathway_concordance)

#9 Matched-factor heatmap
df_plot <- matched_pathway_concordance %>%
  mutate(
    Factor = factor_ref,
    RowLabel = case_when(view == "Proteomics" & sign == "positive" ~ "Proteomics (+)",
                         view == "Proteomics" & sign == "negative" ~ "Proteomics (-)",
                         view == "Transcriptomics" & sign == "positive" ~ "Transcriptomics (+)",
                         view == "Transcriptomics" & sign == "negative" ~ "Transcriptomics (-)",
                         TRUE ~ paste(view, sign))) %>%
  select(RowLabel, Factor, jaccard, n_ref, n_scot, n_intersection, n_union)

# Preserve row order
df_plot$RowLabel <- factor(df_plot$RowLabel, levels = c("Proteomics (+)",
                                                        "Proteomics (-)",
                                                        "Transcriptomics (+)",
                                                        "Transcriptomics (-)"))

# Preserve factor order
df_plot$Factor <- factor(df_plot$Factor, levels = unique(df_plot$Factor))

# Label NAs clearly for plotting
df_plot <- df_plot %>%
  mutate(label = ifelse(is.na(jaccard), "NA", sprintf("%.2f", jaccard)),
         jaccard_plot = jaccard)

p <- ggplot(df_plot, aes(x = Factor, y = RowLabel, fill = jaccard_plot)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1),
                      na.value = "grey85", name = "Jaccard\nsimilarity") +
  theme_bw(base_size = 12) +
  labs(title = "Matched-factor pathway concordance between reference and SCOT+ models",
       x = "Matched factor", y = NULL) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5))

print(p)

ggsave(filename = file.path(outdir, "Figure_PathwayConcordance_matched_clean.svg"),
       plot = p, device = "svg", width = 8, height = 4.5)

#10 Summary statistics
cat("\nMean matched-factor pathway Jaccard:", round(mean(matched_pathway_concordance$jaccard, na.rm = TRUE), 3), "\n")

for (v in views) {
  for (s in signs) {
    df_sub <- matched_pathway_concordance %>%
      filter(view == v, sign == s)
    
    cat("Mean pathway Jaccard for", v, "/", s, ":",
        round(mean(df_sub$jaccard, na.rm = TRUE), 3), "\n")
  }
}