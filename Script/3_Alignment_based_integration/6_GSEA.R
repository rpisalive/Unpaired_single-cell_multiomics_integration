library(data.table)
library(ggplot2)
library(MOFAdata)
library(dplyr)
library(MOFA2)
library(tidyr)
library(svglite)
library(extrafont)
library(readr)

set.seed(1)

#1 Paths and global settings
setwd("C:/Users/49152/Downloads/Multi-omics/")

outdir <- "C:/Users/49152/Downloads/Multi-omics/SCOT_plus/MOFA_output/graphs/MOFA+complied_SCOTplus/GSEA/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

loadfonts(device = "win")
theme_set(theme_bw(base_size = 10, base_family = "Arial"))

# Parameters
top_genes_n  <- 3
max_pathways <- 10
max_genes    <- 8

#2 Model loading
model <- load_model("SCOT_plus/MOFA_output/SCOTplus_trained_model_MOFAcomp.hdf5")

#3 Metadata loading and alignment
# For SCOT+ workflow, use the aligned metadata exported after SCOT+,
aligned_metadata <- read.csv("SCOT_plus/AGW_results/aligned_metadata.csv",
                             stringsAsFactors = FALSE, check.names = FALSE)

if (colnames(aligned_metadata)[1] %in% c("", "X", "Unnamed: 0")) {
  colnames(aligned_metadata)[1] <- "SampleID"
}

model_samples <- unname(unlist(samples_names(model)))

metadata <- aligned_metadata %>%
  filter(SampleID %in% model_samples) %>%
  distinct(SampleID, .keep_all = TRUE)

# Align order to model sample order
metadata <- metadata[match(model_samples, metadata$SampleID), , drop = FALSE]

# Safety checks
if (nrow(metadata) != length(model_samples)) {
  stop("Mismatch between model samples and aligned metadata rows after alignment.")
}
if (!all(metadata$SampleID == model_samples)) {
  stop("Sample order mismatch between model and aligned metadata.")
}

metadata_for_join <- metadata %>%
  rename(sample = SampleID)

model@samples_metadata <- model@samples_metadata %>%
  left_join(metadata_for_join, by = "sample")

# Inspect
print(head(model@samples_metadata))
print(colnames(model@samples_metadata))

#4 Gene set data
data("MSigDB_v6.0_C5_mouse")

# Settings
views <- c("Transcriptomics", "Proteomics")
factors_all <- seq_len(get_dimensions(model)$K)

#5 Helper functions
empty_plot <- function(label) {
  ggplot() +
    annotate("text", 0.5, 0.5, label = label, size = 5, hjust = 0.5, vjust = 0.5) +
    theme_void()
}

clean_feature_names <- function(x) {
  sub("_(rna|prot)$", "", x)
}

build_dotplot_df <- function(enrichment, sign_label, factors_vec, top_n = 10) {
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
    
    top_paths <- valid_gene_paths
    scores_subset <- scores[top_paths]
    
    score_signed <- if (sign_label == "negative") -as.numeric(scores_subset) else as.numeric(scores_subset)
    
    df_tmp <- data.frame(Factor = paste0("Factor", fac), FactorNumber = fac,
                         Pathway = top_paths, Score = score_signed,
                         Genes = as.numeric(gene_counts), Sign = sign_label,
                         stringsAsFactors = FALSE)
    
    df_list[[length(df_list) + 1]] <- df_tmp
  }
  
  bind_rows(df_list)
}

extract_factor_weights <- function(model_obj, view, factor) {
  w <- get_weights(model_obj, views = view, factors = factor, as.data.frame = TRUE)
  
  if (!all(c("feature", "value") %in% colnames(w))) {
    stop("Expected columns 'feature' and 'value' not found in weights table.")
  }
  
  w %>%
    mutate(feature_clean = clean_feature_names(feature), abs_value = abs(value),
           sign = case_when(value > 0 ~ "positive", value < 0 ~ "negative",
                            TRUE ~ "zero")) %>%
    arrange(desc(abs_value))
}

# Preparation of a separate model copy for GSEA
model_gsea <- model

for (v in views) {
  feat <- features_names(model_gsea)[[v]]
  
  # Remove modality suffix and convert to uppercase
  feat_clean <- toupper(sub("_(rna|prot)$", "", feat))
  
  if (any(duplicated(feat_clean))) {
    warning(
      "Duplicated feature names detected after cleaning in view: ", v,
      ". This may affect feature-set matching."
    )
  }
  
  features_names(model_gsea)[[v]] <- feat_clean
}

for (v in views) {
  model_feats <- features_names(model_gsea)[[v]]
  gs_feats <- colnames(MSigDB_v6.0_C5_mouse)
  
  cat("\nView:", v, "\n")
  cat("Model features:", length(model_feats), "\n")
  cat("Gene-set features:", length(gs_feats), "\n")
  cat("Overlap:", sum(model_feats %in% gs_feats), "\n")
}

#7 Main loop over views
for (view in views) {
  
  message("Processing view: ", view)
  
  view_outdir <- file.path(outdir, view)
  dir.create(view_outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Run enrichment
  enrich_pos <- run_enrichment(model_gsea, view = view, factors = factors_all,
                               feature.sets = MSigDB_v6.0_C5_mouse, sign = "positive",
                               statistical.test = "parametric")
  
  enrich_neg <- run_enrichment(model_gsea, view = view, factors = factors_all,
                               feature.sets = MSigDB_v6.0_C5_mouse, sign = "negative",
                               statistical.test = "parametric")
  
  # Combined dotplot Data
  df_pos <- build_dotplot_df(enrich_pos, "positive", factors_all, max_pathways)
  df_neg <- build_dotplot_df(enrich_neg, "negative", factors_all, max_pathways)
  df_dot <- bind_rows(df_pos, df_neg)
  
  if (nrow(df_dot) == 0) {
    message("  No significant pathways for ", view)
    next
  }
  
  df_dot <- df_dot %>%
    arrange(Score) %>%
    mutate(Pathway = factor(Pathway, levels = unique(Pathway)))
  
  p_dot <- ggplot(df_dot, aes(x = Factor, y = Pathway, size = Genes, color = Score)) +
    geom_point() +
    scale_x_discrete(drop = FALSE) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          guide = guide_colorbar(order = 1)) +
    scale_size(range = c(2, 6), guide = guide_legend(order = 2)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 8, face = "bold"),
          axis.title.y = element_text(size = 8, face = "bold"),
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 8)) +
    labs(x = "MOFA+ Factor", y = "Enriched Pathway", color = "Signed enrichment score",
         size = "# Genes")
  
  dot_file <- file.path(view_outdir, paste0(view, "_combined_dotplot.svg"))
  ggsave(dot_file, p_dot, device = svglite, width = 10, height = 9)
  
  fwrite(df_dot, file.path(view_outdir, paste0(view, "_combined_dotdf.csv")))
  message("  Saved combined dot plot: ", dot_file)

  # Precompute weights once per factor
  weights_by_factor <- lapply(factors_all, function(fac) extract_factor_weights(model, view, fac))
  names(weights_by_factor) <- paste0("Factor", factors_all)
  
  # Per-sign/ per-factor detailed ouput
  for (sign in c("positive", "negative")) {
    
    enrichment <- if (sign == "positive") enrich_pos else enrich_neg
    
    valid_factors <- factors_all[
      sapply(factors_all, function(f) {
        w <- weights_by_factor[[paste0("Factor", f)]]
        if (nrow(w) == 0) return(FALSE)
        if (sign == "positive") {
          sum(w$value > 0, na.rm = TRUE) > 5
        } else {
          sum(w$value < 0, na.rm = TRUE) > 5
        }
      })
    ]
    
    for (fac in valid_factors) {
      
      fac_dir <- file.path(view_outdir, paste0(sign, "_Factor", fac))
      dir.create(fac_dir, recursive = TRUE, showWarnings = FALSE)
      
      sig_pathways <- enrichment$sigPathways[[fac]]
      
      # Top enrichment plot
      p_top <- if (is.null(sig_pathways) || length(sig_pathways) == 0) {
        empty_plot("No significant pathways")
      } else {
        plot_enrichment(enrichment, factor = fac, max.pathways = max_pathways,
                        text_size = 0.7)
      }
      
      ggsave(file.path(fac_dir, "top_enrichment.svg"), p_top, device = svglite,
             width = 6, height = 5)
      
      # Detailed enrichment plot
      p_det <- if (is.null(sig_pathways) || length(sig_pathways) == 0) {
        empty_plot("No significant pathways")
      } else {
        plot_enrichment_detailed(enrichment, factor = fac, max.genes = max_genes,
                                 max.pathways = 5)
      }
      
      ggsave(file.path(fac_dir, "detailed_enrichment.svg"), p_det, device = svglite,
             width = 7, height = 6)
      
      # Export weights table for this factor
      w <- weights_by_factor[[paste0("Factor", fac)]]
      
      w_sign <- if (sign == "positive") {
        w %>% filter(value > 0)
      } else {
        w %>% filter(value < 0)
      }
      
      write.csv(w_sign, file = file.path(fac_dir, paste0(view, "_Factor", fac, "_", sign, "_weights.csv")),
                row.names = FALSE, quote = FALSE)
      
      # Export enrichment tables
      if (!is.null(sig_pathways) && length(sig_pathways) > 0) {
        
        scores_all <- enrichment$set.statistics[, fac]
        sig_paths_valid <- intersect(sig_pathways, names(scores_all))
        
        if (length(sig_paths_valid) > 0) {
          enrich_tbl <- data.frame(Pathway = sig_paths_valid,
                                   Score = as.numeric(scores_all[sig_paths_valid]),
                                   Sign = sign, Factor = fac, View = view,
                                   stringsAsFactors = FALSE) %>% 
            arrange(desc(abs(Score)))
          
          write.csv(enrich_tbl, file = file.path(fac_dir, paste0(view, "_Factor", fac, "_", sign, "_enrichment_table.csv")),
                    row.names = FALSE, quote = FALSE)
        }
      }
    }
  }
}