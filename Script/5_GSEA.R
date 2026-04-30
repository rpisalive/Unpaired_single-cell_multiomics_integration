library(data.table)
library(ggplot2)
library(MOFAdata)
library(dplyr)
library(MOFA2)
library(tidyr)
library(svglite)
library(extrafont)

set.seed(1)

# Formatting and path settings
loadfonts(device = "win")
theme_set(theme_bw(base_size = 10, base_family = "Arial"))

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
model_dir <- file.path(script_dir, "trained_model")
meta_dir <- file.path(script_dir, "processed_data")
align_dir <- file.path(script_dir, "aligned_data")

# User setting
model_name <- "Unpaired"  # Choose from "Paired", "Unpaired", "SCOT+"

allowed_models <- c("Paired", "Unpaired", "SCOT+")
if (!(model_name %in% allowed_models)) {
  stop(sprintf("Warning: invalid model '%s'. Choose from: Paired, Unpaired, or SCOT+.", model_name))
}
cat("Model input correct!\n")

outdir <- file.path(script_dir, "GSEA_output", model_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Global settings
top_genes_n  <- 3
max_pathways <- 10
max_genes    <- 8

views <- c("Transcriptomics", "Proteomics")

# Helper functions
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
    
    df_tmp <- data.frame(Factor = paste0("Factor", fac), FactorNumber = fac, Pathway = top_paths,
                         Score = score_signed, Genes = as.numeric(gene_counts), Sign = sign_label,
                         stringsAsFactors = FALSE)
    
    df_list[[length(df_list) + 1]] <- df_tmp
  }
  
  bind_rows(df_list)
}

extract_factor_weights <- function(mofa_model, view, factor) {
  w <- get_weights(mofa_model, views = view, factors = factor, as.data.frame = TRUE)
  
  if (!all(c("feature", "value") %in% colnames(w))) {
    stop("Expected columns 'feature' and 'value' not found in weights table.")
  }
  
  w %>%
    mutate(feature_clean = clean_feature_names(feature), abs_value = abs(value),
           sign = case_when(value > 0 ~ "positive", value < 0 ~ "negative", TRUE ~ "zero")
    ) %>%
    arrange(desc(abs_value))
}

prepare_model_for_gsea <- function(mofa_model, views, feature_sets) {
  model_gsea <- mofa_model
  
  for (v in views) {
    feat <- features_names(model_gsea)[[v]]
    feat_clean <- toupper(sub("_(rna|prot)$", "", feat))
    
    if (any(duplicated(feat_clean))) {
      warning("Duplicated feature names detected after cleaning in view: ", v,
              ". This may affect feature-set matching.")
    }
    
    features_names(model_gsea)[[v]] <- feat_clean
  }
  
  for (v in views) {
    model_feats <- features_names(model_gsea)[[v]]
    gs_feats <- colnames(feature_sets)
    
    cat("\nView:", v, "\n")
    cat("Model features:", length(model_feats), "\n")
    cat("Gene-set features:", length(gs_feats), "\n")
    cat("Overlap:", sum(model_feats %in% gs_feats), "\n")
  }
  
  return(model_gsea)
}

run_gsea_for_model <- function(mofa_model, outdir, views, feature_sets,
                               max_pathways = 10, max_genes = 8) {
  factors_all <- seq_len(get_dimensions(mofa_model)$K)
  model_gsea <- prepare_model_for_gsea(mofa_model, views, feature_sets)
  
  for (view in views) {
    message("Processing view: ", view)
    
    view_outdir <- file.path(outdir, view)
    dir.create(view_outdir, recursive = TRUE, showWarnings = FALSE)
    
    enrich_pos <- run_enrichment(model_gsea, view = view, factors = factors_all,
                                 feature.sets = feature_sets, sign = "positive",
                                 statistical.test = "parametric")
    
    enrich_neg <- run_enrichment(model_gsea, view = view, factors = factors_all,
                                 feature.sets = feature_sets, sign = "negative",
                                 statistical.test = "parametric")
    
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
      labs(x = "MOFA+ Factor", y = "Enriched Pathway",color = "Signed enrichment score",
           size = "# Genes")
    
    dot_file <- file.path(view_outdir, paste0(view, "_combined_dotplot.svg"))
    ggsave(dot_file, p_dot, device = svglite, width = 10, height = 9)
    fwrite(df_dot, file.path(view_outdir, paste0(view, "_combined_dotdf.csv")))
    message("  Saved combined dot plot: ", dot_file)
    
    weights_by_factor <- lapply(factors_all, function(fac) extract_factor_weights(mofa_model, view, fac))
    names(weights_by_factor) <- paste0("Factor", factors_all)
    
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
        
        p_top <- if (is.null(sig_pathways) || length(sig_pathways) == 0) {
          empty_plot("No significant pathways")
        } else {
          plot_enrichment(enrichment, factor = fac, max.pathways = max_pathways, text_size = 0.7)
        }
        
        ggsave(file.path(fac_dir, "top_enrichment.svg"), p_top,
               device = svglite, width = 6, height = 5)
        
        p_det <- if (is.null(sig_pathways) || length(sig_pathways) == 0) {
          empty_plot("No significant pathways")
        } else {
          plot_enrichment_detailed(enrichment, factor = fac, max.genes = max_genes, max.pathways = 5)
        }
        
        ggsave(file.path(fac_dir, "detailed_enrichment.svg"), p_det,
               device = svglite, width = 7, height = 6)
        
        w <- weights_by_factor[[paste0("Factor", fac)]]
        w_sign <- if (sign == "positive") {
          w %>% filter(value > 0)
        } else {
          w %>% filter(value < 0)
        }
        
        write.csv(w_sign, file = file.path(fac_dir,
                                           paste0(view, "_Factor", fac, "_", sign, "_weights.csv")),
                  row.names = FALSE, quote = FALSE)
        
        if (!is.null(sig_pathways) && length(sig_pathways) > 0) {
          scores_all <- enrichment$set.statistics[, fac]
          sig_paths_valid <- intersect(sig_pathways, names(scores_all))
          
          if (length(sig_paths_valid) > 0) {
            enrich_tbl <- data.frame(Pathway = sig_paths_valid,
                                     Score = as.numeric(scores_all[sig_paths_valid]),
                                     Sign = sign, Factor = fac, View = view,
                                     stringsAsFactors = FALSE) %>%
              arrange(desc(abs(Score)))
            
            write.csv(enrich_tbl, file = file.path(fac_dir,
                                                   paste0(view, "_Factor", fac, "_", sign, "_enrichment_table.csv")),
                      row.names = FALSE, quote = FALSE)
          }
        }
      }
    }
  }
}

# Load gene set data
data("MSigDB_v6.0_C5_mouse")

# Main
if (model_name == "Paired") {
  mofa_model <- load_model(file.path(model_dir, "Paired_model.hdf5"))
  
  rna_metadata <- read.csv( file.path(meta_dir, "Processed_transcriptomics_metadata.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  metadata <- rna_metadata %>%
    distinct(SampleID, .keep_all = TRUE)
  
  mofa_model <- attach_metadata_to_model(mofa_model, metadata = metadata, sample_col = "SampleID")
  
} else if (model_name == "SCOT+") {
  mofa_model <- load_model(file.path(model_dir, "SCOT+-aligned_model.hdf5"))
  
  aligned_metadata <- read.csv(file.path(align_dir, "aligned_metadata.csv"),
                               stringsAsFactors = FALSE, check.names = FALSE)
  
  mofa_model <- attach_metadata_to_model(mofa_model, metadata = aligned_metadata,
                                         sample_col = "SampleID")
  
} else if (model_name == "Unpaired") {
  mofa_model <- load_model(file.path(model_dir, "Unpaired_model.hdf5"))
  
  rna_metadata <- read.csv(file.path(meta_dir, "Processed_transcriptomics_metadata.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)
  prot_metadata <- read.csv(file.path(meta_dir, "Processed_proteomics_metadata.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
  
  if ("CellLines" %in% colnames(rna_metadata) && !"CellLine" %in% colnames(rna_metadata)) {
    rna_metadata <- rna_metadata %>% rename(CellLine = CellLines)
  }
  if ("CellLines" %in% colnames(prot_metadata) && !"CellLine" %in% colnames(prot_metadata)) {
    prot_metadata <- prot_metadata %>% rename(CellLine = CellLines)
  }
  
  rna_metadata$SampleID  <- paste0(as.character(rna_metadata$SampleID), "_rna")
  prot_metadata$SampleID <- paste0(as.character(prot_metadata$SampleID), "_prot")
  rna_metadata$Modality  <- "Transcriptomics"
  prot_metadata$Modality <- "Proteomics"
  
  metadata <- bind_rows(prot_metadata, rna_metadata) %>%
    distinct(SampleID, .keep_all = TRUE)
  
  mofa_model <- attach_metadata_to_model(mofa_model, metadata = metadata, sample_col = "SampleID")
}

run_gsea_for_model(mofa_model = mofa_model, outdir = outdir, views = views,
                   feature_sets = MSigDB_v6.0_C5_mouse, max_pathways = max_pathways,
                   max_genes = max_genes)