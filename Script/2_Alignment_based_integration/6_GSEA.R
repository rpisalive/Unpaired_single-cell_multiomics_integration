library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(MOFAdata)
library(dplyr)
library(MOFA2)

set.seed(1)
outdir <- "PATH_TO_OUTPUT_DIRECTORY"

#Global plotting defaults
loadfonts(device = "win")
theme_set(theme_bw(base_size = 12, base_family = "Arial"))

#Parameters
top_genes_n  <- 3
max_pathways <- 10
max_genes    <- 8
setwd("PATH_TO_WORKING_DIRECTORY")
model <- load_model("SCOTplus_trained_model.hdf5")
model@samples_metadata <- model@samples_metadata %>%
  mutate(cell_line = sub("^[^_]+_([^_]+).*", "\\1", sample))
data("MSigDB_v6.0_C5_mouse")

#Settings
views  <- c("Transcriptomics","Proteomics")
factors_all <- seq_len(get_dimensions(model)$K)

#Normalize feature names once
for(v in views){
  suffix <- if(v=="Transcriptomics") "_RNA$" else "_PROT$"
  features_names(model)[[v]] <- toupper(features_names(model)[[v]])
  features_names(model)[[v]] <- sub(suffix, "", features_names(model)[[v]])
}


#Helper functions
empty_plot <- function(label){
  ggplot() +
    annotate("text", 0.5, 0.5, label=label, size=5, hjust=0.5, vjust=0.5) +
    theme_void()
}

#Build dotplot dataframe
build_dotplot_df <- function(enrichment, sign_label, top_n = 10){
  
  df_list <- list()
  
  for(fac in seq_along(enrichment$sigPathways)){
    
    sig_paths <- enrichment$sigPathways[[fac]]
    if(length(sig_paths) == 0) next
    
    # Extract all scores for this factor
    scores_all <- enrichment$set.statistics[, fac]
    
    # Keep only pathways that exist in the statistics
    sig_paths_valid <- intersect(sig_paths, names(scores_all))
    if(length(sig_paths_valid) == 0) next
    
    scores <- scores_all[sig_paths_valid]
    
    # Select top pathways by absolute score
    top_paths <- names(sort(abs(scores), decreasing = TRUE))[1:min(top_n, length(scores))]
    if(length(top_paths) == 0) next
    
    # Gene counts (ensure matching rows exist)
    valid_gene_paths <- intersect(top_paths, rownames(enrichment$feature.sets))
    if(length(valid_gene_paths) == 0) next
    
    gene_counts <- rowSums(enrichment$feature.sets[valid_gene_paths, , drop = FALSE])
    
    # Align everything to the same pathways
    top_paths <- valid_gene_paths
    scores_subset <- scores[top_paths]
    
    # Apply sign correction
    if(sign_label == "negative") scores_subset <- -scores_subset
    
    # Final safety check
    if(length(top_paths) == 0 || length(scores_subset) == 0 || length(gene_counts) == 0) next
    
    df_tmp <- data.frame(
      Factor  = paste0("Factor", fac),
      Pathway = top_paths,
      Score   = as.numeric(scores_subset),
      Genes   = as.numeric(gene_counts),
      stringsAsFactors = FALSE
    )
    
    df_list[[length(df_list) + 1]] <- df_tmp
  }
  
  df <- bind_rows(df_list)
  return(df)
}

#Loop per view
for(view in views){
  message("Processing view: ", view)
  view_outdir <- file.path(outdir, view)
  dir.create(view_outdir, recursive = TRUE, showWarnings = FALSE)
  #Run enrichment
  enrich_pos <- run_enrichment(model, view = view, factors = factors_all, feature.sets = MSigDB_v6.0_C5_mouse,
                               sign = "positive", statistical.test = "parametric")
  
  enrich_neg <- run_enrichment(model, view = view, factors = factors_all, feature.sets = MSigDB_v6.0_C5_mouse,
                               sign = "negative", statistical.test = "parametric")
  
  #Build combined dotplot dataframe
  df_pos <- build_dotplot_df(enrich_pos, "positive", max_pathways)
  df_neg <- build_dotplot_df(enrich_neg, "negative", max_pathways)
  df_dot <- bind_rows(df_pos, df_neg)
  
  if(nrow(df_dot) == 0){
    message("  No significant pathways")
    next}
  df_dot <- df_dot[order(df_dot$Score), ]
  df_dot$Pathway <- factor(df_dot$Pathway, levels = unique(df_dot$Pathway))
  p_dot <- ggplot(df_dot, aes(x = Factor, y = Pathway, size = Genes, color = Score)) +
    geom_point() +
    scale_x_discrete(drop = FALSE) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          guide = guide_colorbar(order = 1)) +
    scale_size(range = c(2,6), guide = guide_legend(order = 2)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 8, face = "bold"),
          axis.title.y = element_text(size = 8, face = "bold"),
          legend.title = element_text(size = 8, face = "bold"),
          legend.text  = element_text(size = 8)) +
    labs(x = "MOFA+ Factor", y = "Enriched Pathway",
         color = "Signed enrichment score", size = "# Genes")
  
  dot_file <- file.path(view_outdir, paste0(view, "_combined_dotplot.svg"))
  ggsave(dot_file, p_dot, device=svglite, width = 10, height = 9)
  fwrite(df_dot, file.path(view_outdir, paste0(view, "_combined_dotdf.csv")))
  message("  Saved combined dot plot: ", dot_file)
  
  #Per factor detailed plots
  for(sign in c("positive","negative")){
    enrichment <- if(sign=="positive") enrich_pos else enrich_neg
    # select valid factors
    valid_factors <- factors_all[
      sapply(factors_all, function(f){
        w <- get_weights(model, views=view, factors=f, as.data.frame=TRUE)
        if(nrow(w)==0) return(FALSE)
        if(sign=="positive") sum(w$value>0) > 5 else sum(w$value<0) > 5
      })
    ]
    
    for(fac in valid_factors){
      fac_dir <- file.path(view_outdir, paste0(sign,"_Factor", fac))
      dir.create(fac_dir, recursive = TRUE, showWarnings = FALSE)
      sig_pathways <- enrichment$sigPathways[[fac]]
      #Top pathways
      p_top <- if(length(sig_pathways)==0){
        empty_plot("No significant pathways")
      } else {
        plot_enrichment(enrichment, factor=fac, max.pathways=max_pathways, text_size=0.7)
      }
      
      ggsave(file.path(fac_dir,"top_enrichment.svg"), p_top, device=svglite, width=6, height=5)
      
      #Detailed gene-level
      p_det <- if(length(sig_pathways)==0){
        empty_plot("No significant pathways")
      } else {
        plot_enrichment_detailed(enrichment, factor=fac, max.genes=max_genes, max.pathways=5)
      }
      ggsave(file.path(fac_dir,"detailed_enrichment.svg"), p_det, device=svglite, width=7, height=6)
    }
  }
}