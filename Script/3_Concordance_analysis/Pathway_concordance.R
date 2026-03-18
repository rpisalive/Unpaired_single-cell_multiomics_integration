library(MOFA2)
library(dplyr)
library(MOFAdata)
library(ggplot2)
library(tidyr)
library(purrr)
library(pheatmap)
library(tibble)

outdir <- "PATH_TO_OUTPUT_DIRECTORY"
m_single <- load_model("PATH/TO/single_cell_trained_model.hdf5")
m_scot   <- load_model("PATH/TO/SCOTplus_trained_model.hdf5")
data("MSigDB_v6.0_C5_mouse")
gene_sets <- MSigDB_v6.0_C5_mouse

#Normalise feature names
normalize_features <- function(model){
  for(v in names(features_names(model))){
    feats <- features_names(model)[[v]]
    feats <- toupper(feats)
    feats <- sub("_RNA$", "", feats)
    feats <- sub("_PROT$", "", feats)
    features_names(model)[[v]] <- feats
  }
  return(model)
}

m_single <- normalize_features(m_single)
m_scot   <- normalize_features(m_scot)

#Parameters
views <- c("Transcriptomics", "Proteomics")
signs <- c("positive", "negative")
n_features <- 100
K_single <- get_dimensions(m_single)$K
K_scot   <- get_dimensions(m_scot)$K

#Helper function to get top features per factor
top_features_factor <- function(model, view, factor, n_features = 100, sign = "positive"){
  w <- get_weights(model, views = view, factors = factor, as.data.frame = TRUE)
  if(nrow(w) == 0) return(character(0))
  
  w$value <- as.numeric(w$value)
  
  if(sign == "positive") w <- w[w$value > 0, ]
  if(sign == "negative") w <- w[w$value < 0, ]
  
  if(nrow(w) == 0) return(character(0))
  
  ord <- order(abs(w$value), decreasing = TRUE)
  top_feats <- w[ord[seq_len(min(n_features, nrow(w)))], ]
  
  return(toupper(top_feats$feature))
}

#Function to run enrichment per model
run_model_enrichment <- function(model, model_name){
  enrichment_list <- list()
  K <- get_dimensions(model)$K
  
  for(view in views){
    for(sign in signs){
      for(f in 1:K){
        features <- top_features_factor(model, view, f, n_features = n_features, sign = sign)
        if(length(features) == 0) next
        
        enrich_res <- tryCatch({
          run_enrichment(
            model,
            view = view,
            factors = f,
            feature.sets = gene_sets,
            sign = sign,
            statistical.test = "parametric"
          )
        }, error=function(e) NULL)
        
        if(!is.null(enrich_res)){
          enrichment_list[[paste(model_name, view, sign, f, sep="_")]] <- enrich_res
        }
      }
    }
  }
  return(enrichment_list)
}

enrich_single <- run_model_enrichment(m_single, "Single")
enrich_scot   <- run_model_enrichment(m_scot, "SCOTplus")
names(enrich_single)
names(enrich_scot)

#Get significant pathways from enrichment object
get_sig_pathways <- function(enrich_obj) {
  if(is.null(enrich_obj)) return(character(0))
  sig_list <- enrich_obj$sigPathways
  unique(unlist(sig_list))
}

#Compute Jaccard index between two sets
jaccard_index <- function(set1, set2) {
  if(length(set1) == 0 & length(set2) == 0) return(1)
  if(length(set1) == 0 | length(set2) == 0) return(0)
  length(intersect(set1, set2)) / length(union(set1, set2))
}

#Flatten enrichment list
flatten_enrichment <- function(enrich_list) {
  
  df_list <- lapply(names(enrich_list), function(nm) {
    obj <- enrich_list[[nm]]
    
    #Extract model, view, sign, factor from the list name
    parts <- strsplit(nm, "_")[[1]]
    model_name <- parts[1]
    view <- parts[2]
    sign <- parts[3]
    factor <- parts[4]
    
    paths <- get_sig_pathways(obj)
    if(length(paths) == 0) return(NULL)
    
    data.frame(
      model  = model_name,
      view   = view,
      sign   = sign,
      factor = factor,
      pathway = paths,
      stringsAsFactors = FALSE
    )
  })
  
  df <- do.call(rbind, df_list)
  rownames(df) <- NULL
  return(df)
}


#Flatten both enrichment lists
df_single <- flatten_enrichment(enrich_single)
df_scot   <- flatten_enrichment(enrich_scot)

#Combine and prepare for Jaccard calculation
df_combined <- bind_rows(df_single, df_scot)

#Ensure factors are uniquely identified per model
df_combined <- df_combined %>%
  mutate(factor_id = paste(model, view, sign, factor, sep="_"))

#Create pathway presence matrix
pathway_mat <- df_combined %>%
  select(factor_id, pathway) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = pathway, values_from = value, values_fill = 0) %>%
  column_to_rownames("factor_id")

#Compute Jaccard similarity between all factor pairs
jaccard_similarity <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x | y)
  if(union == 0) return(NA)
  intersection / union
}

factor_ids <- rownames(pathway_mat)
jaccard_mat <- matrix(0, nrow = length(factor_ids), ncol = length(factor_ids),
                      dimnames = list(factor_ids, factor_ids))

for(i in seq_along(factor_ids)) {
  for(j in seq_along(factor_ids)) {
    jaccard_mat[i,j] <- jaccard_similarity(pathway_mat[i,], pathway_mat[j,])
  }
}

cross_mat <- jaccard_mat[grep("Single", rownames(jaccard_mat)), grep("SCOTplus", colnames(jaccard_mat))]

shorten_factor_names <- function(fnames) {
  fnames <- gsub("^Single", "Ref", fnames)
  fnames <- gsub("^SCOTplus", "SCOT", fnames)
  fnames <- gsub("Transcriptomics", "RNA", fnames)
  fnames <- gsub("Proteomics", "PROT", fnames)
  fnames <- gsub("positive", "+", fnames)
  fnames <- gsub("negative", "-", fnames)
  #Convert underscores to structured format, e.g. Ref_RNA_+_1 → Ref-RNA-F1(+)
  fnames <- gsub("_([+-])_([0-9]+)", "-F\\2(\\1)", fnames)
  fnames <- gsub("_", "-", fnames)
  return(fnames)
}

#Apply to cross-model matrix
rownames(cross_mat) <- shorten_factor_names(rownames(cross_mat))
colnames(cross_mat) <- shorten_factor_names(colnames(cross_mat))
write.csv(cross_mat, file.path(outdir, "pathway_concordancesistency.csv"), row.names = TRUE)

#Fig8
svg(file.path(outdir, "pathway_concordance_heatmap.svg"), width = 10, height = 10)
pheatmap(
  cross_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "red"))(50),
  main = paste0(
    "Pathway concordance (Jaccard similarity)\n",
    "Abbreviations: Ref = Reference, SCOT = SCOT+ aligned\n",
    "RNA = Transcriptomic view, PROT = Proteomic view"
  ),
  #main = "Pathway concordance (Jaccard similarity)\nRows: Reference MOFA+   Columns: SCOT+",
  fontsize = 14,
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col = 45,
  border_color = NA,
  legend = TRUE,
  legend_breaks = c(0, 0.5, 1),
  legend_labels = c("0", "0.5", "1")
)

dev.off()