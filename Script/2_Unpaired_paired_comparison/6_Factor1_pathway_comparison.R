library(MOFA2)
library(MOFAdata)
library(dplyr)
library(readr)

outdir <- "C:/Users/49152/Downloads/Multi-omics/MOFA/output/graphs/Unpaired_single_cell/Minimal_baseline/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

m_ref <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/single_cell_trained_model_MOFAcomplied.hdf5")
m_unp <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/Unpaired_single_cell_model.hdf5")

data("MSigDB_v6.0_C5_mouse")

matched_pairs <- read.csv(file.path(outdir, "matched_factor_pairs_reference_vs_unpaired.csv"),
                          stringsAsFactors = FALSE)

# Use the best matched factor only
top_pair <- matched_pairs %>%
  arrange(desc(MeanAbsCorrelation)) %>%
  slice(1)

ref_factor <- as.integer(sub("Reference_F", "", top_pair$ReferenceFactor))
unp_factor <- as.integer(sub("Unpaired_F", "", top_pair$UnpairedFactor))

normalize_features <- function(model) {
  for (v in names(features_names(model))) {
    feats <- features_names(model)[[v]]
    feats <- toupper(feats)
    feats <- sub("_(RNA|PROT)$", "", feats)
    features_names(model)[[v]] <- feats
  }
  model
}

m_ref <- normalize_features(m_ref)
m_unp <- normalize_features(m_unp)

views <- c("Transcriptomics", "Proteomics")
signs <- c("positive", "negative")

get_sig_paths <- function(model, view, factor, sign) {
  enrich <- tryCatch(
    run_enrichment(
      model,
      view = view,
      factors = factor,
      feature.sets = MSigDB_v6.0_C5_mouse,
      sign = sign,
      statistical.test = "parametric"
    ),
    error = function(e) NULL
  )
  if (is.null(enrich) || is.null(enrich$sigPathways)) return(character(0))
  unique(unlist(enrich$sigPathways))
}

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

res <- list()

for (v in views) {
  for (s in signs) {
    ref_paths <- get_sig_paths(m_ref, v, ref_factor, s)
    unp_paths <- get_sig_paths(m_unp, v, unp_factor, s)
    
    res[[paste(v, s, sep = "_")]] <- data.frame(
      View = v,
      Sign = s,
      ReferenceFactor = top_pair$ReferenceFactor,
      UnpairedFactor = top_pair$UnpairedFactor,
      MeanAbsCorrelation = top_pair$MeanAbsCorrelation,
      n_ref = length(ref_paths),
      n_unpaired = length(unp_paths),
      n_intersection = length(intersect(ref_paths, unp_paths)),
      n_union = length(union(ref_paths, unp_paths)),
      Jaccard = jaccard(ref_paths, unp_paths),
      stringsAsFactors = FALSE
    )
  }
}

pathway_summary <- bind_rows(res)

write.csv(pathway_summary, file.path(outdir, "TopFactor_pathway_comparison_reference_vs_unpaired.csv"),
          row.names = FALSE)

print(pathway_summary)