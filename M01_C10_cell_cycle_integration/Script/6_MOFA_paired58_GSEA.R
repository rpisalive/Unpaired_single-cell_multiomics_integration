# =========================================================
# GSEA of treatment-associated factors from native paired58 MOFA+ NanoSPLITS C10 RO-3306 experiment: RNA, protein and phosphopeptide
# =========================================================
#
# Run 5_MOFA_paired58_exploratory_downstream.R first.
#
# This script uses the condition-oriented MOFA+ weights exported by that script. Therefore:
#   positive NES = enrichment toward RO-3306-treated/G2M-arrested cells;
#   negative NES = enrichment toward untreated Control cells.
#
# RNA features are already genes. Protein features are mapped to genes and duplicate protein rows for a gene are averaged.
# Phosphopeptide features are mapped to genes and the peptide with the largest absolute loading represents each gene.
# The phosphopeptide GSEA is exploratory because site-specific and even opposing phosphoregulation can be lost during gene-level collapse.
# Preserve phosphopeptide loading tables for mechanistic/site-level validation.
#
# Required packages:
#   tidyverse, fgsea, msigdbr, svglite
#
# =========================================================
# 1. Packages and paths
# =========================================================

required_packages <- c("tidyverse", "fgsea", "msigdbr", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Install the following required package(s) before running this script: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
  library(svglite)
})

get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    context <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)

    if (!is.null(context) && !is.null(context$path) && nzchar(context$path)) {
      return(normalizePath(dirname(context$path), winslash = "/", mustWork = TRUE))
    }
  }

  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)

  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    return(normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE))
  }

  warning("Could not identify the script directory; using the working directory.")
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
processed_dir <- file.path(script_dir, "processed_data")
exploratory_dir <- file.path(processed_dir, "mofa_downstream", "native_paired58", "exploratory")
gsea_dir <- file.path(processed_dir, "mofa_downstream", "native_paired58", "gsea")

dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

weights_file <- file.path(exploratory_dir, "mofa_weights_oriented_annotated.csv")
association_file <- file.path(exploratory_dir, "mofa_factor_condition_associations.csv")

required_files <- c(weights_file, association_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Required input file(s) do not exist. Run the exploratory script first:\n", paste(missing_files, collapse = "\n"))
}

# =========================================================
# 2. Analysis settings
# =========================================================

set.seed(1)

expected_views <- c("RNA", "Protein", "Phosphopeptide")

# Restrict primary GSEA to factors associated with condition. This reduces
# exploratory multiplicity and keeps the manuscript analysis focused.
run_all_factors <- FALSE
condition_factor_fdr <- 0.05

# Hallmark is the compact primary library. Reactome adds pathway detail.
gene_set_libraries <- tribble(~Library, ~Collection, ~Subcollection, "Hallmark", "H", NA_character_, "Reactome", "C2", "CP:REACTOME")

min_pathway_size <- 10L
max_pathway_size <- 500L
gsea_fdr_threshold <- 0.05
top_pathways_per_direction <- 5L

base_family <- "Arial"
theme_set(theme_bw(base_size = 10, base_family = base_family))

# =========================================================
# 3. Helper functions
# =========================================================

require_columns <- function(data, required, object_name) {
  missing <- setdiff(required, colnames(data))

  if (length(missing) > 0) {
    stop(object_name, " is missing required column(s): ", paste(missing, collapse = ", "))
  }

  invisible(TRUE)
}

load_msigdb_collection <- function(collection, subcollection = NA_character_) {
  available_arguments <- names(formals(msigdbr::msigdbr))
  arguments <- list(species = "Mus musculus")

  # H and C2 are human MSigDB collection codes. Request the human database explicitly and let msigdbr map its gene sets to mouse orthologues.
  if ("db_species" %in% available_arguments) {
    arguments$db_species <- "HS"
  }

  if ("collection" %in% available_arguments) {
    arguments$collection <- collection
  } else {
    arguments$category <- collection
  }

  if (!is.na(subcollection) && nzchar(subcollection)) {
    if ("subcollection" %in% available_arguments) {
      arguments$subcollection <- subcollection
    } else {
      arguments$subcategory <- subcollection
    }
  }

  output <- do.call(msigdbr::msigdbr, arguments)

  require_columns(output, c("gs_name", "gene_symbol"), "msigdbr output")

  output %>%
    transmute(
      Pathway = as.character(gs_name),
      Gene = trimws(as.character(gene_symbol))
    ) %>%
    filter(!is.na(Gene), nzchar(Gene)) %>%
    distinct(Pathway, Gene)
}

collapse_weights_to_genes <- function(weights) {
  expanded <- weights %>% mutate(Gene = as.character(Gene)) %>% tidyr::separate_rows(Gene, sep = "[;,]") %>% mutate(Gene = trimws(Gene)) %>%
    filter(!is.na(Gene), nzchar(Gene), is.finite(OrientedWeight))

  rna_protein <- expanded %>% filter(View %in% c("RNA", "Protein")) %>% group_by(View, Factor, Gene) %>% 
    summarise(RankScore = mean(OrientedWeight), SourceFeatureCount = n_distinct(SourceFeatureID),
              CollapseMethod = if_else(first(View) == "RNA", "identity_or_mean_if_duplicated", "mean_across_protein_rows"), .groups = "drop")

  phosphopeptide <- expanded %>% filter(View == "Phosphopeptide") %>% group_by(View, Factor, Gene) %>%
    mutate(SourceFeatureCountBeforeCollapse = n_distinct(SourceFeatureID)) %>%
    slice_max(order_by = abs(OrientedWeight), n = 1, with_ties = FALSE) %>%
    summarise(RankScore = first(OrientedWeight), SourceFeatureCount = first(SourceFeatureCountBeforeCollapse),
              CollapseMethod = "maximum_absolute_phosphopeptide_weight", RepresentativeFeatureID = first(SourceFeatureID),
              RepresentativeModifiedSequence = first(ModifiedSequence), .groups = "drop")

  bind_rows(rna_protein, phosphopeptide) %>% filter(is.finite(RankScore)) %>% arrange(View, Factor, desc(RankScore))
}

make_rank_vector <- function(rank_table) {
  if (anyDuplicated(rank_table$Gene)) {
    stop("Gene-level rank table contains duplicated gene symbols.")
  }

  ranks <- rank_table$RankScore
  names(ranks) <- rank_table$Gene
  ranks <- ranks[is.finite(ranks) & !is.na(names(ranks)) & nzchar(names(ranks))]
  sort(ranks, decreasing = TRUE)
}

run_one_gsea <- function(rank_table, pathways, view, factor, library_name) {
  ranks <- make_rank_vector(rank_table)
  genes_in_sets <- unique(unlist(pathways, use.names = FALSE))
  overlap_n <- sum(names(ranks) %in% genes_in_sets)

  overlap_qc <- tibble(View = view, Factor = factor, Library = library_name, RankedGenes = length(ranks),
                       GenesRepresentedInLibrary = overlap_n, RankOverlapFraction = overlap_n / length(ranks))

  if (overlap_n < min_pathway_size) {
    warning("Skipping ", view, " / ", factor, " / ", library_name, ": only ", overlap_n, " ranked genes overlap the gene-set library.")

    return(list(results = tibble(), overlap = overlap_qc))
  }

  result <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks, minSize = min_pathway_size, maxSize = max_pathway_size, eps = 0, nproc = 1) %>%
    as_tibble() %>%  mutate(leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"), View = view, Factor = factor,
                            Library = library_name, Direction = case_when(NES > 0 ~ "Treated-associated", NES < 0 ~ "Control-associated",
                                                                          TRUE ~ "Neutral")) %>%
    select(View, Factor, Library, pathway, NES, ES, pval, padj, size, Direction, leadingEdge, everything()) %>% arrange(padj, desc(abs(NES)))

  list(results = result, overlap = overlap_qc)
}

make_readable_pathway <- function(pathway) {
  pathway %>% str_remove("^HALLMARK_") %>% str_remove("^REACTOME_") %>% str_replace_all("_", " ") %>% str_to_sentence()
}

save_svg <- function(filename, plot, width, height) {
  ggsave(filename = file.path(gsea_dir, filename), plot = plot, device = svglite::svglite, width = width, height = height)
}

# =========================================================
# 4. Load oriented weights and select factors
# =========================================================

weights <- readr::read_csv(weights_file, show_col_types = FALSE)
associations <- readr::read_csv(association_file, show_col_types = FALSE)

require_columns(
  weights,
  c(
    "View", "Factor", "SourceFeatureID", "Gene", "ModifiedSequence",
    "OrientedWeight"
  ),
  "Oriented MOFA+ weights"
)
require_columns(
  associations,
  c("Factor", "PValue", "FDR", "PrioritizedForInterpretation"),
  "Factor-condition association table"
)

if (!setequal(unique(weights$View), expected_views)) {
  stop(
    "Unexpected view names in the weight table: ",
    paste(sort(unique(weights$View)), collapse = ", ")
  )
}

if (run_all_factors) {
  factors_to_test <- associations %>%
    arrange(PValue) %>%
    pull(Factor)
} else {
  factors_to_test <- associations %>%
    filter(FDR <= condition_factor_fdr) %>%
    arrange(PValue) %>%
    pull(Factor)
}

if (length(factors_to_test) == 0) {
  factors_to_test <- associations %>%
    slice_min(PValue, n = 1, with_ties = FALSE) %>%
    pull(Factor)

  warning(
    "No factor passed condition FDR <= ", condition_factor_fdr,
    ". Running exploratory GSEA only for the lowest-P factor: ",
    factors_to_test
  )
}

cat("Factors selected for GSEA: ", paste(factors_to_test, collapse = ", "), "\n", sep = "")

# =========================================================
# 5. Build one ranked gene list per view and factor
# =========================================================

gene_ranks <- weights %>%
  filter(Factor %in% factors_to_test) %>%
  collapse_weights_to_genes()

if (nrow(gene_ranks) == 0) {
  stop("No gene-level ranks were produced.")
}

rank_counts <- gene_ranks %>%
  count(View, Factor, name = "RankedGenes")

missing_view_factor <- tidyr::expand_grid(
  View = expected_views,
  Factor = factors_to_test
) %>%
  anti_join(rank_counts, by = c("View", "Factor"))

if (nrow(missing_view_factor) > 0) {
  stop("At least one selected view/factor combination has no gene-level ranks.")
}

write_csv(
  gene_ranks,
  file.path(gsea_dir, "mofa_gene_level_ranked_weights.csv"),
  na = ""
)

# =========================================================
# 6. Load current mouse MSigDB libraries
# =========================================================

gene_sets <- list()
gene_set_membership <- list()

for (i in seq_len(nrow(gene_set_libraries))) {
  library_name <- gene_set_libraries$Library[[i]]
  collection <- gene_set_libraries$Collection[[i]]
  subcollection <- gene_set_libraries$Subcollection[[i]]

  cat("Loading MSigDB library: ", library_name, "\n", sep = "")

  membership <- load_msigdb_collection(collection, subcollection)

  if (nrow(membership) == 0) {
    stop("No gene sets were returned for library: ", library_name)
  }

  gene_set_membership[[library_name]] <- membership %>%
    mutate(Library = library_name, .before = 1)

  gene_sets[[library_name]] <- split(membership$Gene, membership$Pathway)
}

write_csv(
  bind_rows(gene_set_membership),
  file.path(gsea_dir, "msigdb_mouse_gene_set_membership_used.csv")
)

# =========================================================
# 7. Run preranked GSEA
# =========================================================

gsea_results <- list()
overlap_results <- list()
result_index <- 1L

for (view in expected_views) {
  for (factor in factors_to_test) {
    rank_table <- gene_ranks %>%
      filter(View == view, Factor == factor)

    for (library_name in names(gene_sets)) {
      cat("GSEA: ", view, " | ", factor, " | ", library_name, "\n", sep = "")

      one_result <- run_one_gsea(
        rank_table = rank_table,
        pathways = gene_sets[[library_name]],
        view = view,
        factor = factor,
        library_name = library_name
      )

      gsea_results[[result_index]] <- one_result$results
      overlap_results[[result_index]] <- one_result$overlap
      result_index <- result_index + 1L
    }
  }
}

gsea_table <- bind_rows(gsea_results)
overlap_table <- bind_rows(overlap_results)

if (nrow(gsea_table) == 0) {
  stop("GSEA completed but produced no pathway results. Check gene-symbol overlap.")
}

write_csv(
  gsea_table,
  file.path(gsea_dir, "mofa_gsea_all_results.csv"),
  na = ""
)
write_csv(
  gsea_table %>% filter(!is.na(padj), padj <= gsea_fdr_threshold),
  file.path(gsea_dir, "mofa_gsea_significant_fdr05.csv"),
  na = ""
)
write_csv(
  overlap_table,
  file.path(gsea_dir, "mofa_gsea_gene_overlap_qc.csv")
)

# =========================================================
# 8. Compact enrichment plots and cross-view table
# =========================================================

significant_results <- gsea_table %>%
  filter(!is.na(padj), padj <= gsea_fdr_threshold)

top_positive <- significant_results %>%
  filter(NES > 0) %>%
  group_by(View, Factor, Library) %>%
  slice_min(padj, n = top_pathways_per_direction, with_ties = FALSE) %>%
  ungroup()

top_negative <- significant_results %>%
  filter(NES < 0) %>%
  group_by(View, Factor, Library) %>%
  slice_min(padj, n = top_pathways_per_direction, with_ties = FALSE) %>%
  ungroup()

plot_table <- bind_rows(top_positive, top_negative) %>%
  distinct(View, Factor, Library, pathway, .keep_all = TRUE) %>%
  mutate(
    PathwayLabel = make_readable_pathway(pathway),
    PlotID = paste(Library, pathway, sep = "___"),
    View = factor(View, levels = expected_views),
    Factor = factor(Factor, levels = factors_to_test)
  )

if (nrow(plot_table) > 0) {
  pathway_labels <- setNames(plot_table$PathwayLabel, plot_table$PlotID)

  p_gsea <- plot_table %>%
    ggplot(
      aes(
        x = NES,
        y = reorder(PlotID, NES),
        colour = NES,
        size = -log10(pmax(padj, .Machine$double.xmin))
      )
    ) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
    geom_point(alpha = 0.9) +
    facet_grid(View ~ Factor, scales = "free_y", space = "free_y") +
    scale_colour_gradient2(
      low = "#0072B2",
      mid = "white",
      high = "#D55E00",
      midpoint = 0
    ) +
    scale_y_discrete(labels = function(x) unname(pathway_labels[x])) +
    labs(
      x = "Normalized enrichment score",
      y = NULL,
      colour = "NES",
      size = expression(-log[10](FDR))
    ) +
    theme(
      legend.position = "top",
      axis.text.y = element_text(size = 7),
      panel.grid.major.y = element_blank()
    )

  plot_width <- max(9, 5 * length(factors_to_test))
  plot_height <- max(8, 2.8 * length(expected_views))

  save_svg("mofa_gsea_top_pathways.svg", p_gsea, plot_width, plot_height)
} else {
  warning("No pathways passed GSEA FDR <= ", gsea_fdr_threshold, "; no dot plot was created.")
}

cross_view_nes <- gsea_table %>%
  select(Library, Factor, pathway, View, NES, padj) %>%
  pivot_wider(
    names_from = View,
    values_from = c(NES, padj),
    names_glue = "{View}_{.value}"
  ) %>%
  arrange(Library, Factor, pathway)

write_csv(
  cross_view_nes,
  file.path(gsea_dir, "mofa_gsea_pathway_NES_across_views.csv"),
  na = ""
)

# =========================================================
# 9. Analysis notes and reproducibility records
# =========================================================

factor_selection_note <- if (any(associations$FDR <= condition_factor_fdr, na.rm = TRUE)) {
  paste0(
    "GSEA factors were selected at condition-association FDR <= ",
    condition_factor_fdr, "."
  )
} else {
  paste0(
    "No factor passed condition-association FDR <= ", condition_factor_fdr,
    "; only the lowest-P factor was analysed as explicitly exploratory."
  )
}

analysis_notes <- c(
  "MOFA+ factor signs were oriented in the upstream exploratory script.",
  "Positive NES denotes the RO-3306-treated/G2M-arrested direction; negative NES denotes Control.",
  factor_selection_note,
  "RNA ranks use gene-level weights.",
  "Protein ranks average weights for duplicate protein rows mapping to the same gene.",
  paste0(
    "Phosphopeptide ranks use the peptide with the maximum absolute weight per gene. ",
    "Treat phosphopeptide GSEA as exploratory and retain site-level loadings for mechanistic validation."
  ),
  "No parent-protein adjustment is performed in this GSEA script."
)

write_lines(analysis_notes, file.path(gsea_dir, "GSEA_analysis_notes.txt"))
write_lines(
  capture.output(sessionInfo()),
  file.path(gsea_dir, "GSEA_sessionInfo.txt")
)

summary_report <- tibble(
  Item = c(
    "Factors tested",
    "Views tested",
    "Gene-set libraries",
    "Minimum pathway size",
    "Maximum pathway size",
    "Significance threshold",
    "Significant results"
  ),
  Value = c(
    paste(factors_to_test, collapse = ";"),
    paste(expected_views, collapse = ";"),
    paste(names(gene_sets), collapse = ";"),
    as.character(min_pathway_size),
    as.character(max_pathway_size),
    as.character(gsea_fdr_threshold),
    as.character(sum(!is.na(gsea_table$padj) & gsea_table$padj <= gsea_fdr_threshold))
  )
)

write_csv(summary_report, file.path(gsea_dir, "mofa_gsea_summary.csv"))

cat("\nPaired58 MOFA+ GSEA complete.\n")
cat("Output directory:\n", gsea_dir, "\n", sep = "")
