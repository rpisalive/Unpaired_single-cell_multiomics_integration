# =========================================================
# Exploratory downstream analysis of native paired58 MOFA+ NanoSPLITS C10 RO-3306 experiment: RNA, protein and phosphopeptide
# =========================================================
#
# This script is intentionally restricted to the native paired 58-cell model.
# It:
#   1. loads the trained MOFA+ model and paired metadata;
#   2. validates sample and view identities;
#   3. quantifies variance explained by every factor and view;
#   4. tests factor scores for Control versus Treated association;
#   5. orients factor scores and weights so positive values point toward Treated;
#   6. annotates protein and phosphopeptide weights with gene/feature metadata;
#   7. exports compact tables and manuscript-ready exploratory plots.
#
# Factor signs in MOFA+ are arbitrary. Multiplying a factor score and all of its weights by -1 does not change the model.
# The oriented exports below make positive scores/weights consistently represent the RO-3306-treated direction.
# GSEA should use these oriented weights, not the raw MOFA+ signs.
#
# =========================================================
# 1. Packages and paths
# =========================================================

required_packages <- c("MOFA2", "tidyverse", "ggplot2", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Install the following required package(s) before running this script: ", paste(missing_packages, collapse = ", "))
}

library(MOFA2)
library(tidyverse)
library(ggplot2)
library(svglite)

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
model_dir <- file.path(processed_dir, "mofa_models", "native_paired58")
shared_dir <- file.path(processed_dir, "shared")
downstream_dir <- file.path(processed_dir, "mofa_downstream", "native_paired58", "exploratory")
dir.create(downstream_dir, recursive = TRUE, showWarnings = FALSE)

model_file <- file.path(model_dir, "native_paired58_model.hdf5")
metadata_file <- file.path(model_dir, "native_paired58_metadata.csv")
protein_annotation_file <- file.path(shared_dir, "protein_feature_metadata_preprocessed.csv")
phosphopeptide_annotation_file <- file.path(shared_dir, "phosphopeptide_feature_metadata_preprocessed.csv")

required_files <- c(model_file, metadata_file, protein_annotation_file, phosphopeptide_annotation_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Required input file(s) do not exist:\n", paste(missing_files, collapse = "\n"))
}

# =========================================================
# 2. Analysis settings
# =========================================================

expected_views <- c("RNA", "Protein", "Phosphopeptide")
expected_n_samples <- 58L
expected_condition_counts <- c(Control = 27L, Treated = 31L)

condition_levels <- c("Control", "Treated")
condition_colours <- c(Control = "#0072B2", Treated = "#D55E00")
direction_colours <- c("Control-associated" = "#0072B2", "Treated-associated" = "#D55E00")

condition_factor_fdr <- 0.05
top_features_each_direction <- 10L
base_family <- "Arial"

theme_set(theme_bw(base_size = 11, base_family = base_family))

# =========================================================
# 3. Helper functions
# =========================================================

validate_unique_ids <- function(ids, object_name) {
  ids <- as.character(ids)
  
  if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop(object_name, " contains missing or blank identifiers.")
  }
  
  duplicated_ids <- unique(ids[duplicated(ids)])
  
  if (length(duplicated_ids) > 0) {
    stop(object_name, " contains duplicated identifiers: ", paste(head(duplicated_ids, 20), collapse = ", "))
  }
  
  invisible(TRUE)
}

require_columns <- function(data, required, object_name) {
  missing <- setdiff(required, colnames(data))
  
  if (length(missing) > 0) {
    stop(object_name, " is missing required column(s): ", paste(missing, collapse = ", "))
  }
  
  invisible(TRUE)
}

standardize_factor_name <- function(x) {
  x <- as.character(x)
  number <- stringr::str_extract(x, "[0-9]+")
  
  if (anyNA(number)) {
    stop("Could not extract a factor number from: ", paste(unique(x[is.na(number)]), collapse = ", "))
  }
  
  paste0("Factor", as.integer(number))
}

attach_metadata_to_model <- function(mofa_model, metadata, sample_col = "SampleID") {
  model_samples <- unname(unlist(MOFA2::samples_names(mofa_model), use.names = FALSE))
  model_samples <- as.character(model_samples)
  
  validate_unique_ids(model_samples, "MOFA+ model sample names")
  validate_unique_ids(metadata[[sample_col]], "Paired metadata sample IDs")
  
  missing_metadata <- setdiff(model_samples, metadata[[sample_col]])
  extra_metadata <- setdiff(metadata[[sample_col]], model_samples)
  
  if (length(missing_metadata) > 0 || length(extra_metadata) > 0) {
    stop("Model and metadata sample IDs do not match.\n", "Missing from metadata: ",
         ifelse(length(missing_metadata) == 0, "none", paste(missing_metadata, collapse = ", ")), "\nExtra in metadata: ",
         ifelse(length(extra_metadata) == 0, "none", paste(extra_metadata, collapse = ", ")))
  }
  
  metadata_aligned <- metadata[match(model_samples, metadata[[sample_col]]), , drop = FALSE]
  
  if (!identical(as.character(metadata_aligned[[sample_col]]), model_samples)) {
    stop("Internal error while aligning metadata to the model sample order.")
  }
  
  metadata_for_join <- metadata_aligned %>%
    rename(sample = all_of(sample_col)) %>%
    mutate(sample = as.character(sample))
  
  columns_to_replace <- setdiff(colnames(metadata_for_join), "sample")
  
  mofa_model@samples_metadata <- mofa_model@samples_metadata %>%
    mutate(sample = as.character(sample)) %>%
    select(-any_of(columns_to_replace)) %>%
    left_join(metadata_for_join, by = "sample")
  
  list(model = mofa_model, metadata = metadata_aligned)
}

fit_factor_condition_association <- function(data) {
  data <- data %>% mutate(Condition = factor(as.character(Condition), levels = condition_levels))
  
  if (anyNA(data$Condition)) {
    stop("Unexpected condition label in factor-score table.")
  }
  
  fit <- lm(value ~ Condition, data = data)
  coefficient_table <- summary(fit)$coefficients
  
  coefficient_name <- "ConditionTreated"
  
  if (!coefficient_name %in% rownames(coefficient_table)) {
    stop("Could not extract the Treated-versus-Control coefficient.")
  }
  
  beta <- unname(coefficient_table[coefficient_name, "Estimate"])
  standard_error <- unname(coefficient_table[coefficient_name, "Std. Error"])
  p_value <- unname(coefficient_table[coefficient_name, "Pr(>|t|)"])
  critical_value <- qt(0.975, df = df.residual(fit))
  ci_low <- beta - critical_value * standard_error
  ci_high <- beta + critical_value * standard_error
  
  control_scores <- data$value[data$Condition == "Control"]
  treated_scores <- data$value[data$Condition == "Treated"]
  
  pooled_sd <- sqrt(((length(control_scores) - 1) * var(control_scores) +
                       (length(treated_scores) - 1) * var(treated_scores)) /
                      (length(control_scores) + length(treated_scores) - 2))
  
  cohens_d <- if (is.finite(pooled_sd) && pooled_sd > 0) beta / pooled_sd else NA_real_
  orientation_sign <- ifelse(is.finite(beta) && beta < 0, -1, 1)
  
  tibble(ControlN = length(control_scores), TreatedN = length(treated_scores), ControlMean = mean(control_scores),
         TreatedMean = mean(treated_scores), TreatedMinusControl = beta, StandardError = standard_error, CI95Low = ci_low,
         CI95High = ci_high, CohensD = cohens_d, PValue = p_value, OrientationSign = orientation_sign,
         OrientedTreatedMinusControl = beta * orientation_sign,
         OrientedCI95Low = min(ci_low * orientation_sign, ci_high * orientation_sign),
         OrientedCI95High = max(ci_low * orientation_sign, ci_high * orientation_sign),
         OrientedCohensD = cohens_d * orientation_sign)
}

make_feature_annotation <- function(weights, protein_file, phosphopeptide_file) {
  protein_metadata <- readr::read_csv(protein_file, show_col_types = FALSE)
  phosphopeptide_metadata <- readr::read_csv(phosphopeptide_file, show_col_types = FALSE)
  
  require_columns(protein_metadata, c("ProteinRowID", "Protein", "Protein.ID", "Entry.Name", "Gene"),
                  "Protein feature metadata")
  require_columns(phosphopeptide_metadata, c("PeptideRowID", "Modified.Sequence", "Assigned.Modifications", "Protein",
                                             "Entry.Name", "Gene"), "Phosphopeptide feature metadata")
  
  rna_annotation <- weights %>% filter(View == "RNA") %>% distinct(SourceFeatureID) %>%
    transmute(View = "RNA", SourceFeatureID, Gene = SourceFeatureID, EntryName = NA_character_,
              ModifiedSequence = NA_character_, AssignedModifications = NA_character_, FeatureLabel = SourceFeatureID)
  
  protein_annotation <- protein_metadata %>% transmute(View = "Protein", SourceFeatureID = as.character(ProteinRowID),
                                                       Gene = na_if(trimws(as.character(Gene)), ""),
                                                       EntryName = na_if(trimws(as.character(Entry.Name)), ""),
                                                       ModifiedSequence = NA_character_,
                                                       AssignedModifications = NA_character_,
                                                       FeatureLabel = case_when(!is.na(Gene) ~ Gene, !is.na(EntryName) ~ EntryName,
                                                                                TRUE ~ as.character(Protein)))
  
  phosphopeptide_annotation <- phosphopeptide_metadata %>%
    transmute(View = "Phosphopeptide", SourceFeatureID = as.character(PeptideRowID), 
              Gene = na_if(trimws(as.character(Gene)), ""), EntryName = na_if(trimws(as.character(Entry.Name)), ""),
              ModifiedSequence = na_if(trimws(as.character(Modified.Sequence)), ""),
              AssignedModifications = na_if(trimws(as.character(Assigned.Modifications)), ""),
              FeatureLabel = case_when(!is.na(Gene) & !is.na(ModifiedSequence) ~ paste(Gene, ModifiedSequence, sep = " | "),
                                       !is.na(ModifiedSequence) ~ ModifiedSequence, !is.na(Gene) ~ Gene, TRUE ~ SourceFeatureID))
  
  annotation <- bind_rows(rna_annotation, protein_annotation, phosphopeptide_annotation)
  
  duplicated_keys <- annotation %>%
    count(View, SourceFeatureID) %>%
    filter(n > 1)
  
  if (nrow(duplicated_keys) > 0) {
    stop("Feature annotation contains duplicated view/feature-ID keys.")
  }
  
  annotation
}

save_svg <- function(filename, plot, width, height) {
  ggsave(filename = file.path(downstream_dir, filename), plot = plot, device = svglite::svglite, width = width, height = height)
}

# =========================================================
# 4. Load and validate the model and metadata
# =========================================================

cat("Loading paired58 MOFA+ model...\n")
mofa_model <- MOFA2::load_model(model_file)

model_views <- MOFA2::views_names(mofa_model)

if (!setequal(model_views, expected_views)) {
  stop("Unexpected MOFA+ views. Observed: ", paste(model_views, collapse = ", "), "; expected: ",
       paste(expected_views, collapse = ", "))
}

paired_metadata <- readr::read_csv(metadata_file, show_col_types = FALSE)
require_columns(paired_metadata, c("SampleID", "Condition", "ConditionLabel", "Block", "WellIndex"), "Paired58 metadata")
validate_unique_ids(paired_metadata$SampleID, "Paired58 metadata SampleID")

if (nrow(paired_metadata) != expected_n_samples) {
  stop("Expected 58 metadata rows but found ", nrow(paired_metadata), ".")
}

paired_metadata <- paired_metadata %>% mutate(Condition = factor(as.character(Condition), levels = condition_levels))

if (anyNA(paired_metadata$Condition)) {
  stop("Paired metadata contains an unexpected or missing condition label.")
}

observed_condition_counts <- table(paired_metadata$Condition)

if (!identical(as.integer(observed_condition_counts[condition_levels]), as.integer(expected_condition_counts[condition_levels]))) {
  stop("Unexpected paired58 condition composition. Observed Control = ", observed_condition_counts[["Control"]], ", Treated = ",
       observed_condition_counts[["Treated"]], "; expected Control = 27, Treated = 31.")
}

attached <- attach_metadata_to_model(mofa_model, paired_metadata)
mofa_model <- attached$model
paired_metadata <- attached$metadata

cat("Views: ", paste(model_views, collapse = ", "), "\n", sep = "")
cat("Samples: ", nrow(paired_metadata), " (27 Control; 31 Treated)\n", sep = "")
cat("Factors: ", MOFA2::get_dimensions(mofa_model)$K, "\n\n", sep = "")

# =========================================================
# 5. Variance explained
# =========================================================

variance_explained <- MOFA2::get_variance_explained(mofa_model, as.data.frame = TRUE)

variance_per_factor <- variance_explained$r2_per_factor %>% as_tibble() %>%
  transmute(Group = as.character(group), View = as.character(view), Factor = standardize_factor_name(factor),
            VarianceExplainedPercent = as.numeric(value))

variance_total <- variance_explained$r2_total %>% as_tibble() %>%
  transmute(Group = as.character(group), View = as.character(view), TotalVarianceExplainedPercent = as.numeric(value))

write_csv(variance_per_factor, file.path(downstream_dir, "mofa_variance_explained_per_factor_view.csv"))
write_csv(variance_total, file.path(downstream_dir, "mofa_total_variance_explained_by_view.csv"))

p_variance <- variance_per_factor %>% mutate(View = factor(View, levels = expected_views),
                                             Factor = factor(Factor, 
                                                             levels = unique(standardize_factor_name(MOFA2::factors_names(mofa_model))))) %>%
  ggplot(aes(x = View, y = Factor, fill = VarianceExplainedPercent)) + geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", VarianceExplainedPercent)), size = 3) + scale_fill_gradient(low = "white", high = "#3B4CC0") +
  labs(x = NULL, y = NULL, fill = "Variance (%)") + theme(panel.grid = element_blank())

save_svg("variance_explained_factor_by_view.svg", p_variance, 7.5, 4.5)

p_variance_total <- variance_total %>% mutate(View = factor(View, levels = expected_views)) %>%
  ggplot(aes(x = View, y = TotalVarianceExplainedPercent, fill = View)) + geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", TotalVarianceExplainedPercent)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c(RNA = "#009E73", Protein = "#E69F00", Phosphopeptide = "#CC79A7")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + labs(x = NULL, y = "Total variance explained (%)")

save_svg("total_variance_explained_by_view.svg", p_variance_total, 6.5, 4.5)

p_data_overview <- MOFA2::plot_data_overview(mofa_model)
save_svg("data_overview.svg", p_data_overview, 7.5, 4.5)

# =========================================================
# 6. Factor scores and association with treatment
# =========================================================

factor_scores <- MOFA2::get_factors(mofa_model, factors = "all", as.data.frame = TRUE) %>% as_tibble() %>%
  transmute(SampleID = as.character(sample), Group = as.character(group), Factor = standardize_factor_name(factor), value = as.numeric(value)) %>%
  left_join(paired_metadata %>% mutate(Condition = as.character(Condition)) %>% 
              select(SampleID, Condition, ConditionLabel, Block, WellIndex), by = "SampleID")

if (anyNA(factor_scores$Condition)) {
  stop("One or more factor-score rows could not be matched to condition metadata.")
}

factor_associations <- factor_scores %>% group_by(Factor) %>% group_modify(~ fit_factor_condition_association(.x)) %>% ungroup() %>%
  mutate(FDR = p.adjust(PValue, method = "BH"), ConditionAssociated = FDR <= condition_factor_fdr) %>% arrange(PValue)

prioritized_factors <- factor_associations %>% filter(ConditionAssociated) %>% pull(Factor)

if (length(prioritized_factors) == 0) {
  prioritized_factors <- factor_associations %>% slice_min(PValue, n = 1, with_ties = FALSE) %>% pull(Factor)
  
  warning("No factor passed condition FDR <= ", condition_factor_fdr, ". The lowest-P factor is retained for exploratory plots only.")
}

factor_associations <- factor_associations %>% mutate(PrioritizedForInterpretation = Factor %in% prioritized_factors)
factor_scores_oriented <- factor_scores %>% left_join(factor_associations %>%
                                                        select(Factor, OrientationSign, FDR, PrioritizedForInterpretation), by = "Factor") %>%
  mutate(OrientedValue = value * OrientationSign, Condition = factor(Condition, levels = condition_levels))

write_csv(factor_associations, file.path(downstream_dir, "mofa_factor_condition_associations.csv"))
write_csv(factor_scores_oriented, file.path(downstream_dir, "mofa_factor_scores_oriented_long.csv"))

p_factor_condition <- factor_scores_oriented %>% ggplot(aes(x = Condition, y = OrientedValue, colour = Condition, fill = Condition)) +
  geom_violin(alpha = 0.16, width = 0.9, trim = FALSE) + geom_boxplot(width = 0.24, outlier.shape = NA, alpha = 0.35) +
  geom_jitter(width = 0.10, height = 0, size = 1.4, alpha = 0.8) + facet_wrap(~ Factor, scales = "free_y", nrow = 1) + 
  scale_colour_manual(values = condition_colours) + scale_fill_manual(values = condition_colours) +
  labs(x = NULL, y = "Oriented factor score", colour = "Condition", fill = "Condition") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 35, hjust = 1))

save_svg("factor_scores_by_condition.svg", p_factor_condition, 11, 4.5)

p_factor_effect <- factor_associations %>% mutate(Factor = forcats::fct_reorder(Factor, OrientedTreatedMinusControl)) %>%
  ggplot(aes(x = OrientedTreatedMinusControl, y = Factor)) + geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_errorbarh(aes(xmin = OrientedCI95Low, xmax = OrientedCI95High), height = 0.18, colour = "grey35") + 
  geom_point(aes(colour = FDR <= condition_factor_fdr), size = 2.8) + 
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey45"), labels = c(`TRUE` = "FDR <= 0.05", `FALSE` = "FDR > 0.05")) +
  labs(x = "Oriented Treated - Control factor-score difference (95% CI)", y = NULL, colour = NULL) + theme(legend.position = "top")

save_svg("factor_treatment_effects.svg", p_factor_effect, 7.2, 4.5)

factor_pair <- factor_associations %>% arrange(PValue) %>% slice_head(n = 2) %>% pull(Factor)

if (length(factor_pair) == 2) {
  factor_pair_data <- factor_scores_oriented %>% filter(Factor %in% factor_pair) %>% select(SampleID, Condition, Factor, OrientedValue) %>%
    pivot_wider(names_from = Factor, values_from = OrientedValue)
  
  p_factor_pair <- ggplot(factor_pair_data, aes(x = .data[[factor_pair[1]]], y = .data[[factor_pair[2]]], colour = Condition)) +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") + geom_point(size = 2.4, alpha = 0.85) +
    scale_colour_manual(values = condition_colours) + labs(x = factor_pair[1], y = factor_pair[2], colour = "Condition") +
    theme(legend.position = "top")
  
  save_svg("top_condition_factors_scatter.svg", p_factor_pair, 5.8, 5.2)
}

# =========================================================
# 7. Orient and annotate feature weights
# =========================================================

weights <- MOFA2::get_weights(mofa_model, views = "all", factors = "all", as.data.frame = TRUE) %>% as_tibble() %>%
  transmute(MOFAFeature = as.character(feature), View = as.character(view), Factor = standardize_factor_name(factor),
            RawWeight = as.numeric(value)) %>% mutate(SourceFeatureID = sub("^[^:]+::", "", MOFAFeature))

bad_prefix <- !mapply(startsWith, weights$MOFAFeature, paste0(weights$View, "::"))

if (any(bad_prefix)) {
  stop("Unexpected MOFA feature prefix. Examples: ", paste(head(weights$MOFAFeature[bad_prefix], 20), collapse = ", "))
}

weights <- weights %>% left_join(factor_associations %>% select(Factor, OrientationSign, FDR, PrioritizedForInterpretation), by = "Factor") %>%
  mutate(OrientedWeight = RawWeight * OrientationSign)

feature_annotation <- make_feature_annotation(weights, protein_annotation_file, phosphopeptide_annotation_file)

weights_annotated <- weights %>% left_join(feature_annotation, by = c("View", "SourceFeatureID"))

unmatched_annotation <- weights_annotated %>% filter(is.na(FeatureLabel)) %>% distinct(View, SourceFeatureID)

if (nrow(unmatched_annotation) > 0) {
  stop("Feature annotation failed for ", nrow(unmatched_annotation), " model feature(s). Examples: ",
       paste(head(unmatched_annotation$SourceFeatureID, 20), collapse = ", "))
}

write_csv(weights_annotated, file.path(downstream_dir, "mofa_weights_oriented_annotated.csv"), na = "")

top_positive <- weights_annotated %>% filter(OrientedWeight > 0) %>% group_by(View, Factor) %>%
  slice_max(order_by = OrientedWeight, n = top_features_each_direction, with_ties = FALSE) %>% ungroup() %>% mutate(Direction = "Treated-associated")

top_negative <- weights_annotated %>% filter(OrientedWeight < 0) %>% group_by(View, Factor) %>% 
  slice_min(order_by = OrientedWeight, n = top_features_each_direction, with_ties = FALSE) %>% ungroup() %>%
  mutate(Direction = "Control-associated")

top_weights <- bind_rows(top_positive, top_negative) %>% arrange(Factor, View, desc(abs(OrientedWeight)))

write_csv(top_weights, file.path(downstream_dir, "mofa_top_weights_by_view_factor.csv"), na = "")

top_weights_plot_data <- top_weights %>% filter(Factor %in% prioritized_factors) %>% 
  mutate(View = factor(View, levels = expected_views), Factor = factor(Factor, levels = prioritized_factors),
         FeaturePlotID = paste(View, Factor, SourceFeatureID, sep = "___"))

if (nrow(top_weights_plot_data) > 0) {
  feature_labels <- setNames(top_weights_plot_data$FeatureLabel, top_weights_plot_data$FeaturePlotID)
  
  p_top_weights <- top_weights_plot_data %>% ggplot(aes(x = OrientedWeight, y = reorder(FeaturePlotID, OrientedWeight), fill = Direction)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) + geom_col(width = 0.72) +
    facet_grid(View ~ Factor, scales = "free_y", space = "free_y") + scale_fill_manual(values = direction_colours) +
    scale_y_discrete(labels = function(x) unname(feature_labels[x])) + labs(x = "Oriented MOFA+ weight", y = NULL, fill = NULL) +
    theme(legend.position = "top", axis.text.y = element_text(size = 7), panel.grid.major.y = element_blank())
  
  plot_height <- max(7, 2.4 * length(expected_views))
  plot_width <- max(8, 5 * length(prioritized_factors))
  
  save_svg("top_oriented_weights_condition_factors.svg", p_top_weights, plot_width, plot_height)
}

# =========================================================
# 8. Compact summary and reproducibility records
# =========================================================

summary_report <- tibble(Item = c("Model file", "Samples", "Control samples", "Treated samples", "Factors", "Views",
                                  "Condition factor FDR threshold", "Prioritized factors"),
                         Value = c(model_file, as.character(nrow(paired_metadata)), as.character(expected_condition_counts[["Control"]]),
                                   as.character(expected_condition_counts[["Treated"]]), as.character(MOFA2::get_dimensions(mofa_model)$K),
                                   paste(model_views, collapse = ";"), as.character(condition_factor_fdr), paste(prioritized_factors, collapse = ";")))

write_csv(summary_report, file.path(downstream_dir, "mofa_exploratory_summary.csv"))
write_lines(capture.output(sessionInfo()), file.path(downstream_dir, "mofa_exploratory_sessionInfo.txt"))

cat("\nPaired58 MOFA+ exploratory analysis complete.\n")
cat("Output directory:\n", downstream_dir, "\n", sep = "")
cat("Prioritized factor(s) for interpretation: ", paste(prioritized_factors, collapse = ", "), "\n", sep = "")