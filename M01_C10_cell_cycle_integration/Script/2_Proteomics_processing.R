# =========================================================
# Proteomics preprocessing for NanoSPLITS C10 RO-3306 experiment
# =========================================================
# Purpose:
#   Starting from:
#     C10Treated_Protein_intensities.tsv
#
#   Generate:
#     1. Shared log2-transformed, proDA-style median-normalized protein matrices
#     2. MOFA+-ready protein matrices
#        - all 68 author-QC proteomics cells
#        - strict paired RNA/protein-QC benchmark cells (expected n = 58)
#     3. SCOT+-ready complete protein matrices
#        - missing values imputed only in the SCOT+ branch
#     4. SCOT+-ready protein PCA embeddings
#     5. Protein sample QC, feature QC, feature-selection, and preprocessing metadata
#
# Modelling decisions:
#   - Mouse proteins are identified using "_MOUSE" in Entry.Name.
#   - Zero intensity means not detected and is converted to NA
#   - Positive intensities are log2 transformed
#   - On the log2 scale, each sample is normalized by subtracting its median deviation 
#     from the across-sample protein row mean
#   - Protein features must be detected in at least 50% of cells in the
#     relevant cohort, matching the NanoSPLITS heatmap logic:
#       missing <= 34 among 68 cells
#   - The top variable proteins are selected separately in the all-68 and
#     paired-58 cohorts
#   - MOFA+ matrices retain NA values
#   - MOFA+ matrices are feature-wise z-scored using observed values only
#   - SCOT+ uses the same z-scored feature space, with missing values replaced
#     by 0, which is the observed feature mean after z-scoring
#   - MOFA+ outputs are features x samples
#   - SCOT+ outputs are samples x features or samples x PCs
#

library(tidyverse)
library(readr)
library(rstudioapi)

# =========================================================
# 1. Path setup
# =========================================================

get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && !is.null(ctx$path) && nzchar(ctx$path)) {
      return(normalizePath(dirname(ctx$path), winslash = "/", mustWork = TRUE))
    }
  }
  
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    return(normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE))
  }
  
  warning("Could not determine the script directory from RStudio or Rscript. Using the current working directory.")
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
raw_data_dir <- "C:/Users/49152/Downloads/2nd_paper/Github_data"
processed_dir <- file.path(script_dir, "processed_data")
shared_dir <- file.path(processed_dir, "shared")
mofa_dir <- file.path(processed_dir, "mofa_input")
scot_dir <- file.path(processed_dir, "scot_input")
qc_dir <- file.path(processed_dir, "qc")

walk(c(processed_dir, shared_dir, mofa_dir, scot_dir, qc_dir), dir.create, showWarnings = FALSE, recursive = TRUE)

protein_file <- file.path(raw_data_dir, "C10Treated_Protein_intensities.tsv")
metadata_file <- file.path(processed_dir, "metadata_master_samples.csv")
protein_feature_metadata_file <- file.path(processed_dir, "metadata_protein_features.csv")

if (!file.exists(protein_file)) stop("Protein input file does not exist:\n", protein_file)
if (!file.exists(metadata_file)) {
  stop("Master metadata file does not exist:\n", metadata_file, "\nRun the metadata preparation script first.")
}

cat("Script directory:   ", script_dir, "\n")
cat("Raw data directory: ", raw_data_dir, "\n")
cat("Output directory:   ", processed_dir, "\n\n")

# =========================================================
# 2. User-adjustable parameters
# =========================================================

min_detected_fraction <- 0.50
n_variable_proteins <- 2000
n_scot_pcs <- 30
expected_n_protein_qc <- 68
expected_n_paired <- 58

# =========================================================
# 3. Helper functions
# =========================================================

safe_numeric_matrix_allow_na <- function(df, object_name = "matrix") {
  original <- as.matrix(df)
  original_na <- is.na(original)
  numeric_values <- suppressWarnings(as.numeric(original))
  mat <- matrix(numeric_values, nrow = nrow(df), ncol = ncol(df), dimnames = list(rownames(df), colnames(df)))
  introduced_na <- is.na(mat) & !original_na
  
  if (any(introduced_na)) {
    stop(object_name, " contains non-numeric values in presumed sample columns. ", sum(introduced_na),
         " value(s) were converted to NA.")
  }
  if (any(!is.finite(mat[!is.na(mat)]))) {
    stop(object_name, " contains non-finite non-missing values.")
  }
  mat
}

validate_unique_ids <- function(ids, object_name) {
  if (anyNA(ids)) stop(object_name, " contains missing identifiers.")
  if (any(!nzchar(ids))) stop(object_name, " contains blank identifiers.")
  duplicated_ids <- unique(ids[duplicated(ids)])
  if (length(duplicated_ids) > 0) {
    stop(object_name, " contains duplicated identifiers: ", paste(duplicated_ids, collapse = ", "))
  }
  invisible(TRUE)
}

order_metadata_to_samples <- function(metadata, sample_ids) {
  index <- match(sample_ids, metadata$SampleID)
  if (anyNA(index)) {
    stop("Metadata could not be matched for sample(s): ", paste(sample_ids[is.na(index)], collapse = ", "))
  }
  output <- metadata[index, , drop = FALSE]
  if (!identical(output$SampleID, sample_ids)) {
    stop("Internal error while ordering metadata to sample IDs.")
  }
  output
}

write_matrix <- function(mat, file, id_column) {
  as.data.frame(mat, check.names = FALSE) %>%
    rownames_to_column(var = id_column) %>%
    write_csv(file, na = "")
}

proda_median_normalize <- function(mat) {
  
  # Mean abundance of each protein across observed samples.
  row_average <- rowMeans(mat, na.rm = TRUE)
  
  # For each sample, calculate the median deviation from the
  # across-sample protein mean. This reproduces the core calculation
  # used by proDA::median_normalization().
  sample_offsets <- vapply(seq_len(ncol(mat)), function(j) {
      median(mat[, j] - row_average, na.rm = TRUE)
    }, numeric(1))
  
  names(sample_offsets) <- colnames(mat)
  
  if (any(!is.finite(sample_offsets))) {
    bad_samples <- names(sample_offsets)[!is.finite(sample_offsets)]
    
    stop("Cannot median-normalize sample(s) with no usable detected values: ", paste(bad_samples, collapse = ", "))
  }
  
  normalized_mat <- sweep(mat, 2, sample_offsets, FUN = "-")
  
  if (any(!is.finite(normalized_mat[!is.na(normalized_mat)]))) {
    stop("Median normalization produced non-finite observed values.")
  }
  
  list(matrix = normalized_mat, sample_offsets = sample_offsets)
}

calculate_protein_qc <- function(mat, cohort_name, min_detected_fraction) {
  detected_n <- rowSums(!is.na(mat))
  detected_fraction <- rowMeans(!is.na(mat))
  observed_variance <- apply(mat, 1, var, na.rm = TRUE)
  
  tibble(ProteinRowID = rownames(mat), Cohort = cohort_name, DetectedN = detected_n, MissingN = ncol(mat) - detected_n,
         DetectedFraction = detected_fraction, MissingFraction = 1 - detected_fraction,
         MeanObserved = rowMeans(mat, na.rm = TRUE), MedianObserved = apply(mat, 1, median, na.rm = TRUE),
         VarianceObserved = observed_variance, PassDetectionFilter = DetectedFraction >= min_detected_fraction,
         PassVarianceFilter = is.finite(VarianceObserved) & VarianceObserved > 0,
         PassCombinedFilter = PassDetectionFilter & PassVarianceFilter)
}

select_top_variable_proteins <- function(mat, n_features) {
  observed_variance <- apply(mat, 1, var, na.rm = TRUE)
  valid <- is.finite(observed_variance) & observed_variance > 0
  observed_variance <- observed_variance[valid]
  if (length(observed_variance) == 0) {
    stop("No variable proteins are available for feature selection.")
  }
  n_use <- min(n_features, length(observed_variance))
  if (n_use < n_features) {
    warning("Requested ", n_features, " variable proteins, but only ", n_use, " proteins were available after filtering.")
  }
  selected <- names(sort(observed_variance, decreasing = TRUE))[seq_len(n_use)]
  list(selected = selected, variance = observed_variance)
}

zscore_rows_preserve_na <- function(mat) {
  observed_n <- rowSums(!is.na(mat))
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- apply(mat, 1, sd, na.rm = TRUE)
  keep <- observed_n >= 2 & is.finite(row_means) & is.finite(row_sds) & row_sds > 0
  dropped_features <- rownames(mat)[!keep]
  mat_keep <- mat[keep, , drop = FALSE]
  means_keep <- row_means[keep]
  sds_keep <- row_sds[keep]
  if (nrow(mat_keep) == 0) stop("No proteins remain after scaling filters.")
  z <- sweep(mat_keep, 1, means_keep, FUN = "-")
  z <- sweep(z, 1, sds_keep, FUN = "/")
  if (any(!is.finite(z[!is.na(z)]))) {
    stop("Non-finite observed values remain after protein z-scoring.")
  }
  list(matrix = z, means = means_keep, sds = sds_keep, dropped_features = dropped_features)
}

impute_scot_feature_means <- function(z_mat) {
  imputed <- z_mat
  imputed[is.na(imputed)] <- 0
  if (anyNA(imputed) || any(!is.finite(imputed))) {
    stop("SCOT+ protein matrix is not complete and finite after imputation.")
  }
  imputed
}

make_scot_pca <- function(samples_by_features, n_pcs = 30, center = TRUE, scale_features = FALSE) {
  if (anyNA(samples_by_features) || any(!is.finite(samples_by_features))) {
    stop("SCOT+ PCA input contains missing or non-finite values.")
  }
  n_pcs_use <- min(n_pcs, nrow(samples_by_features) - 1, ncol(samples_by_features))
  if (n_pcs_use < 2) stop("Too few samples or features to calculate PCA.")
  pca <- prcomp(samples_by_features, center = center, scale. = scale_features, rank. = n_pcs_use)
  embedding <- pca$x[, seq_len(n_pcs_use), drop = FALSE]
  all_variance_fraction <- pca$sdev^2 / sum(pca$sdev^2)
  variance_fraction <- all_variance_fraction[seq_len(n_pcs_use)]
  variance_explained <- tibble(PC = colnames(embedding), VarianceExplained = variance_fraction,
                               CumulativeVarianceExplained = cumsum(variance_fraction))
  list(embedding = embedding, variance_explained = variance_explained, rotation = pca$rotation)
}

safe_cor <- function(x, y, method = "pearson") {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3) return(NA_real_)
  sd_x <- sd(x)
  sd_y <- sd(y)
  if (!is.finite(sd_x) || !is.finite(sd_y) || sd_x == 0 || sd_y == 0) {
    return(NA_real_)
  }
  cor(x, y, method = method)
}

assert_mofa_export <- function(mat, metadata) {
  if (!identical(colnames(mat), metadata$SampleID)) {
    stop("MOFA+ protein matrix columns do not match metadata order.")
  }
  if (any(!is.finite(mat[!is.na(mat)]))) {
    stop("MOFA+ protein matrix contains non-finite observed values.")
  }
  invisible(TRUE)
}

assert_scot_export <- function(mat, metadata) {
  if (!identical(rownames(mat), metadata$SampleID)) {
    stop("SCOT+ protein matrix rows do not match metadata order.")
  }
  if (anyNA(mat) || any(!is.finite(mat))) {
    stop("SCOT+ protein matrix must be complete and finite.")
  }
  invisible(TRUE)
}

# =========================================================
# 4. Load and validate master metadata
# =========================================================

cat("Loading master metadata...\n")
master_metadata <- read_csv(metadata_file, show_col_types = FALSE)

required_meta_cols <- c("SampleID", "TreatmentCode", "Condition", "ConditionLabel", "Block", "WellIndex",
                        "Protein_available", "Protein_detected_n", "Protein_missing_n", "Protein_missing_fraction",
                        "Protein_QC_pass_author", "ProteinQC68", "Paired_RNA_ProteinQC")
missing_meta_cols <- setdiff(required_meta_cols, colnames(master_metadata))
if (length(missing_meta_cols) > 0) {
  stop("metadata_master_samples.csv is missing required column(s): ", paste(missing_meta_cols, collapse = ", "))
}

validate_unique_ids(master_metadata$SampleID, "Master metadata SampleID column")
for (column_name in c("Protein_available", "Protein_QC_pass_author", "ProteinQC68", "Paired_RNA_ProteinQC")) {
  if (!is.logical(master_metadata[[column_name]])) {
    stop(column_name, " must be imported as a logical TRUE/FALSE column.")
  }
}

protein_qc_sample_ids <- master_metadata %>% filter(ProteinQC68) %>% pull(SampleID)
paired_sample_ids <- master_metadata %>% filter(Paired_RNA_ProteinQC) %>% pull(SampleID)

validate_unique_ids(protein_qc_sample_ids, "Protein-QC sample-ID set")
validate_unique_ids(paired_sample_ids, "Paired RNA/protein-QC sample-ID set")

if (!all(paired_sample_ids %in% protein_qc_sample_ids)) {
  stop("The paired RNA/protein-QC sample set is not a subset of the protein-QC sample set.")
}
if (length(protein_qc_sample_ids) != expected_n_protein_qc) {
  stop("Expected ", expected_n_protein_qc, " protein-QC samples but found ", length(protein_qc_sample_ids), ".")
}
if (length(paired_sample_ids) != expected_n_paired) {
  stop("Expected ", expected_n_paired, " paired samples but found ", length(paired_sample_ids), ".")
}

cat("Author-QC protein samples:            ", length(protein_qc_sample_ids), "\n")
cat("Strict paired RNA/protein-QC samples: ", length(paired_sample_ids), "\n\n")

# =========================================================
# 5. Load and validate protein matrix
# =========================================================

cat("Loading protein intensity matrix...\n")
protein_raw <- read_tsv(protein_file, show_col_types = FALSE, progress = FALSE, name_repair = "minimal")

if (anyDuplicated(colnames(protein_raw))) {
  duplicated_columns <- unique(colnames(protein_raw)[duplicated(colnames(protein_raw))])
  
  stop("Duplicated column names detected in protein input: ", paste(duplicated_columns, collapse = ", "))
}
required_protein_annot <- c("Protein", "Protein.ID", "Entry.Name", "Gene")
missing_protein_annot <- setdiff(required_protein_annot, colnames(protein_raw))
if (length(missing_protein_annot) > 0) {
  stop("Protein file is missing expected annotation column(s): ", paste(missing_protein_annot, collapse = ", "))
}

# Define protein sample columns from the metadata rather than assuming that every non-annotation column
# in the TSV is a sample column.
protein_sample_cols <- master_metadata %>%
  filter(Protein_available) %>%
  pull(SampleID)

validate_unique_ids(protein_sample_cols, "Protein matrix sample columns defined by metadata")

missing_protein_sample_cols <- setdiff(protein_sample_cols, colnames(protein_raw))

if (length(missing_protein_sample_cols) > 0) {
  stop("Protein samples marked as available in metadata are missing from the ",
       "protein intensity file: ", paste(missing_protein_sample_cols, collapse = ", "))
}

protein_mouse <- protein_raw %>%
  mutate(ProteinRowID = sprintf("PROT_%06d", row_number())) %>%
  filter(!is.na(Entry.Name), grepl("_MOUSE", Entry.Name, fixed = TRUE))

if (nrow(protein_mouse) == 0) {
  stop("No mouse protein rows were detected using Entry.Name.")
}

validate_unique_ids(protein_mouse$ProteinRowID, "Mouse protein row IDs")

duplicated_entry_names <- unique(protein_mouse$Entry.Name[duplicated(protein_mouse$Entry.Name)])

if (length(duplicated_entry_names) > 0) {
  stop("Duplicated mouse Entry.Name identifiers were detected. ",
       "Do not silently treat these as independent protein features. Examples: ",
       paste(head(duplicated_entry_names, 20), collapse = ", "))
}

protein_mat_raw <- protein_mouse %>%
  select(all_of(protein_sample_cols)) %>%
  safe_numeric_matrix_allow_na(object_name = "Raw protein intensity matrix")
rownames(protein_mat_raw) <- protein_mouse$ProteinRowID

if (any(protein_mat_raw < 0, na.rm = TRUE)) {
  stop("Raw protein intensity matrix contains negative values.")
}
protein_mat_raw[protein_mat_raw == 0] <- NA_real_

missing_in_matrix <- setdiff(protein_qc_sample_ids, colnames(protein_mat_raw))
if (length(missing_in_matrix) > 0) {
  stop("Protein-QC metadata samples missing from the protein matrix: ", paste(missing_in_matrix, collapse = ", "))
}

# =========================================================
# 6. Validate sample QC against metadata
# =========================================================

protein_detected_recalculated <- vapply(seq_len(ncol(protein_mat_raw)), function(j) {
    detected <- !is.na(protein_mat_raw[, j])
    
    dplyr::n_distinct(protein_mouse$Entry.Name[detected], na.rm = TRUE)
  }, integer(1))

names(protein_detected_recalculated) <- colnames(protein_mat_raw)

protein_qc_validation <- master_metadata %>%
  filter(Protein_available) %>%
  select(SampleID, Condition, Protein_detected_n_metadata = Protein_detected_n,
         Protein_QC_pass_author_metadata = Protein_QC_pass_author) %>%
  mutate(Protein_detected_n_recalculated = unname(protein_detected_recalculated[SampleID]),
         Protein_QC_pass_author_recalculated = case_when(
           Condition == "Control" & Protein_detected_n_recalculated >= 2000 ~ TRUE,
           Condition == "Treated" & Protein_detected_n_recalculated >= 2500 ~ TRUE, TRUE ~ FALSE))

if (!isTRUE(all.equal(unname(protein_qc_validation$Protein_detected_n_recalculated),
                      unname(protein_qc_validation$Protein_detected_n_metadata), tolerance = 0, check.attributes = FALSE))) {
  stop("Recalculated protein detection counts do not match metadata.")
}
if (!identical(protein_qc_validation$Protein_QC_pass_author_recalculated,
               protein_qc_validation$Protein_QC_pass_author_metadata)) {
  stop("Recalculated author protein-QC decisions do not match metadata.")
}
write_csv(protein_qc_validation, file.path(qc_dir, "protein_sample_qc_validation.csv"))

protein_raw_all68 <- protein_mat_raw[, protein_qc_sample_ids, drop = FALSE]
protein_raw_paired58 <- protein_mat_raw[, paired_sample_ids, drop = FALSE]

# =========================================================
# 7. Log2 transformation and proDA-style median normalization
# =========================================================

protein_log2_all68 <- log2(protein_raw_all68)
protein_log2_paired58 <- log2(protein_raw_paired58)

normalized_all68 <- proda_median_normalize(protein_log2_all68)
normalized_paired58 <- proda_median_normalize(protein_log2_paired58)
protein_norm_all68 <- normalized_all68$matrix
protein_norm_paired58 <- normalized_paired58$matrix

sample_qc_all68 <- tibble(SampleID = colnames(protein_norm_all68), DetectedProteins = colSums(!is.na(protein_norm_all68)),
                          MissingProteins = colSums(is.na(protein_norm_all68)),
                          MissingFraction = colMeans(is.na(protein_norm_all68)),
                          RawLog2Median = apply(protein_log2_all68, 2, median, na.rm = TRUE),
                          NormalizationOffset = unname(normalized_all68$sample_offsets),
                          NormalizedLog2Median = apply(protein_norm_all68, 2, median, na.rm = TRUE)) %>%
  left_join(master_metadata %>%
              select(SampleID, Condition, ConditionLabel, Block, WellIndex, Paired_RNA_ProteinQC), by = "SampleID")

if (anyNA(sample_qc_all68$Condition)) {
  stop("One or more protein-QC samples lack condition metadata.")
}
write_csv(sample_qc_all68, file.path(qc_dir, "protein_sample_qc_metrics_all68.csv"))
write_matrix(protein_norm_all68, file.path(shared_dir, "protein_log2_proda_median_normalized_all68.csv"), "ProteinRowID")
write_matrix(protein_norm_paired58, file.path(shared_dir, "protein_log2_proda_median_normalized_paired58.csv"), "ProteinRowID")

# =========================================================
# 8. Cohort-specific protein QC and filtering
# =========================================================

protein_qc_all68 <- calculate_protein_qc(protein_norm_all68, "ProteinQC68", min_detected_fraction)
protein_qc_paired58 <- calculate_protein_qc(protein_norm_paired58, "Paired58", min_detected_fraction)

proteins_keep_all68 <- protein_qc_all68 %>%
  filter(PassCombinedFilter) %>% pull(ProteinRowID)
proteins_keep_paired58 <- protein_qc_paired58 %>%
  filter(PassCombinedFilter) %>% pull(ProteinRowID)

if (length(proteins_keep_all68) == 0) stop("No proteins passed the all-68 filter.")
if (length(proteins_keep_paired58) == 0) stop("No proteins passed the paired-58 filter.")

protein_filtered_all68 <- protein_norm_all68[proteins_keep_all68, protein_qc_sample_ids, drop = FALSE]
protein_filtered_paired58 <- protein_norm_paired58[proteins_keep_paired58, paired_sample_ids, drop = FALSE]

write_csv(bind_rows(protein_qc_all68, protein_qc_paired58), file.path(qc_dir, "protein_feature_qc_metrics_by_cohort.csv"))
write_lines(proteins_keep_all68, file.path(shared_dir, "protein_features_after_filter_all68.txt"))
write_lines(proteins_keep_paired58, file.path(shared_dir, "protein_features_after_filter_paired58.txt"))

# =========================================================
# 9. Cohort-specific variable-protein selection
# =========================================================

selection_all68 <- select_top_variable_proteins(protein_filtered_all68, n_variable_proteins)
selection_paired58 <- select_top_variable_proteins(protein_filtered_paired58, n_variable_proteins)

write_lines(selection_all68$selected, file.path(shared_dir, "protein_variable_features_all68.txt"))
write_lines(selection_paired58$selected, file.path(shared_dir, "protein_variable_features_paired58.txt"))

all_feature_ids <- union(rownames(protein_filtered_all68), rownames(protein_filtered_paired58))
variable_protein_table <- tibble(ProteinRowID = all_feature_ids,
                                 Variance_all68 = selection_all68$variance[
                                   match(ProteinRowID, names(selection_all68$variance))],
                                 Variance_paired58 = selection_paired58$variance[
                                   match(ProteinRowID, names(selection_paired58$variance))],
                                 Selected_all68 = ProteinRowID %in% selection_all68$selected,
                                 Selected_paired58 = ProteinRowID %in% selection_paired58$selected) %>%
  left_join(protein_mouse %>%
              select(ProteinRowID, all_of(required_protein_annot)), by = "ProteinRowID") %>%
  arrange(desc(Variance_all68))
write_csv(variable_protein_table, file.path(qc_dir, "protein_variable_feature_selection_table.csv"))

# =========================================================
# 10. Feature-wise scaling while preserving missing values
# =========================================================

protein_all68_selected <- protein_filtered_all68[selection_all68$selected, protein_qc_sample_ids, drop = FALSE]
protein_paired58_selected <- protein_filtered_paired58[selection_paired58$selected, paired_sample_ids, drop = FALSE]

scaled_all68 <- zscore_rows_preserve_na(protein_all68_selected)
scaled_paired58 <- zscore_rows_preserve_na(protein_paired58_selected)
protein_all68_z <- scaled_all68$matrix
protein_paired58_z <- scaled_paired58$matrix

write_lines(scaled_all68$dropped_features, file.path(qc_dir, "protein_all68_features_dropped_during_scaling.txt"))
write_lines(scaled_paired58$dropped_features, file.path(qc_dir, "protein_paired58_features_dropped_during_scaling.txt"))
write_csv(tibble(ProteinRowID = names(scaled_all68$means), ObservedMean = unname(scaled_all68$means),
                 ObservedSD = unname(scaled_all68$sds)), file.path(qc_dir, "protein_scaling_parameters_all68.csv"))
write_csv(tibble(ProteinRowID = names(scaled_paired58$means), ObservedMean = unname(scaled_paired58$means),
                 ObservedSD = unname(scaled_paired58$sds)), file.path(qc_dir, "protein_scaling_parameters_paired58.csv"))

metadata_all68 <- order_metadata_to_samples(master_metadata, colnames(protein_all68_z))
metadata_paired58 <- order_metadata_to_samples(master_metadata, colnames(protein_paired58_z))
assert_mofa_export(protein_all68_z, metadata_all68)
assert_mofa_export(protein_paired58_z, metadata_paired58)

# =========================================================
# 11. MOFA+ exports: missing values retained
# =========================================================

write_matrix(protein_all68_z, file.path(mofa_dir, "mofa_protein_all68_z_with_na.csv"), "ProteinRowID")
write_matrix(protein_paired58_z, file.path(mofa_dir, "mofa_protein_paired58_z_with_na.csv"), "ProteinRowID")
write_matrix(protein_all68_selected, file.path(mofa_dir, "mofa_protein_all68_log2_proda_median_normalized_with_na.csv"),
             "ProteinRowID")
write_matrix(protein_paired58_selected, file.path(mofa_dir, "mofa_protein_paired58_log2_proda_median_normalized_with_na.csv"),
             "ProteinRowID")
write_csv(metadata_all68, file.path(mofa_dir, "mofa_protein_all68_metadata.csv"))
write_csv(metadata_paired58, file.path(mofa_dir, "mofa_protein_paired58_metadata.csv"))

# =========================================================
# 12. SCOT+ exports: imputation only for alignment geometry
# =========================================================

protein_all68_z_imputed <- impute_scot_feature_means(protein_all68_z)
protein_paired58_z_imputed <- impute_scot_feature_means(protein_paired58_z)
scot_protein_all68 <- t(protein_all68_z_imputed)
scot_protein_paired58 <- t(protein_paired58_z_imputed)

assert_scot_export(scot_protein_all68, metadata_all68)
assert_scot_export(scot_protein_paired58, metadata_paired58)

write_matrix(scot_protein_all68, file.path(scot_dir, "scot_protein_all68_samples_by_features_z_mean_imputed.csv"),
             "SampleID")
write_matrix(scot_protein_paired58, file.path(scot_dir, "scot_protein_paired58_samples_by_features_z_mean_imputed.csv"),
             "SampleID")
write_csv(metadata_all68, file.path(scot_dir, "scot_protein_all68_metadata.csv"))
write_csv(metadata_paired58, file.path(scot_dir, "scot_protein_paired58_metadata.csv"))
scot_imputation_summary <- bind_rows(tibble(Cohort = "all68", ProteinRowID = rownames(protein_all68_z),
                                            MissingBeforeImputation = rowSums(is.na(protein_all68_z))),
                                     tibble(Cohort = "paired58", ProteinRowID = rownames(protein_paired58_z),
                                            MissingBeforeImputation = rowSums(is.na(protein_paired58_z)))) %>%
  mutate(ImputationMethod = "zero_after_feature_zscore")
write_csv(scot_imputation_summary, file.path(qc_dir, "protein_scot_imputation_summary.csv"))

# =========================================================
# 13. SCOT+-ready PCA embeddings
# =========================================================

cat("Computing SCOT+ protein PCA embeddings...\n")
scot_pca_all68 <- make_scot_pca(scot_protein_all68, n_pcs = n_scot_pcs)
scot_pca_paired58 <- make_scot_pca(scot_protein_paired58, n_pcs = n_scot_pcs)

if (!identical(rownames(scot_pca_all68$embedding), metadata_all68$SampleID)) {
  stop("All-68 protein PCA sample order does not match metadata.")
}
if (!identical(rownames(scot_pca_paired58$embedding), metadata_paired58$SampleID)) {
  stop("Paired-58 protein PCA sample order does not match metadata.")
}

write_matrix(scot_pca_all68$embedding, file.path(scot_dir, "scot_protein_all68_pca_embedding.csv"), "SampleID")
write_matrix(scot_pca_paired58$embedding, file.path(scot_dir, "scot_protein_paired58_pca_embedding.csv"), "SampleID")
write_csv(scot_pca_all68$variance_explained, file.path(scot_dir, "scot_protein_all68_pca_variance_explained.csv"))
write_csv(scot_pca_paired58$variance_explained, file.path(scot_dir, "scot_protein_paired58_pca_variance_explained.csv"))

# =========================================================
# 14. Paired protein PCA diagnostics
# =========================================================

diagnostic_scores <- as.data.frame(scot_pca_paired58$embedding, check.names = FALSE) %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(metadata_paired58 %>%
              select(SampleID, Condition, ConditionLabel, Block, WellIndex, Protein_detected_n, Protein_missing_fraction),
            by = "SampleID")

if (anyNA(diagnostic_scores$Condition)) {
  stop("Protein PCA scores could not be fully matched to metadata.")
}
write_csv(diagnostic_scores, file.path(qc_dir, "protein_paired58_pca_scores.csv"))
write_csv(scot_pca_paired58$variance_explained, file.path(qc_dir, "protein_paired58_pca_variance_explained.csv"))

available_pcs <- intersect(paste0("PC", 1:5), colnames(diagnostic_scores))
technical_metrics <- c("Protein_detected_n", "Protein_missing_fraction")
technical_correlations <- map_dfr(available_pcs, function(pc_name) {
  map_dfr(technical_metrics, function(metric_name) {
    tibble(PC = pc_name, Metric = metric_name, Correlation = safe_cor(diagnostic_scores[[pc_name]],
                                                                      diagnostic_scores[[metric_name]]))
    })})
write_csv(technical_correlations, file.path(qc_dir, "protein_pca_technical_metric_correlations.csv"))

condition_pca_summary <- diagnostic_scores %>%
  group_by(Condition) %>%
  summarise(n = n(), across(all_of(available_pcs), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
            .groups = "drop")
block_pca_summary <- diagnostic_scores %>%
  group_by(Block) %>%
  summarise(n = n(), across(all_of(available_pcs), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
            .groups = "drop")
write_csv(condition_pca_summary, file.path(qc_dir, "protein_pca_scores_summary_by_condition.csv"))
write_csv(block_pca_summary, file.path(qc_dir, "protein_pca_scores_summary_by_block.csv"))

if (interactive() && all(c("PC1", "PC2") %in% colnames(diagnostic_scores))) {
  plot(diagnostic_scores$PC1, diagnostic_scores$PC2, col = ifelse(diagnostic_scores$Condition == "Control", "blue", "red"),
       pch = 16, xlab = "PC1", ylab = "PC2", main = "Protein PCA: paired SCOT+ input")
  legend("topright", legend = c("Control", "Treated"), col = c("blue", "red"), pch = 16)
}

# =========================================================
# 15. Feature metadata exports
# =========================================================

protein_feature_metadata_current <- protein_mouse %>%
  select(ProteinRowID, all_of(required_protein_annot))

if (file.exists(protein_feature_metadata_file)) {
  protein_feature_metadata_saved <- read_csv(protein_feature_metadata_file, show_col_types = FALSE)
  if ("ProteinRowID" %in% colnames(protein_feature_metadata_saved)) {
    if (!identical(protein_feature_metadata_current$ProteinRowID, protein_feature_metadata_saved$ProteinRowID)) {
      stop("ProteinRowID values do not match the metadata preparation output.")
    }
  }
}

selected_feature_metadata <- protein_feature_metadata_current %>%
  mutate(Selected_all68 = ProteinRowID %in% rownames(protein_all68_z),
         Selected_paired58 = ProteinRowID %in% rownames(protein_paired58_z)) %>%
  left_join(protein_qc_all68 %>%
              select(ProteinRowID, DetectedFraction_all68 = DetectedFraction, VarianceObserved_all68 = VarianceObserved),
            by = "ProteinRowID") %>%
  left_join(protein_qc_paired58 %>%
              select(ProteinRowID, DetectedFraction_paired58 = DetectedFraction,
                     VarianceObserved_paired58 = VarianceObserved), by = "ProteinRowID")
write_csv(selected_feature_metadata, file.path(shared_dir, "protein_feature_metadata_preprocessed.csv"))

# =========================================================
# 16. Final checks and summary
# =========================================================

stopifnot(identical(colnames(protein_all68_z), metadata_all68$SampleID),
          identical(colnames(protein_paired58_z), metadata_paired58$SampleID),
          identical(rownames(scot_protein_all68), metadata_all68$SampleID),
          identical(rownames(scot_protein_paired58), metadata_paired58$SampleID), all(is.finite(scot_protein_all68)),
          all(is.finite(scot_protein_paired58)), all(is.finite(scot_pca_all68$embedding)),
          all(is.finite(scot_pca_paired58$embedding)))

summary_report <- tibble(Item = c("Input protein samples", "Input mouse protein rows", "Author-QC protein samples",
                                  "Paired benchmark samples", "Minimum detected fraction", "Variable proteins requested",
                                  "Proteins after all68 filter", "Proteins after paired58 filter", "All68 MOFA features",
                                  "Paired58 MOFA features", "All68 MOFA missing fraction", "Paired58 MOFA missing fraction",
                                  "SCOT all68 samples", "SCOT paired58 samples", "SCOT PCs exported"),
                         Value = as.character(c(ncol(protein_mat_raw), nrow(protein_mat_raw), length(protein_qc_sample_ids),
                                                length(paired_sample_ids), min_detected_fraction, n_variable_proteins,
                                                nrow(protein_filtered_all68), nrow(protein_filtered_paired58),
                                                nrow(protein_all68_z), nrow(protein_paired58_z),
                                                mean(is.na(protein_all68_z)), mean(is.na(protein_paired58_z)),
                                                nrow(scot_protein_all68), nrow(scot_protein_paired58),
                                                ncol(scot_pca_all68$embedding))))
write_csv(summary_report, file.path(qc_dir, "protein_preprocessing_summary.csv"))

parameter_report <- tibble(Parameter = c("min_detected_fraction", "n_variable_proteins", "n_scot_pcs",
                                         "zero_represents_missing", "normalization", "mofa_missing_values",
                                         "scot_imputation"), Value = c(as.character(min_detected_fraction),
                                                                       as.character(n_variable_proteins),
                                                                       as.character(n_scot_pcs), "TRUE",
                                                                       "log2_then_proda_median_normalization", "retained",
                                                                       "zero_after_feature_zscore"))
write_csv(parameter_report, file.path(qc_dir, "protein_preprocessing_parameters.csv"))
write_lines(capture.output(sessionInfo()), file.path(qc_dir, "protein_preprocessing_sessionInfo.txt"))
cat("\nProteomics preprocessing complete.\n\n")
cat("Primary MOFA+ outputs:\n")
cat("  ", file.path(mofa_dir, "mofa_protein_all68_z_with_na.csv"), "\n")
cat("  ", file.path(mofa_dir, "mofa_protein_paired58_z_with_na.csv"), "\n\n")
cat("Primary SCOT+ PCA outputs:\n")
cat("  ", file.path(scot_dir, "scot_protein_all68_pca_embedding.csv"), "\n")
cat("  ", file.path(scot_dir, "scot_protein_paired58_pca_embedding.csv"), "\n\n")
cat("QC summary:\n")
cat("  ", file.path(qc_dir, "protein_preprocessing_summary.csv"), "\n")
