# =========================================================
# Phosphoproteomics preprocessing for NanoSPLITS C10 RO-3306 experiment
# =========================================================
# Purpose:
#   Starting from:
#     C10Treated_Peptide_intensities.tsv
#
#   Generate:
#     1. Phosphopeptide sample and feature QC
#     2. MOFA+-ready phosphopeptide matrices
#        - all 68 author-QC proteomics cells
#        - strict paired RNA/protein-QC/peptide cells (expected n = 58)
#     3. Mechanistic-validation inputs using the all-68 cohort
#     4. Feature metadata linking every matrix row to its modified sequence, assigned modification, parent protein, and gene
#
# Modelling decisions:
#   - Mouse peptide rows are identified using "_MOUSE" in Entry.Name.
#   - Phosphopeptides are identified by the phosphorylation mass shift in Modified.Sequence: grepl("79\\.9", Modified.Sequence).
#   - Zero intensity means not detected and is converted to NA.
#   - Positive intensities are log2 transformed.
#   - proDA-style normalization offsets are estimated from all mouse peptides, not only the sparse phosphopeptide subset.
#   - Features detected in at least 20% of a cohort and with nonzero observed variance are retained for that cohort's MOFA+ view.
#   - The >=5 detected cells per condition rule is recorded separately for mechanistic testing;
#     it is not part of the MOFA filter because doing so could remove treatment-specific phosphorylation features.
#   - All phosphopeptides passing the MOFA filter are retained; no additional top-variable-feature selection
#     is applied to this small view.
#   - MOFA+ matrices retain NA values and are feature-wise z-scored using observed values only.
#   - No phosphoproteomic SCOT+ input is generated. Primary SCOT+ alignment is RNA <-> total protein;
#     phosphoproteomics is a MOFA+ view and mechanistic validation layer.
#   - Treatment is not regressed out.
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

  warning("Could not determine the script directory from RStudio or Rscript. ",
          "Using the current working directory.")

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()

raw_data_dir <- "C:/Users/49152/Downloads/2nd_paper/Github_data"

processed_dir <- file.path(script_dir, "processed_data")
shared_dir <- file.path(processed_dir, "shared")
mofa_dir <- file.path(processed_dir, "mofa_input")
mechanistic_dir <- file.path(processed_dir, "mechanistic_input")
qc_dir <- file.path(processed_dir, "qc")

walk(c(processed_dir, shared_dir, mofa_dir, mechanistic_dir, qc_dir), dir.create, showWarnings = FALSE, recursive = TRUE)

peptide_file <- file.path(raw_data_dir, "C10Treated_Peptide_intensities.tsv")
metadata_file <- file.path(processed_dir, "metadata_master_samples.csv")
phosphopeptide_feature_metadata_file <- file.path(processed_dir, "metadata_phosphopeptide_features.csv")

for (input_file in c(peptide_file, metadata_file, phosphopeptide_feature_metadata_file)) {
  if (!file.exists(input_file)) {
    stop("Required input file does not exist:\n", input_file,
         "\nRun the metadata preparation script first and check the input paths.")
  }
}

cat("Script directory:   ", script_dir, "\n")
cat("Raw data directory: ", raw_data_dir, "\n")
cat("Output directory:   ", processed_dir, "\n\n")

# =========================================================
# 2. User-adjustable parameters
# =========================================================

# Overall detection requirement for MOFA+ feature inclusion.
min_detected_fraction <- 0.20

# Minimum observed cells in each condition for downstream mechanistic tests.
# This annotates eligibility but is not part of the MOFA+ feature filter.
min_detected_per_condition <- 5L

# Experiment-specific validation values.
expected_n_peptide_samples <- 96L
expected_n_protein_qc <- 68L
expected_n_paired <- 58L
expected_n_phosphopeptide_rows <- 286L

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
  if (anyNA(ids)) {
    stop(object_name, " contains missing identifiers.")
  }

  if (any(!nzchar(trimws(as.character(ids))))) {
    stop(object_name, " contains blank identifiers.")
  }

  duplicated_ids <- unique(ids[duplicated(ids)])

  if (length(duplicated_ids) > 0) {
    stop(object_name, " contains duplicated identifiers. Examples: ", paste(head(duplicated_ids, 20), collapse = ", "))
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
  row_average <- rowMeans(mat, na.rm = TRUE)

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

calculate_phosphopeptide_qc <- function(mat, metadata, cohort_name, min_detected_fraction, min_detected_per_condition) {

  if (!identical(colnames(mat), metadata$SampleID)) {
    stop("Phosphopeptide QC matrix and metadata sample orders do not match.")
  }

  control_samples <- metadata$SampleID[metadata$Condition == "Control"]
  treated_samples <- metadata$SampleID[metadata$Condition == "Treated"]

  if (length(control_samples) == 0 || length(treated_samples) == 0) {
    stop("Both Control and Treated samples are required for phosphopeptide QC.")
  }

  detected_n <- rowSums(!is.na(mat))
  detected_fraction <- rowMeans(!is.na(mat))
  control_detected_n <- rowSums(!is.na(mat[, control_samples, drop = FALSE]))
  treated_detected_n <- rowSums(!is.na(mat[, treated_samples, drop = FALSE]))
  observed_variance <- apply(mat, 1, var, na.rm = TRUE)

  tibble(PeptideRowID = rownames(mat), Cohort = cohort_name, DetectedN = detected_n, MissingN = ncol(mat) - detected_n,
         DetectedFraction = detected_fraction, MissingFraction = 1 - detected_fraction, ControlN = length(control_samples),
         ControlDetectedN = control_detected_n, ControlDetectedFraction = control_detected_n / length(control_samples),
         TreatedN = length(treated_samples), TreatedDetectedN = treated_detected_n,
         TreatedDetectedFraction = treated_detected_n / length(treated_samples),
         DetectionFractionDifference_TreatedMinusControl = TreatedDetectedFraction - ControlDetectedFraction,
         MeanObserved = rowMeans(mat, na.rm = TRUE), MedianObserved = apply(mat, 1, median, na.rm = TRUE),
         VarianceObserved = observed_variance, PassDetectionFilter = DetectedFraction >= min_detected_fraction,
         PassVarianceFilter = is.finite(VarianceObserved) & VarianceObserved > 0,
         PassMOFAFilter = PassDetectionFilter & PassVarianceFilter, 
         MechanisticTestEligible = ControlDetectedN >= min_detected_per_condition &
           TreatedDetectedN >= min_detected_per_condition & PassVarianceFilter)
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

  if (nrow(mat_keep) == 0) {
    stop("No phosphopeptides remain after scaling filters.")
  }

  z <- sweep(mat_keep, 1, means_keep, FUN = "-")
  z <- sweep(z, 1, sds_keep, FUN = "/")

  if (any(!is.finite(z[!is.na(z)]))) {
    stop("Non-finite observed values remain after phosphopeptide z-scoring.")
  }

  list(matrix = z, means = means_keep, sds = sds_keep, dropped_features = dropped_features)
}

assert_mofa_export <- function(mat, metadata) {
  if (!identical(colnames(mat), metadata$SampleID)) {
    stop("MOFA+ phosphopeptide matrix columns do not match metadata order.")
  }

  if (any(!is.finite(mat[!is.na(mat)]))) {
    stop("MOFA+ phosphopeptide matrix contains non-finite observed values.")
  }

  invisible(TRUE)
}

make_observed_long_table <- function(log2_mat, z_mat, feature_metadata, sample_metadata) {
  if (!identical(rownames(log2_mat), rownames(z_mat)) || !identical(colnames(log2_mat), colnames(z_mat))) {
    stop("Mechanistic log2 and z-score matrices do not have matching dimensions.")
  }

  log2_long <- as.data.frame(log2_mat, check.names = FALSE) %>%
    rownames_to_column(var = "PeptideRowID") %>%
    pivot_longer(cols = -PeptideRowID, names_to = "SampleID", values_to = "Log2MedianNormalizedIntensity")

  z_long <- as.data.frame(z_mat, check.names = FALSE) %>%
    rownames_to_column(var = "PeptideRowID") %>%
    pivot_longer(cols = -PeptideRowID, names_to = "SampleID", values_to = "FeatureZScore")

  log2_long %>%
    left_join(z_long, by = c("PeptideRowID", "SampleID")) %>%
    filter(!is.na(Log2MedianNormalizedIntensity)) %>%
    left_join(feature_metadata, by = "PeptideRowID") %>%
    left_join(sample_metadata, by = "SampleID")
}

# =========================================================
# 4. Load and validate master metadata
# =========================================================

cat("Loading master metadata...\n")

master_metadata <- read_csv(metadata_file, show_col_types = FALSE)

required_meta_cols <- c("SampleID", "TreatmentCode", "Condition", "ConditionLabel", "Block", "WellIndex", "Peptide_available",
                        "Phosphopeptide_available", "Peptide_detected_n", "Phosphopeptide_detected_n", "ProteinQC68",
                        "Phosphopeptide_on_ProteinQC68", "Paired_RNA_ProteinQC_Peptide")

missing_meta_cols <- setdiff(required_meta_cols, colnames(master_metadata))

if (length(missing_meta_cols) > 0) {
  stop("metadata_master_samples.csv is missing required column(s): ", paste(missing_meta_cols, collapse = ", "))
}

validate_unique_ids(master_metadata$SampleID, "Master metadata SampleID column")

logical_meta_cols <- c("Peptide_available", "Phosphopeptide_available", "ProteinQC68", "Phosphopeptide_on_ProteinQC68",
                       "Paired_RNA_ProteinQC_Peptide")

for (column_name in logical_meta_cols) {
  if (!is.logical(master_metadata[[column_name]])) {
    stop(column_name, " must be imported as a logical TRUE/FALSE column.")
  }
}

peptide_sample_ids <- master_metadata %>%
  filter(Peptide_available) %>%
  pull(SampleID)

phosphopeptide_all68_ids <- master_metadata %>%
  filter(Phosphopeptide_on_ProteinQC68) %>%
  pull(SampleID)

phosphopeptide_paired58_ids <- master_metadata %>%
  filter(Paired_RNA_ProteinQC_Peptide) %>%
  pull(SampleID)

validate_unique_ids(peptide_sample_ids, "Peptide sample-ID set")
validate_unique_ids(phosphopeptide_all68_ids, "Phosphopeptide all-68 sample-ID set")
validate_unique_ids(phosphopeptide_paired58_ids, "Phosphopeptide paired-58 sample-ID set")

if (!all(phosphopeptide_all68_ids %in% peptide_sample_ids)) {
  stop("The phosphopeptide all-68 set is not a subset of peptide-available samples.")
}

if (!all(phosphopeptide_paired58_ids %in% phosphopeptide_all68_ids)) {
  stop("The paired-58 phosphopeptide set is not a subset of the all-68 set.")
}

observed_sample_counts <- c(PeptideSamples = length(peptide_sample_ids), ProteinQC68 = length(phosphopeptide_all68_ids),
                            Paired58 = length(phosphopeptide_paired58_ids))

expected_sample_counts <- c(PeptideSamples = expected_n_peptide_samples, ProteinQC68 = expected_n_protein_qc,
                            Paired58 = expected_n_paired)

if (any(observed_sample_counts != expected_sample_counts)) {
  stop("Unexpected phosphoproteomics sample-set sizes. Observed: ",
       paste(names(observed_sample_counts), observed_sample_counts, sep = "=", collapse = ", "), ". Expected: ",
       paste(names(expected_sample_counts), expected_sample_counts, sep = "=", collapse = ", "))
}

metadata_all68 <- order_metadata_to_samples(master_metadata, phosphopeptide_all68_ids)
metadata_paired58 <- order_metadata_to_samples(master_metadata, phosphopeptide_paired58_ids)

if (anyNA(metadata_all68$Condition) || anyNA(metadata_paired58$Condition)) {
  stop("One or more phosphoproteomics samples lack condition metadata.")
}

cat("Peptide-available samples:             ", length(peptide_sample_ids), "\n")
cat("Phosphopeptide all-68 samples:        ", length(phosphopeptide_all68_ids), "\n")
cat("Strict paired phosphopeptide samples:", length(phosphopeptide_paired58_ids), "\n\n")

# =========================================================
# 5. Load and validate peptide matrix
# =========================================================

cat("Loading peptide intensity matrix...\n")

peptide_raw <- read_tsv(peptide_file, show_col_types = FALSE, progress = FALSE, name_repair = "minimal")

if (anyDuplicated(colnames(peptide_raw))) {
  duplicated_columns <- unique(colnames(peptide_raw)[duplicated(colnames(peptide_raw))])

  stop("Duplicated column names detected in peptide input: ", paste(duplicated_columns, collapse = ", "))
}

required_peptide_annot <- c("Peptide.Sequence", "Modified.Sequence", "Assigned.Modifications", "Protein", "ID", "Entry.Name",
                            "Gene", "Protein.Description")

missing_peptide_annot <- setdiff(required_peptide_annot, colnames(peptide_raw))

if (length(missing_peptide_annot) > 0) {
  stop("Peptide file is missing expected annotation column(s): ", paste(missing_peptide_annot, collapse = ", "))
}

input_sample_cols <- setdiff(colnames(peptide_raw), required_peptide_annot)

missing_sample_cols <- setdiff(peptide_sample_ids, input_sample_cols)
unexpected_sample_cols <- setdiff(input_sample_cols, peptide_sample_ids)

if (length(missing_sample_cols) > 0 || length(unexpected_sample_cols) > 0) {
  stop("Peptide sample columns do not match metadata.\nMissing from peptide file: ",
       ifelse(length(missing_sample_cols) == 0, "none", paste(missing_sample_cols, collapse = ", ")),
       "\nUnexpected non-annotation columns: ", ifelse(length(unexpected_sample_cols) == 0, "none",
                                                       paste(unexpected_sample_cols, collapse = ", ")))
}

# Add the persistent row identifier before filtering, matching the metadata script.
peptide_mouse <- peptide_raw %>%
  mutate(PeptideRowID = sprintf("PEP_%06d", row_number())) %>%
  filter(!is.na(Entry.Name), grepl("_MOUSE", Entry.Name, fixed = TRUE))

if (nrow(peptide_mouse) == 0) {
  stop("No mouse peptide rows were detected using Entry.Name.")
}

validate_unique_ids(peptide_mouse$PeptideRowID, "Mouse peptide row IDs")

phosphopeptide_mouse <- peptide_mouse %>%
  filter(!is.na(Modified.Sequence), grepl("79\\.9", Modified.Sequence))

if (nrow(phosphopeptide_mouse) == 0) {
  stop("No phosphopeptides were detected using the 79.9 mass-shift pattern.")
}

if (nrow(phosphopeptide_mouse) != expected_n_phosphopeptide_rows) {
  stop("Expected ", expected_n_phosphopeptide_rows, " mouse phosphopeptide rows but detected ", nrow(phosphopeptide_mouse),
       ". Check the input file and phosphorylation pattern.")
}

validate_unique_ids(phosphopeptide_mouse$Modified.Sequence, "Mouse phosphopeptide Modified.Sequence identifiers")

phosphopeptide_feature_metadata <- phosphopeptide_mouse %>%
  select(PeptideRowID, all_of(required_peptide_annot)) %>%
  mutate(IsPhosphopeptide = TRUE)

saved_phosphopeptide_feature_metadata <- read_csv(phosphopeptide_feature_metadata_file, show_col_types = FALSE)

if (!"PeptideRowID" %in% colnames(saved_phosphopeptide_feature_metadata)) {
  stop("Saved phosphopeptide feature metadata lacks PeptideRowID.")
}

if (!identical(phosphopeptide_feature_metadata$PeptideRowID, saved_phosphopeptide_feature_metadata$PeptideRowID)) {
  stop("Phosphopeptide row IDs do not match the metadata preparation output.")
}

# Use all mouse peptides to estimate stable normalization offsets.
peptide_mat_raw <- peptide_mouse %>%
  select(all_of(peptide_sample_ids)) %>%
  safe_numeric_matrix_allow_na(object_name = "Raw mouse peptide intensity matrix")

rownames(peptide_mat_raw) <- peptide_mouse$PeptideRowID

if (any(peptide_mat_raw < 0, na.rm = TRUE)) {
  stop("Raw peptide intensity matrix contains negative values.")
}

peptide_mat_raw[peptide_mat_raw == 0] <- NA_real_

phosphopeptide_row_ids <- phosphopeptide_mouse$PeptideRowID
phosphopeptide_mat_raw <- peptide_mat_raw[phosphopeptide_row_ids, peptide_sample_ids, drop = FALSE]

# =========================================================
# 6. Validate peptide detection counts against metadata
# =========================================================

peptide_detected_recalculated <- colSums(!is.na(peptide_mat_raw))
phosphopeptide_detected_recalculated <- colSums(!is.na(phosphopeptide_mat_raw))

detection_validation <- master_metadata %>%
  filter(Peptide_available) %>%
  select(SampleID, Peptide_detected_n_metadata = Peptide_detected_n,
         Phosphopeptide_detected_n_metadata = Phosphopeptide_detected_n) %>%
  mutate(Peptide_detected_n_recalculated = unname(peptide_detected_recalculated[SampleID]),
         Phosphopeptide_detected_n_recalculated = unname(phosphopeptide_detected_recalculated[SampleID]))

if (!isTRUE(all.equal(unname(detection_validation$Peptide_detected_n_recalculated),
                      unname(detection_validation$Peptide_detected_n_metadata), tolerance = 0, check.attributes = FALSE))) {
  stop("Recalculated mouse peptide detection counts do not match metadata.")
}

if (!isTRUE(all.equal(unname(detection_validation$Phosphopeptide_detected_n_recalculated),
                      unname(detection_validation$Phosphopeptide_detected_n_metadata), tolerance = 0, check.attributes = FALSE))) {
  stop("Recalculated phosphopeptide detection counts do not match metadata.")
}

write_csv(detection_validation, file.path(qc_dir, "phosphopeptide_sample_detection_validation.csv"))

# =========================================================
# 7. Log2 transformation and proDA-style median normalization
# =========================================================

peptide_raw_all68 <- peptide_mat_raw[ , phosphopeptide_all68_ids, drop = FALSE]

peptide_raw_paired58 <- peptide_mat_raw[ , phosphopeptide_paired58_ids, drop = FALSE]

peptide_log2_all68 <- log2(peptide_raw_all68)
peptide_log2_paired58 <- log2(peptide_raw_paired58)

normalized_all68 <- proda_median_normalize(peptide_log2_all68)
normalized_paired58 <- proda_median_normalize(peptide_log2_paired58)

peptide_norm_all68 <- normalized_all68$matrix
peptide_norm_paired58 <- normalized_paired58$matrix

phosphopeptide_norm_all68 <- peptide_norm_all68[phosphopeptide_row_ids, phosphopeptide_all68_ids, drop = FALSE]

phosphopeptide_norm_paired58 <- peptide_norm_paired58[phosphopeptide_row_ids, phosphopeptide_paired58_ids, drop = FALSE]

write_matrix(phosphopeptide_norm_all68,
             file.path(shared_dir, "phosphopeptide_log2_proda_median_normalized_all68_all_features.csv"), "PeptideRowID")

write_matrix(phosphopeptide_norm_paired58, 
             file.path(shared_dir, "phosphopeptide_log2_proda_median_normalized_paired58_all_features.csv"), "PeptideRowID")

normalization_offsets <- bind_rows(tibble(Cohort = "all68", SampleID = names(normalized_all68$sample_offsets),
                                          NormalizationOffset = unname(normalized_all68$sample_offsets)),
                                   tibble(Cohort = "paired58", SampleID = names(normalized_paired58$sample_offsets),
                                          NormalizationOffset = unname(normalized_paired58$sample_offsets)))

write_csv(normalization_offsets, file.path(qc_dir, "phosphoproteomics_normalization_offsets.csv"))

# =========================================================
# 8. Sample QC
# =========================================================

sample_qc_all68 <- tibble(Cohort = "all68", SampleID = phosphopeptide_all68_ids,
                          AllMousePeptidesDetected = colSums(!is.na(peptide_raw_all68)),
                          PhosphopeptidesDetected = colSums(!is.na(phosphopeptide_norm_all68)),
                          PhosphopeptidesMissing = colSums(is.na(phosphopeptide_norm_all68)),
                          PhosphopeptideMissingFraction = colMeans(is.na(phosphopeptide_norm_all68)),
                          NormalizationOffset = unname(normalized_all68$sample_offsets)) %>%
  left_join(metadata_all68 %>% select(SampleID, TreatmentCode, Condition, ConditionLabel, Block, WellIndex), by = "SampleID")

sample_qc_paired58 <- tibble(Cohort = "paired58", SampleID = phosphopeptide_paired58_ids,
                             AllMousePeptidesDetected = colSums(!is.na(peptide_raw_paired58)),
                             PhosphopeptidesDetected = colSums(!is.na(phosphopeptide_norm_paired58)),
                             PhosphopeptidesMissing = colSums(is.na(phosphopeptide_norm_paired58)),
                             PhosphopeptideMissingFraction = colMeans(is.na(phosphopeptide_norm_paired58)),
                             NormalizationOffset = unname(normalized_paired58$sample_offsets)) %>%
  left_join(metadata_paired58 %>% select(SampleID, TreatmentCode, Condition, ConditionLabel, Block, WellIndex), by = "SampleID")

sample_qc_combined <- bind_rows(sample_qc_all68, sample_qc_paired58)

if (anyNA(sample_qc_combined$Condition)) {
  stop("One or more phosphoproteomics QC rows lack condition metadata.")
}

write_csv(sample_qc_combined, file.path(qc_dir, "phosphopeptide_sample_qc_metrics_by_cohort.csv"))

cat("Phosphopeptide sample QC summary:\n")
print(sample_qc_combined %>% group_by(Cohort, Condition) %>%
        summarise(n = n(), median_detected = median(PhosphopeptidesDetected),
                  median_missing_fraction = median(PhosphopeptideMissingFraction), .groups = "drop"))

# =========================================================
# 9. Cohort-specific feature QC and MOFA filtering
# =========================================================

phosphopeptide_qc_all68 <- calculate_phosphopeptide_qc(phosphopeptide_norm_all68, metadata_all68, "all68",
                                                       min_detected_fraction, min_detected_per_condition)

phosphopeptide_qc_paired58 <- calculate_phosphopeptide_qc(phosphopeptide_norm_paired58, metadata_paired58, "paired58",
                                                          min_detected_fraction, min_detected_per_condition)

phosphopeptide_qc_combined <- bind_rows(phosphopeptide_qc_all68, phosphopeptide_qc_paired58) %>%
  left_join(phosphopeptide_feature_metadata, by = "PeptideRowID")

write_csv(phosphopeptide_qc_combined, file.path(qc_dir, "phosphopeptide_feature_qc_metrics_by_cohort.csv"))

features_keep_all68 <- phosphopeptide_qc_all68 %>% filter(PassMOFAFilter) %>% pull(PeptideRowID)

features_keep_paired58 <- phosphopeptide_qc_paired58 %>% filter(PassMOFAFilter) %>% pull(PeptideRowID)

if (length(features_keep_all68) == 0) {
  stop("No phosphopeptides passed the all-68 MOFA filter.")
}

if (length(features_keep_paired58) == 0) {
  stop("No phosphopeptides passed the paired-58 MOFA filter.")
}

write_lines(features_keep_all68, file.path(shared_dir, "phosphopeptide_features_after_filter_all68.txt"))

write_lines(features_keep_paired58, file.path(shared_dir, "phosphopeptide_features_after_filter_paired58.txt"))

phosphopeptide_filtered_all68 <- phosphopeptide_norm_all68[features_keep_all68, phosphopeptide_all68_ids, drop = FALSE]

phosphopeptide_filtered_paired58 <- phosphopeptide_norm_paired58[features_keep_paired58, phosphopeptide_paired58_ids,
                                                                 drop = FALSE]

# =========================================================
# 10. Feature-wise scaling with missing values retained
# =========================================================

scaled_all68 <- zscore_rows_preserve_na(phosphopeptide_filtered_all68)
scaled_paired58 <- zscore_rows_preserve_na(phosphopeptide_filtered_paired58)

phosphopeptide_all68_z <- scaled_all68$matrix
phosphopeptide_paired58_z <- scaled_paired58$matrix

write_lines(scaled_all68$dropped_features, file.path(qc_dir, "phosphopeptide_all68_features_dropped_during_scaling.txt"))

write_lines(scaled_paired58$dropped_features, file.path(qc_dir, "phosphopeptide_paired58_features_dropped_during_scaling.txt"))

write_csv(tibble(PeptideRowID = names(scaled_all68$means), ObservedMean = unname(scaled_all68$means),
                 ObservedSD = unname(scaled_all68$sds)), file.path(qc_dir, "phosphopeptide_scaling_parameters_all68.csv"))

write_csv(tibble(PeptideRowID = names(scaled_paired58$means), ObservedMean = unname(scaled_paired58$means),
                 ObservedSD = unname(scaled_paired58$sds)), file.path(qc_dir, "phosphopeptide_scaling_parameters_paired58.csv"))

assert_mofa_export(phosphopeptide_all68_z, metadata_all68)
assert_mofa_export(phosphopeptide_paired58_z, metadata_paired58)

# =========================================================
# 11. MOFA+ exports
# =========================================================

write_matrix(phosphopeptide_all68_z, file.path(mofa_dir, "mofa_phosphopeptide_all68_z_with_na.csv"), "PeptideRowID")

write_matrix(phosphopeptide_paired58_z, file.path(mofa_dir, "mofa_phosphopeptide_paired58_z_with_na.csv"), "PeptideRowID")

write_matrix(phosphopeptide_filtered_all68, file.path(mofa_dir, "mofa_phosphopeptide_all68_log2_proda_median_normalized_with_na.csv"),
             "PeptideRowID")

write_matrix(phosphopeptide_filtered_paired58, file.path(mofa_dir, "mofa_phosphopeptide_paired58_log2_proda_median_normalized_with_na.csv"),
             "PeptideRowID")

write_csv(metadata_all68, file.path(mofa_dir, "mofa_phosphopeptide_all68_metadata.csv"))

write_csv(metadata_paired58, file.path(mofa_dir, "mofa_phosphopeptide_paired58_metadata.csv"))

# =========================================================
# 12. Mechanistic-validation inputs
# =========================================================

# Use all 68 QC cells for phosphopeptide-only mechanistic validation because
# RNA pairing is not required for this analysis.
mechanistic_feature_ids_all68 <- phosphopeptide_qc_all68 %>%
  filter(MechanisticTestEligible) %>%
  pull(PeptideRowID)

if (length(mechanistic_feature_ids_all68) == 0) {
  stop("No phosphopeptides met the mechanistic-test eligibility rule.")
}

mechanistic_log2_all68 <- phosphopeptide_norm_all68[mechanistic_feature_ids_all68, phosphopeptide_all68_ids, drop = FALSE]

mechanistic_scaled_all68 <- zscore_rows_preserve_na(mechanistic_log2_all68)
mechanistic_z_all68 <- mechanistic_scaled_all68$matrix

write_lines(mechanistic_feature_ids_all68, file.path(mechanistic_dir, "phosphopeptide_features_min5_per_condition_all68.txt"))

write_matrix(mechanistic_log2_all68, file.path(mechanistic_dir, "phosphopeptide_mechanistic_all68_log2_proda_median_normalized_with_na.csv"),
             "PeptideRowID")

write_matrix(mechanistic_z_all68, file.path(mechanistic_dir, "phosphopeptide_mechanistic_all68_z_with_na.csv"), "PeptideRowID")

mechanistic_long_all68 <- make_observed_long_table(mechanistic_log2_all68, mechanistic_z_all68, phosphopeptide_feature_metadata,
                                                   metadata_all68 %>%
                                                     select(SampleID, TreatmentCode, Condition, ConditionLabel, Block, WellIndex))

write_csv(mechanistic_long_all68, file.path(mechanistic_dir, "phosphopeptide_mechanistic_all68_observed_long.csv"))

# This script intentionally does not adjust phosphopeptide abundance for parent total-protein abundance.
# That comparison should be performed downstream using matched parent-protein measurements,
# for example by modelling phosphopeptide abundance with condition and parent-protein abundance as separate predictors.

# =========================================================
# 13. Preprocessed feature metadata
# =========================================================

preprocessed_feature_metadata <- phosphopeptide_feature_metadata %>%
  mutate(SelectedForMOFA_all68 = PeptideRowID %in% rownames(phosphopeptide_all68_z),
         SelectedForMOFA_paired58 = PeptideRowID %in% rownames(phosphopeptide_paired58_z),
         MechanisticTestEligible_all68 = PeptideRowID %in% mechanistic_feature_ids_all68) %>%
  left_join(phosphopeptide_qc_all68 %>% select(PeptideRowID, DetectedFraction_all68 = DetectedFraction,
                                               ControlDetectedN_all68 = ControlDetectedN,
                                               TreatedDetectedN_all68 = TreatedDetectedN,
                                               DetectionFractionDifference_all68 = DetectionFractionDifference_TreatedMinusControl,
                                               VarianceObserved_all68 = VarianceObserved), by = "PeptideRowID") %>%
  left_join(phosphopeptide_qc_paired58 %>% select(PeptideRowID, DetectedFraction_paired58 = DetectedFraction,
                                                  ControlDetectedN_paired58 = ControlDetectedN,
                                                  TreatedDetectedN_paired58 = TreatedDetectedN,
                                                  DetectionFractionDifference_paired58 = DetectionFractionDifference_TreatedMinusControl,
                                                  VarianceObserved_paired58 = VarianceObserved), by = "PeptideRowID")

write_csv(preprocessed_feature_metadata, file.path(shared_dir, "phosphopeptide_feature_metadata_preprocessed.csv"))

# =========================================================
# 14. Final checks and summary
# =========================================================

stopifnot(identical(colnames(phosphopeptide_all68_z), metadata_all68$SampleID),
          identical(colnames(phosphopeptide_paired58_z), metadata_paired58$SampleID),
          nrow(phosphopeptide_all68_z) > 0, nrow(phosphopeptide_paired58_z) > 0,
          all(is.finite(phosphopeptide_all68_z[!is.na(phosphopeptide_all68_z)])),
          all(is.finite(phosphopeptide_paired58_z[!is.na(phosphopeptide_paired58_z)])))

summary_report <- tibble(Item = c("Input peptide samples", "Input mouse peptide rows", "Extracted mouse phosphopeptide rows",
                                  "All68 samples", "Paired58 samples", "Minimum detected fraction for MOFA",
                                  "Minimum detected per condition for mechanistic tests", "MOFA features all68", 
                                  "MOFA features paired58", "Mechanistic-test-eligible features all68",
                                  "MOFA missing fraction all68", "MOFA missing fraction paired58", 
                                  "SCOT phosphopeptide export generated"),
                         Value = as.character(c(ncol(peptide_mat_raw), nrow(peptide_mat_raw), nrow(phosphopeptide_mat_raw),
                                                length(phosphopeptide_all68_ids), length(phosphopeptide_paired58_ids),
                                                min_detected_fraction, min_detected_per_condition, nrow(phosphopeptide_all68_z),
                                                nrow(phosphopeptide_paired58_z), length(mechanistic_feature_ids_all68),
                                                mean(is.na(phosphopeptide_all68_z)), mean(is.na(phosphopeptide_paired58_z)),
                                                FALSE)))

write_csv(summary_report, file.path(qc_dir, "phosphoproteomics_preprocessing_summary.csv"))

parameter_report <- tibble(Parameter = c("phosphorylation_pattern", "min_detected_fraction", "min_detected_per_condition",
                                         "normalization_reference", "normalization", "feature_selection", "mofa_missing_values",
                                         "scot_export"),
                           Value = c("79\\.9_in_Modified.Sequence", as.character(min_detected_fraction),
                                     as.character(min_detected_per_condition), "all_mouse_peptides_within_each_cohort",
                                     "log2_then_proda_median_normalization", "all_features_passing_detection_and_variance_filters",
                                     "retained", "not_generated"))

write_csv(parameter_report, file.path(qc_dir, "phosphoproteomics_preprocessing_parameters.csv"))

write_lines(capture.output(sessionInfo()), file.path(qc_dir, "phosphoproteomics_preprocessing_sessionInfo.txt"))

# =========================================================
# 15. Completion message
# =========================================================

cat("\nPhosphoproteomics preprocessing complete.\n\n")

cat("Primary MOFA+ outputs:\n")
cat("  ", file.path(mofa_dir, "mofa_phosphopeptide_all68_z_with_na.csv"), "\n")
cat("  ", file.path(mofa_dir, "mofa_phosphopeptide_paired58_z_with_na.csv"), "\n\n")

cat("Primary mechanistic-validation outputs:\n")
cat("  ", file.path(mechanistic_dir, "phosphopeptide_mechanistic_all68_log2_proda_median_normalized_with_na.csv"), "\n")
cat("  ", file.path(mechanistic_dir, "phosphopeptide_mechanistic_all68_observed_long.csv"), "\n\n")

cat("QC summary:\n")
cat("  ", file.path(qc_dir, "phosphoproteomics_preprocessing_summary.csv"), "\n")
