# =========================================================
# Transcriptomics preprocessing for NanoSPLITS C10 RO-3306 experiment
# =========================================================
# Purpose:
#   Starting from:
#     C10Treated_NormCounts.csv
#
#   Generate:
#     1. Shared filtered and logCPM-transformed RNA matrix
#     2. MOFA+-ready RNA matrices
#        - all 62 RNA cells
#        - strict paired RNA/protein-QC benchmark cells
#     3. SCOT+-ready RNA matrices
#        - all 62 RNA cells
#        - strict paired RNA/protein-QC benchmark cells
#     4. SCOT+-ready RNA PCA embeddings
#     5. RNA QC, feature-selection, and preprocessing metadata
#
# Important modelling decisions:
#   - C10Treated_NormCounts.csv is treated as a non-negative, continuous count-like matrix.
#   - edgeR logCPM is used to remove remaining sample library-size effects
#     and produce approximately Gaussian continuous RNA input.
#   - The matrix is transformed to logCPM using edgeR::cpm(..., log = TRUE, prior.count = 1).
#   - Treatment is not regressed out because CDK1 inhibition / G2M arrest
#     is the biological signal of interest.
#   - MOFA+ outputs are features x samples.
#   - SCOT+ canonical outputs are samples x features.
#   - SCOT+ matrices must be complete and finite.
#   - Paired benchmark feature selection and scaling are performed using
#     only the paired benchmark cells.
#
# =========================================================
# Recommended package versions
# =========================================================
# R >= 4.4
# tidyverse >= 2.0.0
# readr >= 2.1.5
# rstudioapi >= 0.17.1
# edgeR >= 4.0

library(tidyverse)
library(rstudioapi)
library(edgeR)

# =========================================================
# 1. Path setup
# =========================================================

get_script_dir <- function() {
  # Case 1: script is run from RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    
    if (!is.null(ctx) && !is.null(ctx$path) && nzchar(ctx$path)) {
      return(normalizePath(dirname(ctx$path), winslash = "/", mustWork = TRUE))
    }
  }
  
  # Case 2: script is run using Rscript
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    
    return(normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE))
  }
  
  # Case 3: fallback to current working directory
  warning("Could not determine the script directory from RStudio or Rscript. ",
          "Using the current working directory.")
  
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()

# Raw-data directory supplied by user
raw_data_dir <- "C:/Users/49152/Downloads/2nd_paper/Github_data"

# Main processed-data directory
processed_dir <- file.path(script_dir, "processed_data")

# Subdirectories
shared_dir <- file.path(processed_dir, "shared")

mofa_dir <- file.path(processed_dir, "mofa_input")

scot_dir <- file.path(processed_dir, "scot_input")

qc_dir <- file.path(processed_dir, "qc")

walk(c(processed_dir, shared_dir, mofa_dir, scot_dir, qc_dir), dir.create, showWarnings = FALSE,
     recursive = TRUE)

# Input files
rna_file <- file.path(raw_data_dir, "C10Treated_NormCounts.csv")

metadata_file <- file.path(processed_dir, "metadata_master_samples.csv")

# Informative input checks
if (!file.exists(rna_file)) {
  stop("RNA input file does not exist:\n", rna_file)
}

if (!file.exists(metadata_file)) {
  stop("Master metadata file does not exist:\n", metadata_file,
       "\nRun the metadata preparation script first.")
}

cat("Script directory:   ", script_dir, "\n")
cat("Raw data directory: ", raw_data_dir, "\n")
cat("Output directory:   ", processed_dir, "\n\n")

# =========================================================
# 2. User-adjustable parameters
# =========================================================

# Retain genes detected above zero in at least this fraction of cells in the relevant analysis cohort.
min_detected_fraction <- 0.20

# Number of top variable genes used for both MOFA+ and SCOT+
n_variable_genes <- 2000

# Number of PCs exported for SCOT+
n_scot_pcs <- 30

# Average prior count used for the edgeR logCPM transformation
logcpm_prior_count <- 1

# Expected sample-set sizes.
# Deviations produce hard failures because output filenames and downstream
# analyses explicitly assume 62 RNA cells and 58 paired benchmark cells.
expected_n_rna <- 62
expected_n_paired <- 58

# =========================================================
# 3. Helper functions
# =========================================================

safe_numeric_matrix <- function(df, object_name = "matrix") {
  mat <- as.matrix(df)
  
  suppressWarnings(storage.mode(mat) <- "numeric")
  
  if (anyNA(mat)) {
    stop(object_name, " contains missing or non-numeric values.")
  }
  
  if (any(!is.finite(mat))) {
    stop(object_name, " contains non-finite values.")
  }
  
  mat
}

validate_unique_ids <- function(ids, object_name) {
  
  if (anyNA(ids)) {
    stop(object_name, " contains missing identifiers.")
  }
  
  if (any(!nzchar(trimws(as.character(ids))))) {
    stop(object_name, " contains blank or whitespace-only identifiers.")
  }
  
  duplicated_ids <- unique(ids[duplicated(ids)])
  
  if (length(duplicated_ids) > 0) {
    stop(object_name, " contains duplicated identifiers: ",
         paste(duplicated_ids, collapse = ", "))
  }
  
  invisible(TRUE)
}

order_metadata_to_samples <- function(
    metadata,
    sample_ids) {
  
  index <- match(sample_ids, metadata$SampleID)
  
  if (anyNA(index)) {
    stop("Metadata could not be matched for sample(s): ",
         paste(sample_ids[is.na(index)], collapse = ", "))
  }
  
  output <- metadata[index, ,drop = FALSE]
  
  if (!identical(output$SampleID, sample_ids)) {
    stop("Internal error while ordering metadata to sample IDs.")
  }
  
  output
}

write_matrix <- function(mat, file, id_column) {
  as.data.frame(mat, check.names = FALSE) %>%
    rownames_to_column(var = id_column) %>%
    write_csv(file)
}

calculate_gene_qc <- function(mat, cohort_name, min_detected_fraction) {
  
  qc <- tibble(Gene = rownames(mat), Cohort = cohort_name,
               DetectedN_gt0 = rowSums(mat > 0, na.rm = TRUE),
               DetectedFraction_gt0 = rowMeans(mat > 0, na.rm = TRUE),
               ZeroFraction = rowMeans(mat == 0, na.rm = TRUE),
               MeanExpression = rowMeans(mat, na.rm = TRUE),
               VarianceExpression = apply(mat, 1, var, na.rm = TRUE)) %>%
    mutate(PassDetectionFilter = DetectedFraction_gt0 >= min_detected_fraction,
           PassVarianceFilter = is.finite(VarianceExpression) & VarianceExpression > 0,
           PassCombinedFilter = PassDetectionFilter & PassVarianceFilter)
  qc
}

select_top_variable_genes <- function(mat, n_features) {
  
  feature_variance <- apply(mat, 1, var, na.rm = TRUE)
  
  valid <- is.finite(feature_variance) & feature_variance > 0
  
  feature_variance <- feature_variance[valid]
  
  if (length(feature_variance) == 0) {
    stop("No variable genes are available for feature selection.")
  }
  
  n_use <- min(n_features, length(feature_variance))
  
  if (n_use < n_features) {
    warning("Requested ", n_features, " variable genes, but only ", n_use,
            " variable genes were available.")
  }
  
  selected <- names(sort(feature_variance, decreasing = TRUE))[seq_len(n_use)]
  
  list(selected = selected, variance = feature_variance)
}

zscore_rows <- function(mat) {
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  
  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  
  keep <- is.finite(row_sd) & row_sd > 0
  
  dropped_features <- rownames(mat)[!keep]
  
  mat_keep <- mat[keep, , drop = FALSE]
  
  if (nrow(mat_keep) == 0) {
    stop("No features remain after removing zero-variance features.")
  }
  
  row_means <- rowMeans(mat_keep, na.rm = TRUE)
  
  row_sds <- apply(mat_keep, 1, sd, na.rm = TRUE)
  
  z <- sweep(mat_keep, 1, row_means, FUN = "-")
  
  z <- sweep(z, 1, row_sds, FUN = "/")
  
  if (anyNA(z)) {
    stop("Missing values remain after row-wise z-scoring.")
  }
  
  if (any(!is.finite(z))) {
    stop("Non-finite values remain after row-wise z-scoring.")
  }
  
  list(matrix = z, dropped_features = dropped_features)
}

safe_cor <- function(x, y, method = "pearson") {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) < 3) {
    return(NA_real_)
  }
  
  sd_x <- sd(x)
  sd_y <- sd(y)
  
  if (!is.finite(sd_x) || !is.finite(sd_y) || sd_x == 0 || sd_y == 0) {
    return(NA_real_)
  }
  
  cor(x, y, method = method)
}

make_scot_pca <- function(samples_by_features, n_pcs = 30, center = TRUE, scale_features = FALSE) {
  
  if (anyNA(samples_by_features) || any(!is.finite(samples_by_features))) {
    stop("SCOT+ PCA input contains missing or non-finite values.")
  }
  
  n_pcs_use <- min(n_pcs, nrow(samples_by_features) - 1, ncol(samples_by_features))
  
  if (n_pcs_use < 2) {
    stop("Too few samples or features to calculate a PCA embedding.")
  }
  
  pca <- prcomp(samples_by_features, center = center, scale. = scale_features, rank. = n_pcs_use)
  
  embedding <- pca$x[ , seq_len(n_pcs_use), drop = FALSE]
  
  # Calculate fractions relative to total PCA variance, but export
  # only the PCs present in the embedding.
  all_variance_fraction <- pca$sdev^2 / sum(pca$sdev^2)
  
  variance_fraction <- all_variance_fraction[seq_len(n_pcs_use)]
  
  variance_explained <- tibble(PC = colnames(embedding), VarianceExplained = variance_fraction,
                               CumulativeVarianceExplained = cumsum(variance_fraction))
  
  list(embedding = embedding, variance_explained = variance_explained, rotation = pca$rotation)
}

assert_mofa_export <- function(mat, metadata) {
  
  if (!identical(colnames(mat), metadata$SampleID)) {
    stop("MOFA+ matrix column order does not match metadata SampleID order.")
  }
  
  if (anyNA(mat)) {
    stop("This complete RNA MOFA+ export unexpectedly contains missing values.")
  }
  
  if (any(!is.finite(mat))) {
    stop("MOFA+ RNA matrix contains non-finite values.")
  }
  
  invisible(TRUE)
}

assert_scot_export <- function(mat, metadata) {
  
  if (!identical(rownames(mat), metadata$SampleID)) {
    stop("SCOT+ matrix row order does not match metadata SampleID order.")
  }
  
  if (anyNA(mat)) {
    stop("SCOT+ matrix contains missing values.")
  }
  
  if (any(!is.finite(mat))) {
    stop("SCOT+ matrix contains non-finite values.")
  }
  
  invisible(TRUE)
}

# =========================================================
# 4. Load and validate master metadata
# =========================================================

cat("Loading master metadata...\n")

master_metadata <- read_csv(metadata_file, show_col_types = FALSE)

required_meta_cols <- c("SampleID", "TreatmentCode", "Condition", "ConditionLabel", "Block", "WellIndex", "RNA_available",
                        "RNA_norm_sum", "RNA_detected_genes_gt0", "RNA_detected_genes_ge5_norm", "RNA_zero_fraction",
                        "Paired_RNA_ProteinQC")

missing_meta_cols <- setdiff(required_meta_cols, colnames(master_metadata))

if (length(missing_meta_cols) > 0) {
  stop("metadata_master_samples.csv is missing required column(s): ",
       paste(missing_meta_cols,collapse = ", "))
}

validate_unique_ids(master_metadata$SampleID, "Master metadata SampleID column")

if (!is.logical(master_metadata$RNA_available)) {
  stop("RNA_available must be imported as a logical TRUE/FALSE column.")
}

if (!is.logical(master_metadata$Paired_RNA_ProteinQC)) {
  stop("Paired_RNA_ProteinQC must be imported as a logical TRUE/FALSE column.")
}

rna_sample_ids <- master_metadata %>%
  filter(RNA_available) %>% 
  pull(SampleID)

paired_sample_ids <- master_metadata %>%
  filter(Paired_RNA_ProteinQC) %>%
  pull(SampleID)

validate_unique_ids(rna_sample_ids, "RNA sample-ID set")

validate_unique_ids(paired_sample_ids, "Paired sample-ID set")

if (!all(paired_sample_ids %in% rna_sample_ids)) {
  stop("The paired RNA/protein-QC sample set is not a subset of the RNA sample set.")
}

if (length(rna_sample_ids) != expected_n_rna) {
  stop("Expected ", expected_n_rna, " RNA samples but found ", length(rna_sample_ids),
       ". Check metadata and sample IDs before continuing.")
}

if (length(paired_sample_ids) != expected_n_paired) {
  stop("Expected ", expected_n_paired, " paired RNA/protein-QC samples but found ",
       length(paired_sample_ids), ". Check metadata and QC definitions before continuing.")
}

cat("RNA-available samples:               ", length(rna_sample_ids), "\n")

cat("Strict paired RNA/protein-QC samples:", length(paired_sample_ids), "\n\n")

# =========================================================
# 5. Load and validate RNA matrix
# =========================================================

cat("Loading RNA normalized-count matrix...\n")

rna_raw <- read.csv(rna_file, check.names = FALSE, stringsAsFactors = FALSE)

if (ncol(rna_raw) < 2) {
  stop("RNA file must contain one gene column and at least one sample column.")
}

# First column contains gene names
colnames(rna_raw)[1] <- "Gene"

# Check all column names before extracting sample columns.
# Do not use setdiff() here because it removes duplicated names.
if (anyDuplicated(colnames(rna_raw))) {
  duplicated_columns <- unique(colnames(rna_raw)[duplicated(colnames(rna_raw))])
  
  stop("Duplicated RNA column names were detected: ", paste(duplicated_columns, collapse = ", "))
}

if (anyNA(rna_raw$Gene)) {
  stop("RNA matrix contains missing gene names.")
}

if (any(!nzchar(trimws(as.character(rna_raw$Gene))))) {
  stop("RNA matrix contains blank or whitespace-only gene names.")
}

if (anyDuplicated(rna_raw$Gene)) {
  duplicated_original <- unique(rna_raw$Gene[duplicated(rna_raw$Gene)])
  
  stop("Duplicated gene names were detected in the RNA matrix. ",
       "Do not use make.unique(), because suffixed names would interfere ",
       "with later RNA-protein feature matching. ", "Example duplicated genes: ",
       paste(head(duplicated_original, 20), collapse = ", "))
}

rna_matrix_sample_ids <- colnames(rna_raw)[-1]

validate_unique_ids(rna_matrix_sample_ids, "RNA matrix sample columns")

rna_mat <- safe_numeric_matrix(rna_raw[, rna_matrix_sample_ids, drop = FALSE],
                               object_name = "RNA normalized-count matrix")

rownames(rna_mat) <- rna_raw$Gene

if (any(rna_mat < 0)) {
  stop("RNA normalized-count matrix contains negative values.")
}

missing_in_matrix <- setdiff(rna_sample_ids, colnames(rna_mat))

extra_in_matrix <- setdiff(colnames(rna_mat), rna_sample_ids)

if (length(missing_in_matrix) > 0) {
  stop("RNA metadata samples missing from the RNA matrix: ", paste(missing_in_matrix,
                                                                   collapse = ", "))
}

if (length(extra_in_matrix) > 0) {
  warning("RNA matrix sample columns not marked RNA_available in metadata: ", 
          paste(extra_in_matrix, collapse = ", "))
}

# Restrict and reorder matrix using metadata-defined RNA sample order
rna_mat <- rna_mat[ , rna_sample_ids, drop = FALSE]

if (!identical(colnames(rna_mat), rna_sample_ids)) {
  stop("RNA matrix sample order does not match the RNA metadata sample order.")
}

cat("RNA matrix dimensions before filtering: ", nrow(rna_mat), " genes x ", ncol(rna_mat),
    " samples\n")

integer_fraction <- mean(rna_mat == round(rna_mat))

cat("Fraction of matrix entries that are integer-valued: ", round(integer_fraction, 4), "\n\n")

# =========================================================
# 6. RNA sample QC
# =========================================================

rna_qc <- tibble(SampleID = colnames(rna_mat), RNA_norm_sum = colSums(rna_mat, na.rm = TRUE),
                 RNA_detected_genes_gt0 = colSums(rna_mat > 0, na.rm = TRUE),
                 RNA_detected_genes_ge5_norm = colSums(rna_mat >= 5, na.rm = TRUE),
                 RNA_zero_fraction = colMeans(rna_mat == 0, na.rm = TRUE)) %>%
  left_join(master_metadata %>%
              select(SampleID, TreatmentCode, Condition, ConditionLabel, Block,WellIndex,
                     Paired_RNA_ProteinQC), by = "SampleID")

rna_qc_validation <- rna_qc %>%
  select(SampleID,RNA_norm_sum_recalculated = RNA_norm_sum,
         RNA_detected_genes_gt0_recalculated = RNA_detected_genes_gt0,
         RNA_detected_genes_ge5_norm_recalculated = RNA_detected_genes_ge5_norm,
         RNA_zero_fraction_recalculated = RNA_zero_fraction) %>%
  left_join(master_metadata %>%
              select(SampleID, RNA_norm_sum_metadata = RNA_norm_sum,
                     RNA_detected_genes_gt0_metadata = RNA_detected_genes_gt0,
                     RNA_detected_genes_ge5_norm_metadata = RNA_detected_genes_ge5_norm,
                     RNA_zero_fraction_metadata = RNA_zero_fraction), by = "SampleID")

if (!isTRUE(all.equal(unname(rna_qc_validation$RNA_norm_sum_recalculated),
                      unname(rna_qc_validation$RNA_norm_sum_metadata), tolerance = 1e-8, check.attributes = FALSE))) {
  stop("Recalculated RNA normalized sums do not match metadata.")
}

if (!isTRUE(all.equal(unname(rna_qc_validation$RNA_detected_genes_gt0_recalculated),
                      unname(rna_qc_validation$RNA_detected_genes_gt0_metadata), tolerance = 0, check.attributes = FALSE))) {
  stop("Recalculated detected-gene counts do not match metadata.")
}

if (!isTRUE(all.equal(unname(rna_qc_validation$RNA_detected_genes_ge5_norm_recalculated),
                      unname(rna_qc_validation$RNA_detected_genes_ge5_norm_metadata), tolerance = 0, check.attributes = FALSE))) {
  stop("Recalculated normalized-expression >=5 counts do not match metadata.")
}

if (!isTRUE(all.equal(unname(rna_qc_validation$RNA_zero_fraction_recalculated),
                      unname(rna_qc_validation$RNA_zero_fraction_metadata), tolerance = 1e-12, check.attributes = FALSE))) {
  stop("Recalculated RNA zero fractions do not match metadata.")
}

if (anyNA(rna_qc$Condition)) {
  stop("One or more RNA samples lack condition metadata.")
}

if (!identical(rna_qc$SampleID, colnames(rna_mat))) {
  stop("RNA QC table sample order does not match the RNA matrix.")
}

write_csv(rna_qc, file.path(qc_dir, "rna_sample_qc_metrics.csv"))

cat("RNA sample QC summary:\n")

print(rna_qc %>% 
        group_by(Condition) %>%
        summarise(n = n(), median_norm_sum = median(RNA_norm_sum),
                  median_detected_genes_gt0 = median(RNA_detected_genes_gt0),
                  median_detected_genes_ge5_norm = median(RNA_detected_genes_ge5_norm),
                  median_zero_fraction = median(RNA_zero_fraction), .groups = "drop"))

# =========================================================
# 7. Cohort-specific gene QC and filtering
# =========================================================

# ---------------------------------------------------------
# 7a. All RNA cells
# ---------------------------------------------------------

gene_qc_all <- calculate_gene_qc(rna_mat, cohort_name = "RNA_all",
                                 min_detected_fraction = min_detected_fraction)

genes_keep_all <- gene_qc_all %>%
  filter(PassCombinedFilter) %>%
  pull(Gene)

if (length(genes_keep_all) == 0) {
  stop("No genes passed filtering in the all-RNA cohort.")
}

rna_filtered_all <- rna_mat[genes_keep_all, rna_sample_ids, drop = FALSE]

# ---------------------------------------------------------
# 7b. Strict paired benchmark
# ---------------------------------------------------------

rna_mat_paired <- rna_mat[ , paired_sample_ids, drop = FALSE]

gene_qc_paired <- calculate_gene_qc(rna_mat_paired, cohort_name = "RNA_paired",
                                    min_detected_fraction = min_detected_fraction)

genes_keep_paired <- gene_qc_paired %>%
  filter(PassCombinedFilter) %>%
  pull(Gene)

if (length(genes_keep_paired) == 0) {
  stop("No genes passed filtering in the paired benchmark cohort.")
}

rna_filtered_paired <- rna_mat_paired[genes_keep_paired, paired_sample_ids, drop = FALSE]

gene_qc_combined <- bind_rows(gene_qc_all, gene_qc_paired)

write_csv(gene_qc_combined, file.path(qc_dir, "rna_gene_qc_metrics_by_cohort.csv"))

write_lines(genes_keep_all, file.path(shared_dir, "rna_genes_after_detection_filter_all62.txt"))

write_lines(genes_keep_paired, file.path(shared_dir, "rna_genes_after_detection_filter_paired58.txt"))

cat("Genes retained in all-RNA cohort: ", nrow(rna_filtered_all), "\n")

cat("Genes retained in paired cohort:  ", nrow(rna_filtered_paired), "\n\n")

# =========================================================
# 8. LogCPM transformation
# =========================================================

# Calculate logCPM from the complete, unfiltered RNA matrix.
# This ensures that edgeR library sizes are calculated from all genes,
# rather than changing according to the cohort-specific gene filter.
#
# edgeR::cpm() returns log2 counts per million when log = TRUE.
rna_logcpm_full <- edgeR::cpm(rna_mat, log = TRUE, prior.count = logcpm_prior_count)

if (!identical(rownames(rna_logcpm_full), rownames(rna_mat))) {
  stop("Gene order changed unexpectedly during logCPM transformation.")
}

if (!identical(colnames(rna_logcpm_full), colnames(rna_mat))) {
  stop("Sample order changed unexpectedly during logCPM transformation.")
}

if (anyNA(rna_logcpm_full) || any(!is.finite(rna_logcpm_full))) {
  stop("Missing or non-finite values were produced during logCPM transformation.")
}

# Apply the cohort-specific gene and sample subsets after calculating
# logCPM from the complete matrix.
rna_logcpm_all <- rna_logcpm_full[genes_keep_all, rna_sample_ids, drop = FALSE]

rna_logcpm_paired <- rna_logcpm_full[genes_keep_paired, paired_sample_ids, drop = FALSE]

write_matrix(rna_logcpm_all, file.path(shared_dir, "rna_logcpm_all62.csv"), "Gene")

write_matrix(rna_logcpm_paired, file.path(shared_dir, "rna_logcpm_paired58.csv"), "Gene")

cat("Shared all-RNA logCPM matrix: ", nrow(rna_logcpm_all), " genes x ", ncol(rna_logcpm_all), " samples\n")

cat("Shared paired logCPM matrix:  ", nrow(rna_logcpm_paired), " genes x ", ncol(rna_logcpm_paired), " samples\n\n")

# =========================================================
# 9. Cohort-specific variable-gene selection
# =========================================================

selection_all <- select_top_variable_genes(rna_logcpm_all, n_variable_genes)

selection_paired <- select_top_variable_genes(rna_logcpm_paired, n_variable_genes)

write_lines(selection_all$selected, file.path(shared_dir, "rna_variable_genes_all62.txt"))

write_lines(selection_paired$selected, file.path(shared_dir, "rna_variable_genes_paired58.txt"))

all_feature_names <- union(rownames(rna_logcpm_all), rownames(rna_logcpm_paired))

variable_gene_table <- tibble(Gene = all_feature_names,
                              Variance_logCPM_all62 = selection_all$variance[
                                match(Gene, names(selection_all$variance))],
                              Variance_logCPM_paired58 = selection_paired$variance[
                                match(Gene, names(selection_paired$variance))],
                              Selected_all62 = Gene %in% selection_all$selected,
                              Selected_paired58 = Gene %in% selection_paired$selected) %>%
  arrange(desc(Variance_logCPM_all62))

write_csv(variable_gene_table, file.path(qc_dir, "rna_variable_gene_selection_table.csv"))

cat("Variable genes selected, all62:    ", length(selection_all$selected), "\n")

cat("Variable genes selected, paired58: ", length(selection_paired$selected), "\n\n")

# =========================================================
# 10. Scale selected RNA features
# =========================================================

rna_all_logcpm_selected <- rna_logcpm_all[selection_all$selected, rna_sample_ids, drop = FALSE]

rna_paired_logcpm_selected <- rna_logcpm_paired[selection_paired$selected, paired_sample_ids,
                                               drop = FALSE]

scaled_all <- zscore_rows(rna_all_logcpm_selected)
scaled_paired <- zscore_rows(rna_paired_logcpm_selected)

# Features x samples
rna_all_z <- scaled_all$matrix
rna_paired_z <- scaled_paired$matrix

write_lines(scaled_all$dropped_features, 
            file.path(qc_dir, "rna_all62_features_dropped_during_scaling.txt"))

write_lines(scaled_paired$dropped_features,
            file.path(qc_dir, "rna_paired58_features_dropped_during_scaling.txt"))

metadata_all <- order_metadata_to_samples(master_metadata, colnames(rna_all_z))

metadata_paired <- order_metadata_to_samples(master_metadata, colnames(rna_paired_z))

assert_mofa_export(rna_all_z, metadata_all)
assert_mofa_export(rna_paired_z, metadata_paired)

cat("Scaled all62 matrix:    ", nrow(rna_all_z), " genes x ", ncol(rna_all_z), " samples\n")

cat("Scaled paired58 matrix: ", nrow(rna_paired_z), " genes x ", ncol(rna_paired_z), " samples\n\n")

# =========================================================
# 11. MOFA+ and SCOT+ exports
# =========================================================

# Primary MOFA+ inputs: feature-wise z-scored matrices. LogCPM-transformed,
# unscaled matrices are also exported for sensitivity analyses and diagnostics.

# MOFA+ uses features x samples
write_matrix(rna_all_z, file.path(mofa_dir, "mofa_rna_all62_z.csv"), "Gene")

write_matrix(rna_paired_z, file.path(mofa_dir, "mofa_rna_paired58_z.csv"), "Gene")

write_matrix(rna_all_logcpm_selected, file.path(mofa_dir, "mofa_rna_all62_logcpm.csv"), "Gene")

write_matrix(rna_paired_logcpm_selected, file.path(mofa_dir, "mofa_rna_paired58_logcpm.csv"), "Gene")

write_csv(metadata_all, file.path(mofa_dir, "mofa_rna_all62_metadata.csv"))

write_csv(metadata_paired, file.path(mofa_dir, "mofa_rna_paired58_metadata.csv"))

# SCOT+ uses samples x features
scot_all <- t(rna_all_z)
scot_paired <- t(rna_paired_z)

assert_scot_export(scot_all, metadata_all)
assert_scot_export(scot_paired, metadata_paired)

write_matrix(scot_all, file.path(scot_dir, "scot_rna_all62_samples_by_features_z.csv"), "SampleID")

write_matrix(scot_paired, file.path(scot_dir, "scot_rna_paired58_samples_by_features_z.csv"),
             "SampleID")

write_csv(metadata_all, file.path(scot_dir, "scot_rna_all62_metadata.csv"))

write_csv(metadata_paired, file.path(scot_dir, "scot_rna_paired58_metadata.csv"))

cat("SCOT all62 matrix:    ", nrow(scot_all), " samples x ", ncol(scot_all), " genes\n")

cat("SCOT paired58 matrix: ", nrow(scot_paired), " samples x ", ncol(scot_paired), " genes\n\n")

# =========================================================
# 12. SCOT+-ready PCA embeddings
# =========================================================

cat("Computing SCOT RNA PCA embeddings...\n")

scot_pca_all <- make_scot_pca(scot_all, n_pcs = n_scot_pcs)

scot_pca_paired <- make_scot_pca(scot_paired, n_pcs = n_scot_pcs)

if (!identical(rownames(scot_pca_all$embedding), metadata_all$SampleID)) {
  stop("All-RNA PCA embedding sample order does not match metadata.")
}

if (!identical(rownames(scot_pca_paired$embedding), metadata_paired$SampleID)) {
  stop("Paired PCA embedding sample order does not match metadata.")
}

write_matrix(scot_pca_all$embedding, file.path(scot_dir, "scot_rna_all62_pca_embedding.csv"),
             "SampleID")

write_matrix(scot_pca_paired$embedding, file.path(scot_dir, "scot_rna_paired58_pca_embedding.csv"),
             "SampleID")

write_csv(scot_pca_all$variance_explained,
          file.path(scot_dir,"scot_rna_all62_pca_variance_explained.csv"))

write_csv(scot_pca_paired$variance_explained,
          file.path(scot_dir, "scot_rna_paired58_pca_variance_explained.csv"))

cat("SCOT all62 PCA embedding:    ", nrow(scot_pca_all$embedding), " samples x ",
    ncol(scot_pca_all$embedding), " PCs\n")

cat("SCOT paired58 PCA embedding: ", nrow(scot_pca_paired$embedding), " samples x ",
    ncol(scot_pca_paired$embedding), " PCs\n\n")

# =========================================================
# 13. Diagnostic PCA for paired RNA matrix
# =========================================================

cat("Preparing paired RNA PCA diagnostics...\n")

diagnostic_scores <- as.data.frame(scot_pca_paired$embedding, check.names = FALSE) %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(metadata_paired %>%
              select(SampleID, Condition, ConditionLabel, Block, WellIndex, RNA_norm_sum,
                     RNA_detected_genes_ge5_norm, RNA_zero_fraction), by = "SampleID")

if (anyNA(diagnostic_scores$Condition)) {
  stop("Diagnostic PCA scores could not be fully matched to metadata.")
}

write_csv(diagnostic_scores, file.path(qc_dir, "rna_paired58_pca_scores.csv"))

diagnostic_variance <- scot_pca_paired$variance_explained

write_csv(diagnostic_variance, file.path(qc_dir, "rna_paired58_pca_variance_explained.csv"))

# =========================================================
# 14. PCA correlations with technical metrics
# =========================================================

technical_metrics <- c("RNA_norm_sum", "RNA_detected_genes_ge5_norm", "RNA_zero_fraction")

available_pcs <- intersect(paste0("PC", 1:5), colnames(diagnostic_scores))

technical_correlations <- map_dfr(available_pcs, function(pc_name) {
  map_dfr(technical_metrics, function(metric_name) {
    tibble(PC = pc_name, Metric = metric_name, Correlation = safe_cor(diagnostic_scores[[pc_name]],
                                                                      diagnostic_scores[[metric_name]]))
    })
  })

write_csv(technical_correlations, file.path(qc_dir, "rna_pca_technical_metric_correlations.csv"))

cat("PCA correlations with technical metrics:\n")
print(technical_correlations)

# =========================================================
# 15. Condition and block PCA summaries
# =========================================================

condition_pca_summary <- diagnostic_scores %>%
  group_by(Condition) %>%
  summarise(n = n(), across(all_of(available_pcs), list(mean = ~ mean(.x, na.rm = TRUE),
                                                        sd = ~ sd(.x, na.rm = TRUE))),
            .groups = "drop")

write_csv(condition_pca_summary, file.path(qc_dir, "rna_pca_scores_summary_by_condition.csv"))

block_pca_summary <- diagnostic_scores %>%
  group_by(Block) %>%
  summarise(n = n(), across(all_of(available_pcs), list(mean = ~ mean(.x, na.rm = TRUE),
                                                        sd = ~ sd(.x, na.rm = TRUE))),
            .groups = "drop")

write_csv(block_pca_summary, file.path(qc_dir, "rna_pca_scores_summary_by_block.csv"))

# =========================================================
# 16. Optional interactive diagnostic plots
# =========================================================

if (interactive() && all(c("PC1", "PC2") %in% colnames(diagnostic_scores))) {
  plot(diagnostic_scores$PC1, diagnostic_scores$PC2,
       col = ifelse(diagnostic_scores$Condition == "Control", "blue", "red"), pch = 16,
       xlab = "PC1", ylab = "PC2", main = "RNA PCA: paired SCOT+ input")
  
  legend("topright", legend = c("Control", "Treated"), col = c("blue", "red"), pch = 16)
}

# =========================================================
# 17. Final PCA checks
# =========================================================

stopifnot(identical(rownames(scot_pca_all$embedding), metadata_all$SampleID))

stopifnot(identical(rownames(scot_pca_paired$embedding), metadata_paired$SampleID))

stopifnot(all(is.finite(scot_pca_all$embedding)))

stopifnot(all(is.finite(scot_pca_paired$embedding)))

# =========================================================
# 18. Summary report
# =========================================================

summary_report <- tibble(Item = c("Input RNA samples", "Input RNA genes", "Paired benchmark samples",
                                  "Minimum detected fraction", "logCPM prior count",
                                  "Genes after all62 filter", "Genes after paired58 filter",
                                  "Variable genes requested", "All62 scaled features",
                                  "Paired58 scaled features", "SCOT all62 samples",
                                  "SCOT paired58 samples", "SCOT PCs exported"),
                         Value = as.character(c(ncol(rna_mat), nrow(rna_mat), length(paired_sample_ids),
                                                min_detected_fraction, logcpm_prior_count,
                                                nrow(rna_logcpm_all), nrow(rna_logcpm_paired),
                                                n_variable_genes, nrow(rna_all_z), nrow(rna_paired_z),
                                                nrow(scot_all), nrow(scot_paired),
                                                ncol(scot_pca_all$embedding))))

write_csv(summary_report, file.path(qc_dir, "rna_preprocessing_summary.csv"))

# Record parameters separately
parameter_report <- tibble(Parameter = c("min_detected_fraction", "n_variable_genes", "n_scot_pcs",
                                         "logcpm_prior_count"),
                           Value = as.character(c(min_detected_fraction, n_variable_genes,
                                                  n_scot_pcs, logcpm_prior_count)))

write_csv(parameter_report, file.path(qc_dir, "rna_preprocessing_parameters.csv"))

write_lines(capture.output(sessionInfo()), file.path(qc_dir, "rna_preprocessing_sessionInfo.txt"))

# =========================================================
# 19. Completion message
# =========================================================

cat("\nTranscriptomics preprocessing complete.\n\n")

cat("Primary shared outputs:\n")
cat("  ", file.path(shared_dir, "rna_logcpm_all62.csv"), "\n")
cat("  ", file.path(shared_dir, "rna_logcpm_paired58.csv"), "\n\n")

cat("Primary MOFA+ outputs:\n")
cat("  ", file.path(mofa_dir, "mofa_rna_all62_z.csv"), "\n")
cat("  ", file.path(mofa_dir, "mofa_rna_paired58_z.csv"), "\n\n")

cat("Primary SCOT+ feature-matrix outputs:\n")
cat("  ", file.path(scot_dir, "scot_rna_all62_samples_by_features_z.csv"), "\n")
cat("  ", file.path(scot_dir, "scot_rna_paired58_samples_by_features_z.csv"), "\n\n")

cat("Primary SCOT+ PCA outputs:\n")
cat("  ", file.path(scot_dir, "scot_rna_all62_pca_embedding.csv"), "\n")
cat("  ", file.path(scot_dir, "scot_rna_paired58_pca_embedding.csv"), "\n\n")

cat("QC summary:\n")
cat("  ", file.path(qc_dir, "rna_preprocessing_summary.csv"), "\n")
