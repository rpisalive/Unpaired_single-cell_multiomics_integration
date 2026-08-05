# =========================================================
# Metadata preparation for NanoSPLITS C10 RO-3306 experiment
# =========================================================
# Purpose:
#   Starting from the three NanoSPLITS matrices:
#     1. C10Treated_NormCounts.csv
#     2. C10Treated_Protein_intensities.tsv
#     3. C10Treated_Peptide_intensities.tsv
#
#   Prepare:
#     - RNA metadata
#     - Protein metadata
#     - Peptide / phosphopeptide metadata
#     - Master sample metadata
#     - Sample-ID lists for downstream preprocessing, SCOT+, and MOFA+
#
# =========================================================
# Package versions
# =========================================================
# R version recommended: >= 4.4
# tidyverse >= 2.0.0
# readr >= 2.1.5
# rstudioapi >= 0.17.1

library(tidyverse)
library(readr)
library(rstudioapi)

# =========================================================
# 1. Path setup
# =========================================================

get_script_dir <- function() {
  ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
  
  if (!is.null(ctx) && !is.null(ctx$path) && nzchar(ctx$path)) {
    return(dirname(ctx$path))
  } else {
    return(getwd())
  }
}

script_dir <- get_script_dir()

# Raw data directory supplied by user
raw_data_dir <- "C:/Users/49152/Downloads/2nd_paper/Github_data"

# Output directory
processed_dir <- file.path(script_dir, "processed_data")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# Input files
rna_file <- file.path(raw_data_dir, "C10Treated_NormCounts.csv")
protein_file <- file.path(raw_data_dir, "C10Treated_Protein_intensities.tsv")
peptide_file <- file.path(raw_data_dir, "C10Treated_Peptide_intensities.tsv")

# Check that files exist
stopifnot(file.exists(rna_file))
stopifnot(file.exists(protein_file))
stopifnot(file.exists(peptide_file))

cat("Script directory:   ", script_dir, "\n")
cat("Raw data directory:", raw_data_dir, "\n")
cat("Output directory:  ", processed_dir, "\n\n")

# =========================================================
# 2. Helper functions
# =========================================================

parse_condition <- function(sample_id) {
  case_when(str_detect(sample_id, "^Ctrl[A-Za-z]+\\d+$") ~ "Control",
            str_detect(sample_id, "^Trtdt?[A-Za-z]+\\d+$") ~ "Treated", TRUE ~ NA_character_)
}

parse_condition_label <- function(sample_id) {
  case_when(str_detect(sample_id, "^Ctrl[A-Za-z]+\\d+$") ~ "Untreated",
            str_detect(sample_id, "^Trtdt?[A-Za-z]+\\d+$") ~ "G2M_arrested", TRUE ~ NA_character_)
}

parse_treatment_code <- function(sample_id) {
  case_when(str_detect(sample_id, "^Ctrl[A-Za-z]+\\d+$") ~ "Ctrl",
            str_detect(sample_id, "^Trtdt?[A-Za-z]+\\d+$") ~ "Trtd", TRUE ~ NA_character_)
}

# Handles sample names such as:
#   CtrlA1
#   CtrlA10
#   TrtdA1
#   TrtdD9
#   TrtdtA11  # apparent typo in original file; handled without renaming
strip_condition_prefix <- function(sample_id) {
  str_remove(sample_id, "^(Ctrl|Trtdt?)")
}

parse_block <- function(sample_id) {
  stripped <- strip_condition_prefix(sample_id)
  str_extract(stripped, "^[A-Za-z]+")
}

parse_well_index <- function(sample_id) {
  stripped <- strip_condition_prefix(sample_id)
  as.integer(str_extract(stripped, "\\d+"))
}

make_basic_sample_metadata <- function(sample_ids) {
  tibble(SampleID = sample_ids, TreatmentCode = parse_treatment_code(sample_ids),
         Condition = parse_condition(sample_ids), ConditionLabel = parse_condition_label(sample_ids),
         Block = parse_block(sample_ids), WellIndex = parse_well_index(sample_ids))
}

safe_numeric_matrix <- function(df, matrix_name = "matrix") {
  
  # Record missing values already present before conversion.
  original_mat <- as.matrix(df)
  original_na <- is.na(original_mat)
  
  # Convert every value to numeric.
  numeric_values <- suppressWarnings(as.numeric(original_mat))
  
  mat <- matrix(numeric_values, nrow = nrow(df), ncol = ncol(df), dimnames = list(rownames(df), colnames(df)))
  
  # Identify NA values created specifically by numeric conversion.
  introduced_na <- is.na(mat) & !original_na
  
  if (any(introduced_na)) {
    stop(matrix_name, " contains non-numeric values in presumed sample columns. ", sum(introduced_na),
         " value(s) were converted to NA.")
  }
  
  if (any(is.infinite(mat))) {
    stop(matrix_name, " contains infinite values.")
  }
  
  if (any(mat < 0, na.rm = TRUE)) {
    stop(matrix_name, " contains negative values.")
  }
  
  mat
}

# =========================================================
# 3. RNA metadata
# =========================================================

cat("Loading RNA matrix...\n")

rna_raw <- read.csv(rna_file, check.names = FALSE, stringsAsFactors = FALSE)

if (ncol(rna_raw) < 2) {
  stop("RNA file must contain one gene column and at least one sample column.")
}

# The first column contains gene identifiers.
colnames(rna_raw)[1] <- "Gene"

# Check column names before extracting sample columns.
if (anyDuplicated(colnames(rna_raw))) {
  duplicated_ids <- unique(colnames(rna_raw)[duplicated(colnames(rna_raw))])
  
  stop("Duplicated RNA column names detected: ", paste(duplicated_ids, collapse = ", "))
}

# Validate gene identifiers.
if (anyNA(rna_raw$Gene) || any(!nzchar(trimws(as.character(rna_raw$Gene))))) {
  stop("RNA matrix contains missing or empty gene identifiers.")
}

if (anyDuplicated(rna_raw$Gene)) {
  duplicated_genes <- unique(rna_raw$Gene[duplicated(rna_raw$Gene)])
  
  stop("RNA matrix contains duplicated gene identifiers. Examples: ", paste(head(duplicated_genes, 20), collapse = ", "))
}

# All columns after the first column are RNA sample columns.
rna_sample_cols <- colnames(rna_raw)[-1]

rna_mat <- rna_raw[, rna_sample_cols, drop = FALSE]
rna_mat <- safe_numeric_matrix(rna_mat, "RNA matrix")

# RNA input is expected to be complete.
if (anyNA(rna_mat)) {
  stop("RNA matrix contains missing values in sample columns.")
}

rownames(rna_mat) <- rna_raw$Gene

rna_metadata <- make_basic_sample_metadata(rna_sample_cols) %>%
  mutate(RNA_available = TRUE, RNA_norm_sum = colSums(rna_mat, na.rm = TRUE),
         RNA_detected_genes_gt0 = colSums(rna_mat > 0, na.rm = TRUE),
         RNA_detected_genes_ge5_norm = colSums(rna_mat >= 5, na.rm = TRUE),
         RNA_zero_fraction = colMeans(rna_mat == 0, na.rm = TRUE))

cat("RNA samples:", nrow(rna_metadata), "\n")
print(table(rna_metadata$Condition, useNA = "ifany"))

# =========================================================
# 4. Protein metadata
# =========================================================

cat("Loading protein intensity matrix...\n")

protein_raw <- read_tsv(protein_file, show_col_types = FALSE, progress = FALSE, name_repair = "minimal")

# readr otherwise repairs duplicated names automatically.
if (anyDuplicated(colnames(protein_raw))) {
  duplicated_ids <- unique(colnames(protein_raw)[duplicated(colnames(protein_raw))])
  
  stop("Duplicated protein column names detected: ", paste(duplicated_ids, collapse = ", "))
}

required_protein_annot <- c("Protein", "Protein.ID", "Entry.Name", "Gene")

missing_protein_annot <- setdiff(required_protein_annot, colnames(protein_raw))

if (length(missing_protein_annot) > 0) {
  stop("Protein file is missing expected annotation column(s): ", paste(missing_protein_annot, collapse = ", "))
}

# Preserve the original column order.
protein_sample_cols <- colnames(protein_raw)[!colnames(protein_raw) %in% required_protein_annot]

if (length(protein_sample_cols) == 0) {
  stop("No protein sample columns were detected.")
}

# Keep mouse proteins for QC.
# Add a persistent row identifier before filtering so that feature metadata
# can always be mapped back to the corresponding intensity-matrix row.
# The NanoSPLITS author QC uses Entry.Name as the protein identifier.
protein_mouse <- protein_raw %>%
  mutate(ProteinRowID = sprintf("PROT_%06d", row_number())) %>%
  filter(!is.na(Entry.Name), grepl("_MOUSE", Entry.Name))

if (nrow(protein_mouse) == 0) {
  stop("No mouse protein rows were detected using Entry.Name.")
}

protein_mat <- protein_mouse[, protein_sample_cols, drop = FALSE]
protein_mat <- safe_numeric_matrix(protein_mat, "protein matrix")

# Count distinct detected Entry.Name identifiers per sample.
# Zero means not detected; NA is also treated as missing.
protein_detected_n <- vapply(seq_len(ncol(protein_mat)), function(j) {
    
    detected <- (!is.na(protein_mat[, j]) & protein_mat[, j] > 0)
    
    dplyr::n_distinct(protein_mouse$Entry.Name[detected], na.rm = TRUE)
  },
  integer(1))

names(protein_detected_n) <- colnames(protein_mat)

# Total number of distinct mouse protein identifiers.
protein_total_n <- dplyr::n_distinct(protein_mouse$Entry.Name, na.rm = TRUE)

protein_metadata <- make_basic_sample_metadata(protein_sample_cols) %>%
  mutate(Protein_available = TRUE,
    
    # Match detection counts explicitly by sample name.
    Protein_detected_n = unname(protein_detected_n[SampleID]),
    
    Protein_missing_n = protein_total_n - Protein_detected_n,
    
    Protein_missing_fraction = Protein_missing_n / protein_total_n,
    
    Protein_QC_pass_author = case_when(Condition == "Control" & Protein_detected_n >= 2000 ~ TRUE,
                                       Condition == "Treated" & Protein_detected_n >= 2500 ~ TRUE,
                                       TRUE ~ FALSE))

cat("Protein samples:", nrow(protein_metadata), "\n")
print(table(protein_metadata$Condition, useNA = "ifany"))
cat("Protein author-QC samples:", sum(protein_metadata$Protein_QC_pass_author), "\n")
print(table(protein_metadata$Condition, protein_metadata$Protein_QC_pass_author, useNA = "ifany"))

# =========================================================
# 5. Peptide and phosphopeptide metadata
# =========================================================

cat("Loading peptide intensity matrix...\n")

peptide_raw <- read_tsv(peptide_file, show_col_types = FALSE, progress = FALSE, name_repair = "minimal")

if (anyDuplicated(colnames(peptide_raw))) {
  duplicated_ids <- unique(colnames(peptide_raw)[duplicated(colnames(peptide_raw))])
  
  stop("Duplicated peptide column names detected: ", paste(duplicated_ids, collapse = ", "))
}

required_peptide_annot <- c("Peptide.Sequence", "Modified.Sequence", "Assigned.Modifications", "Protein", "ID", "Entry.Name",
                            "Gene", "Protein.Description")

missing_peptide_annot <- setdiff(required_peptide_annot, colnames(peptide_raw))

if (length(missing_peptide_annot) > 0) {
  stop("Peptide file is missing expected annotation column(s): ", paste(missing_peptide_annot, collapse = ", "))
}

# Preserve the original column order.
peptide_sample_cols <- colnames(peptide_raw)[!colnames(peptide_raw) %in% required_peptide_annot]

if (length(peptide_sample_cols) == 0) {
  stop("No peptide sample columns were detected.")
}

protein_only_samples <- setdiff(protein_sample_cols, peptide_sample_cols)

peptide_only_samples <- setdiff(peptide_sample_cols, protein_sample_cols)

if (length(protein_only_samples) > 0 || length(peptide_only_samples) > 0) {
  stop("Protein and peptide sample columns do not match.\n", "Protein-only samples: ",
       ifelse(length(protein_only_samples) == 0, "none", paste(protein_only_samples, collapse = ", ")),
       "\nPeptide-only samples: ", ifelse(length(peptide_only_samples) == 0, "none",
                                          paste(peptide_only_samples, collapse = ", ")))
}

if (anyDuplicated(peptide_sample_cols)) {
  duplicated_ids <- unique(peptide_sample_cols[duplicated(peptide_sample_cols)])
  
  stop("Duplicated peptide sample column names detected: ", paste(duplicated_ids, collapse = ", "))
}

# All mouse peptides.
peptide_mouse <- peptide_raw %>%
  mutate(PeptideRowID = sprintf("PEP_%06d", row_number())) %>%
  filter(grepl("_MOUSE", Entry.Name))

if (nrow(peptide_mouse) == 0) {
  stop("No mouse peptide rows were detected using Entry.Name.")
}

# Phosphopeptides only
# Phosphorylation appears as a ~79.9663 Da mass shift.
phosphopeptide_mouse <- peptide_mouse %>%
  filter(grepl("79\\.9", Modified.Sequence))

if (nrow(phosphopeptide_mouse) == 0) {
  stop("No phosphopeptide rows were detected using the 79.9 mass-shift pattern.")
}

# Use protein sample order so that protein and peptide matrices have
# identical sample-column ordering.
peptide_mat <- peptide_mouse[, protein_sample_cols, drop = FALSE]
peptide_mat <- safe_numeric_matrix(peptide_mat, "peptide matrix")

phosphopeptide_mat <- phosphopeptide_mouse[ , protein_sample_cols, drop = FALSE]
phosphopeptide_mat <- safe_numeric_matrix(phosphopeptide_mat, "phosphopeptide matrix")

peptide_detected_n <- colSums(peptide_mat > 0, na.rm = TRUE)
phosphopeptide_detected_n <- colSums(phosphopeptide_mat > 0, na.rm = TRUE)

peptide_metadata <- make_basic_sample_metadata(protein_sample_cols) %>%
  mutate(Peptide_available = TRUE,
    
    # Match detection counts explicitly by sample name.
    Peptide_detected_n = unname(peptide_detected_n[SampleID]),
    
    Peptide_missing_n = nrow(peptide_mat) - Peptide_detected_n,
    
    Peptide_missing_fraction = Peptide_missing_n / nrow(peptide_mat),
    
    Phosphopeptide_available = TRUE,
    
    Phosphopeptide_detected_n = unname(phosphopeptide_detected_n[SampleID]),
    
    Phosphopeptide_missing_n = nrow(phosphopeptide_mat) - Phosphopeptide_detected_n,
    
    Phosphopeptide_missing_fraction = Phosphopeptide_missing_n / nrow(phosphopeptide_mat))

cat("Peptide samples:", nrow(peptide_metadata), "\n")
print(table(peptide_metadata$Condition, useNA = "ifany"))
cat("Mouse peptide rows:", nrow(peptide_mouse), "\n")
cat("Mouse phosphopeptide rows:", nrow(phosphopeptide_mouse), "\n\n")


# =========================================================
# Validate sample-ID parsing
# =========================================================

all_input_sample_ids <- unique(c(rna_sample_cols, protein_sample_cols, peptide_sample_cols))

parsed_input_ids <- make_basic_sample_metadata(all_input_sample_ids)

bad_sample_ids <- parsed_input_ids %>%
  filter(is.na(TreatmentCode) | is.na(Condition) | is.na(ConditionLabel) | is.na(Block) | is.na(WellIndex)) %>%
  pull(SampleID)

if (length(bad_sample_ids) > 0) {
  stop("The following sample IDs could not be parsed: ", paste(bad_sample_ids, collapse = ", "))
}

# =========================================================
# 6. Master sample metadata
# =========================================================

cat("Constructing master metadata table...\n")

all_sample_ids <- sort(unique(c(rna_metadata$SampleID, protein_metadata$SampleID,
                                peptide_metadata$SampleID)))

master_metadata <- make_basic_sample_metadata(all_sample_ids) %>%
  left_join(rna_metadata %>%
              select(SampleID, RNA_available, RNA_norm_sum, RNA_detected_genes_gt0,
                     RNA_detected_genes_ge5_norm, RNA_zero_fraction), by = "SampleID") %>%
  left_join(protein_metadata %>%
              select(SampleID, Protein_available, Protein_detected_n, Protein_missing_n,
                     Protein_missing_fraction, Protein_QC_pass_author), by = "SampleID") %>%
  left_join(peptide_metadata %>%
              select(SampleID, Peptide_available, Peptide_detected_n, Peptide_missing_n,
                     Peptide_missing_fraction, Phosphopeptide_available, Phosphopeptide_detected_n,
                     Phosphopeptide_missing_n, Phosphopeptide_missing_fraction), by = "SampleID") %>%
  mutate(RNA_available = replace_na(RNA_available, FALSE),
         Protein_available = replace_na(Protein_available, FALSE),
         Peptide_available = replace_na(Peptide_available, FALSE),
         Phosphopeptide_available = replace_na(Phosphopeptide_available, FALSE),
         Protein_QC_pass_author = replace_na(Protein_QC_pass_author, FALSE),
         
         # Core sample sets for downstream analysis
         RNA62 = RNA_available,
         ProteinQC68 = Protein_available & Protein_QC_pass_author,
         
         # Peptide and phosphopeptide data from cells that passed protein-level QC.
         # These are not peptide-specific QC definitions.
         Peptide_on_ProteinQC68 = Peptide_available & Protein_QC_pass_author,
         Phosphopeptide_on_ProteinQC68 = Phosphopeptide_available & Protein_QC_pass_author,
         
         # Strict paired reference set for MOFA+ and SCOT+ benchmark
         Paired_RNA_ProteinQC = RNA_available & Protein_available & Protein_QC_pass_author,
         
         # Since peptide/phosphopeptide samples follow proteomics sample columns,
         # this is the practical paired multiomic set.
         Paired_RNA_ProteinQC_Peptide = RNA_available & Protein_available & Peptide_available &
           Protein_QC_pass_author) %>%
  
  arrange(Condition, Block, WellIndex, SampleID)

cat("Master samples:", nrow(master_metadata), "\n")

if (anyDuplicated(master_metadata$SampleID)) {
  duplicated_ids <- unique(master_metadata$SampleID[duplicated(master_metadata$SampleID)])
  
  stop("Duplicated SampleID values detected in master metadata: ",
       paste(duplicated_ids, collapse = ", "))
}

# Validate the expected sample sets for this specific experiment.
expected_counts <- c(Master_samples = 96L, RNA_available = 62L, Protein_available = 96L, Peptide_available = 96L,
                     Phosphopeptide_available = 96L, Protein_QC = 68L, Paired_RNA_ProteinQC = 58L,
                     Paired_RNA_ProteinQC_Peptide = 58L)

observed_counts <- c(Master_samples = nrow(master_metadata), RNA_available = sum(master_metadata$RNA_available),
                     Protein_available = sum(master_metadata$Protein_available),
                     Peptide_available = sum(master_metadata$Peptide_available),
                     Phosphopeptide_available = sum(master_metadata$Phosphopeptide_available),
                     Protein_QC = sum(master_metadata$ProteinQC68), 
                     Paired_RNA_ProteinQC = sum(master_metadata$Paired_RNA_ProteinQC),
                     Paired_RNA_ProteinQC_Peptide = sum(master_metadata$Paired_RNA_ProteinQC_Peptide))

incorrect_counts <- names(observed_counts)[observed_counts != expected_counts]

if (length(incorrect_counts) > 0) {
  count_messages <- vapply(incorrect_counts, function(x) {
      paste0(x, ": observed = ", observed_counts[[x]], ", expected = ", expected_counts[[x]])
    },
    character(1))
  
  stop("Unexpected experiment sample counts:\n", paste(count_messages, collapse = "\n"))
}

# Validate the condition composition of the author-QC set.
qc_condition_counts <- table(factor(master_metadata$Condition[master_metadata$ProteinQC68],
                                    levels = c("Control", "Treated")))

qc_condition_counts <- setNames(as.integer(qc_condition_counts), c("Control", "Treated"))

expected_qc_condition_counts <- c(Control = 32L, Treated = 36L)

if (any(qc_condition_counts != expected_qc_condition_counts)) {
  stop("Unexpected condition composition in ProteinQC68. ", "Observed: Control = ", qc_condition_counts[["Control"]],
       ", Treated = ", qc_condition_counts[["Treated"]], ". Expected: Control = 32, Treated = 36.")
}

cat("Master samples:", nrow(master_metadata), "\n")
cat("RNA available:", sum(master_metadata$RNA_available), "\n")
cat("Protein available:", sum(master_metadata$Protein_available), "\n")
cat("Protein author-QC:", sum(master_metadata$ProteinQC68), "\n")
cat("Peptide available:", sum(master_metadata$Peptide_available), "\n")
cat("Strict RNA + protein-QC paired:", sum(master_metadata$Paired_RNA_ProteinQC), "\n")
cat("Strict RNA + protein-QC + peptide paired:", sum(master_metadata$Paired_RNA_ProteinQC_Peptide), "\n\n")

cat("Condition composition of strict paired set:\n")
print(table(master_metadata$Condition, master_metadata$Paired_RNA_ProteinQC, useNA = "ifany"))

# =========================================================
# 7. Peptide feature metadata
# =========================================================

# Save peptide/phosphopeptide feature metadata for downstream processing.
# This is not sample metadata, but it will be needed for phosphosite analysis.

peptide_feature_metadata <- peptide_mouse %>%
  select(PeptideRowID, all_of(required_peptide_annot)) %>%
  mutate(IsPhosphopeptide = grepl("79\\.9", Modified.Sequence))

phosphopeptide_feature_metadata <- peptide_feature_metadata %>%
  filter(IsPhosphopeptide)

protein_feature_metadata <- protein_mouse %>%
  select(ProteinRowID, all_of(required_protein_annot))

# =========================================================
# 8. Export metadata and sample sets
# =========================================================

cat("Writing metadata files...\n")

write_csv(rna_metadata, file.path(processed_dir, "metadata_rna_samples.csv"))

write_csv(protein_metadata, file.path(processed_dir, "metadata_protein_samples.csv"))

write_csv(peptide_metadata, file.path(processed_dir, "metadata_peptide_samples.csv"))

write_csv(master_metadata, file.path(processed_dir, "metadata_master_samples.csv"))

write_csv(protein_feature_metadata, file.path(processed_dir, "metadata_protein_features.csv"))

write_csv(peptide_feature_metadata, file.path(processed_dir, "metadata_peptide_features.csv"))

write_csv(phosphopeptide_feature_metadata, file.path(processed_dir, "metadata_phosphopeptide_features.csv"))

# Sample ID lists for downstream scripts
write_lines(master_metadata %>%
              filter(RNA62) %>%
              pull(SampleID), file.path(processed_dir, "sample_ids_rna62.txt"))

write_lines(master_metadata %>%
              filter(ProteinQC68) %>%
              pull(SampleID), file.path(processed_dir, "sample_ids_protein_qc68.txt"))

write_lines(master_metadata %>%
              filter(Paired_RNA_ProteinQC) %>%
              pull(SampleID), file.path(processed_dir, "sample_ids_paired_rna_protein_qc.txt"))

write_lines(master_metadata %>%
              filter(Paired_RNA_ProteinQC_Peptide) %>%
              pull(SampleID), file.path(processed_dir, "sample_ids_paired_rna_protein_qc_peptide.txt"))

# =========================================================
# 9. Sanity checks
# =========================================================

cat("\nSanity checks:\n")

cat("RNA sample IDs not in protein file:\n")
print(setdiff(rna_metadata$SampleID, protein_metadata$SampleID))

cat("Protein sample IDs not in RNA file:\n")
print(setdiff(protein_metadata$SampleID, rna_metadata$SampleID))

cat("\nProtein sample IDs not in peptide file:\n")
print(setdiff(protein_metadata$SampleID, peptide_metadata$SampleID))

cat("\nPeptide sample IDs not in protein file:\n")
print(setdiff(peptide_metadata$SampleID, protein_metadata$SampleID))

cat("\nPotential unusual sample IDs:\n")
print(master_metadata %>%
        filter(is.na(Condition) | is.na(Block) | is.na(WellIndex)) %>%
        select(SampleID, Condition, Block, WellIndex))

cat("\nMetadata preparation complete.\n")
cat("Files written to:\n", processed_dir, "\n")
