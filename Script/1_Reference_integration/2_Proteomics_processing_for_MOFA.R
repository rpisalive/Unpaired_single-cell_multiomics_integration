library(tidyverse)
library(readr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)


setwd("PATH_TO_WORKING_DIRECTORY")

#Read MS intensities
raw <- read_tsv("PATH/TO/C10SVEC_singlecells_Protein_intensities.tsv", 
                show_col_types = FALSE)

#Remove contaminants
raw <- raw[!grepl("contam_sp", raw$PROTID), ]

#Extract sample columns
sample_cols <- setdiff(colnames(raw), "PROTID")

#Replace zeros with NA
raw[sample_cols] <- lapply(raw[sample_cols], function(x) ifelse(x == 0, NA, x))

#Step 2: Total-intensity normalization per cell
expr_mat <- as.matrix(raw[, sample_cols])
rownames(expr_mat) <- raw$PROTID

lib_sizes <- colSums(expr_mat, na.rm = TRUE)
median_lib <- median(lib_sizes, na.rm = TRUE)

#log-normalize (variance-stabilized)
log_norm <- log2(t(t(expr_mat) / lib_sizes * median_lib) + 0.1)

#Convert rownames to a column for mapping
log_norm_df <- as.data.frame(log_norm) %>%
  rownames_to_column(var = "PROTID")

log_norm_long <- log_norm_df %>%
  pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity") %>%
  filter(!is.na(Intensity))

#PROTID → Gene symbol mapping

UniAc <- unique(log_norm_long$PROTID)

geneSymbols <- select(org.Mm.eg.db, keys=UniAc, columns= c("SYMBOL","UNIPROT"),
                      keytype="UNIPROT", multiVals="first")

geneSymbols <- geneSymbols %>%
  mutate(Gene = SYMBOL) %>%
  dplyr::select(Gene,UNIPROT) 

convert <- read_tsv("PATH/TO/convert_uniprot.tsv")

geneSymbols <- full_join(geneSymbols, convert, by = c("Gene","UNIPROT")) %>%
  filter(!is.na(Gene))

log_norm_long <- log_norm_long %>%
  filter(!is.na(Intensity)) %>%
  mutate(UNIPROT = PROTID)

log_norm_long <- full_join(log_norm_long,geneSymbols, by = "UNIPROT") %>%
  mutate(Gene = case_when(is.na(Gene) ~ UNIPROT,
                          TRUE ~ Gene)) %>%
  group_by(SampleID, Gene) %>%
  slice_max(order_by = Intensity, n = 1 , with_ties = FALSE)  %>%
  ungroup()

expr_mat <- log_norm_long %>%
  dplyr::select(SampleID, Gene, Intensity) %>%
  pivot_wider(names_from = SampleID, values_from = Intensity) %>%
  column_to_rownames(var = "Gene")

#Convert to matrix
expr_mat <- as.matrix(expr_mat)

#Metadata
sample_ids <- colnames(expr_mat)
coldata <- data.frame(
  SampleID = sample_ids,
  Chip = substr(sample_ids, 1, 3),
  CellLine = ifelse(grepl("SVEC", sample_ids), "SVEC",
                    ifelse(grepl("C10", sample_ids), "C10", NA_character_)),
  stringsAsFactors = FALSE
)

#Filtering
detected_per_sample <- colSums(!is.na(expr_mat))
keep_samples <- detected_per_sample >= 1000
expr_mat <- expr_mat[, keep_samples, drop = FALSE]
coldata <- coldata[coldata$SampleID %in% colnames(expr_mat), ]

keep_genes <- rowMeans(!is.na(expr_mat)) >= 0.10
expr_mat <- expr_mat[keep_genes, , drop = FALSE]

#Feature selection
feature_var <- apply(expr_mat, 1, var, na.rm = TRUE)
feature_var <- feature_var[!is.na(feature_var) & feature_var > 0]

N_proteins <- min(1000, length(feature_var))
top_features <- names(sort(feature_var, decreasing = TRUE))[seq_len(N_proteins)]
prot_matrix <- expr_mat[top_features, , drop = FALSE]

#Export for MOFA+
write_csv(data.frame(Gene = rownames(prot_matrix), prot_matrix, check.names = FALSE),
          "PATH/TO/proteomics_for_MOFA.csv")

write_csv(coldata,
          "PATH/TO/proteomics_metadata.csv")

#PCA
prot_imputed <- prot_matrix
for (i in 1:nrow(prot_imputed)) {
  row_vals <- prot_imputed[i, ]
  row_mean <- mean(row_vals, na.rm = TRUE)
  row_vals[is.na(row_vals)] <- row_mean
  prot_imputed[i, ] <- row_vals
}

pca <- prcomp(t(prot_imputed), scale. = FALSE)
plot(pca$x[,1], pca$x[,2],
     col = ifelse(coldata$CellLine == "C10", "red", "blue"),
     pch = 16,
     xlab = "PC1", ylab = "PC2",
     main = "Proteomics PCA (imputed): C10 (red) vs SVEC (blue)")