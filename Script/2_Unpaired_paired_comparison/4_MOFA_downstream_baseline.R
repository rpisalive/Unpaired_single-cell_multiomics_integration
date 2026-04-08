library(ggplot2)
library(dplyr)
library(MOFA2)
library(svglite)
library(ggokabeito)
library(extrafont)
library(readr)

setwd("C:/Users/49152/Downloads/Multi-omics/MOFA/")

#1 Model loading
model <- load_model("output/Unpaired_single_cell_model.hdf5")

outdir <- "C:/Users/49152/Downloads/Multi-omics/MOFA/output/graphs/Unpaired_single_cell/Minimal_baseline/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

loadfonts(device = "win")
theme_set(theme_bw(base_size = 16, base_family = "Arial"))

#2 Metadata loading and alignment
rna_metadata <- read.csv("input/4_Unpaired_sc/transcriptomics_metadata.csv",
                         stringsAsFactors = FALSE)

prot_metadata <- read.csv("input/4_Unpaired_sc/proteomics_metadata.csv",
                          stringsAsFactors = FALSE)

if ("CellLines" %in% colnames(rna_metadata) && !"CellLine" %in% colnames(rna_metadata)) {
  rna_metadata <- rna_metadata %>% rename(CellLine = CellLines)
}
if ("CellLines" %in% colnames(prot_metadata) && !"CellLine" %in% colnames(prot_metadata)) {
  prot_metadata <- prot_metadata %>% rename(CellLine = CellLines)
}

rna_metadata$SampleID  <- paste0(as.character(rna_metadata$SampleID), "_rna")
prot_metadata$SampleID <- paste0(as.character(prot_metadata$SampleID), "_prot")

rna_metadata$Modality  <- "Transcriptomics"
prot_metadata$Modality <- "Proteomics"

metadata <- bind_rows(prot_metadata, rna_metadata) %>%
  rename(sample = SampleID)

model_samples <- unname(unlist(samples_names(model)))

metadata <- metadata %>%
  filter(sample %in% model_samples) %>%
  distinct(sample, .keep_all = TRUE)

metadata <- metadata[match(model_samples, metadata$sample), , drop = FALSE]

if (nrow(metadata) != length(model_samples)) {
  stop("Mismatch between model samples and metadata rows after alignment.")
}
if (!all(metadata$sample == model_samples)) {
  stop("Sample order mismatch between model and metadata.")
}

model@samples_metadata <- model@samples_metadata %>%
  left_join(metadata, by = "sample")

#3 Definitions
K <- get_dimensions(model)$K
factors_to_plot <- seq_len(min(5, K))

#4 Factor plots
p_factor_cellline <- plot_factor(model, factors = factors_to_plot, color_by = "CellLine",
                                 dot_size = 1) +
  scale_colour_okabe_ito() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(filename = file.path(outdir, "Factor_by_CellLine.svg"), plot = p_factor_cellline,
       device = svglite, width = 7, height = 5)

p_factor_modality <- plot_factor(model, factors = factors_to_plot, color_by = "Modality",
                                 dot_size = 1) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(filename = file.path(outdir, "Factor_by_Modality.svg"), plot = p_factor_modality,
       device = svglite, width = 7, height = 5)

#5 UMAP
set.seed(42)
factors_umap <- seq_len(min(3, K))
model <- run_umap(model, factors = factors_umap)

p_umap_cellline <- plot_dimred(model, method = "UMAP", color_by = "CellLine",
                               shape_by = "Modality", label = FALSE, dot_size = 2,
                               legend = TRUE) +
  scale_colour_okabe_ito()

ggsave(filename = file.path(outdir, "UMAP_CellLine_shape_Modality.svg"), plot = p_umap_cellline,
       device = svglite, width = 8, height = 6)

p_umap_modality <- plot_dimred(model, method = "UMAP", color_by = "Modality",
                               shape_by = "CellLine", label = FALSE, dot_size = 2,
                               legend = TRUE)

ggsave(filename = file.path(outdir, "UMAP_Modality_shape_CellLine.svg"), plot = p_umap_modality,
       device = svglite, width = 8, height = 6)