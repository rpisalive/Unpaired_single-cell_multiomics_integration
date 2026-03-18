library(ggplot2)
library(dplyr)
library(MOFA2)

setwd("PATH_TO_WORKING_DIRECTORY")
model <- load_model("PATH/TO/SCOTplus_trained_model.hdf5")

plot_data_overview(model)
model@samples_metadata
model@samples_metadata <- model@samples_metadata %>%
  mutate(cell_line = sub("^[^_]+_([^_]+).*", "\\1", sample))
model@samples_metadata

#Total variance explained per view and group
model@cache$variance_explained$r2_total
model@cache$variance_explained$r2_per_factor
plot_variance_explained(model, x="view", y="factor", plot_total = FALSE)
plot_variance_explained(model, x="view", y="factor", plot_total = TRUE)

plot_factor(model, factor = 1:6, color_by = "cell_line")

p <- plot_factor(model, 
                 factors = c(1,2,3,4,5,6),
                 color_by = "cell_line",
                 dot_size = 3,
                 dodge = T,
                 legend = F,
                 add_violin = T,
                 violin_alpha = 0.25)

p <- p + 
  scale_color_manual(values=c("C10"="black", "SVEC"="red")) +
  scale_fill_manual(values=c("C10"="black", "SVEC"="red"))

print(p)

plot_factors(model, 
             factors = 1:6,
             color_by = "cell_line")
views <- c("Transcriptomics", "Proteomics")
K <- get_dimensions(model)$K
factors <- seq_len(K)

outdir <- "PATH_TO_OUTPUT_DIRECTORY"

for (v in views) {
  for (f in factors) {
    
    message("Processing: ", v, " | Factor ", f)
    
    #Weights
    p_weights <- plot_weights(
      model,
      view = v,
      factor = f,
      nfeatures = 10,
      scale = TRUE,
      abs = FALSE,
      text_size = 2)
    
    ggsave(file.path(outdir, paste0(v, "_F", f, "_weights.png")),
           p_weights, width = 6, height = 4, dpi = 300)
    
    
    #Top weights
    p_top <- plot_top_weights(
      model,
      view = v,
      factor = f,
      nfeatures = 20)
    
    ggsave(file.path(outdir, paste0(v, "_F", f, "_topweights.png")),
           p_top, width = 6, height = 4, dpi = 300)
    
    
    #Scatter
    p_scatter <- plot_data_scatter(model, view = v, factor = f, features = 12, 
                                   add_lm = TRUE, color_by = "cell_line")
    
    ggsave(file.path(outdir, paste0(v, "_F", f, "_scatter.png")),
           p_scatter, width = 10, height = 8, dpi = 300)
    
    
    #Heatmap
    png(file.path(outdir, paste0(v, "_F", f, "_heatmap.png")),
        width = 1800, height = 1200, res = 300)
    
    plot_data_heatmap(model, view = v, factor = f, features = 15, cluster_rows = TRUE,
                      cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
                      fontsize_col = 4)
    
    dev.off()
  }
}

#UMAP
set.seed(42)
model <- run_umap(model)
p_umap <- plot_dimred(model,  method = "UMAP",  color_by = "cell_line")

#Save to file
ggsave(filename = file.path(outdir, "UMAP_cellline.png"), plot = p_umap, width = 6,
       height = 5, dpi = 300)