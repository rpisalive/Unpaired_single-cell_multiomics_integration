library(ggplot2)
library(dplyr)
library(MOFA2)
library(svglite)
library(ggokabeito)
library(extrafont)

setwd("PATH_TO_WORKING_DIRECTORY")
model <- load_model("PATH/TO/single_cell_trained_model.hdf5")
outdir <- "PATH_TO_OUTPUT_DIRECTORY"

loadfonts(device = "win")
theme_set(theme_bw(base_size = 16, base_family = "Arial"))

plot_data_overview(model)
model@samples_metadata
model@samples_metadata <- model@samples_metadata %>% 
  mutate(cell_line = sub("^[^_]+_([^_]+).*", "\\1", sample))
model@samples_metadata

# Total variance explained per view and group
model@cache$variance_explained$r2_total
model@cache$variance_explained$r2_per_factor
p_var <- plot_variance_explained(model, x="view", y="factor", plot_total = FALSE)
ggsave(filename = file.path(outdir, "Figure2A.svg"),plot = p_var,
       device = svglite, width = 8,height = 6)
plot_variance_explained(model, x="view", y="factor", plot_total = TRUE)

p_factor <- plot_factor(model, factor = 1:6, color_by = "cell_line")
p_factor <- plot_factor(model, factor = 1:6, color_by = "cell_line", dot_size = 1) +
  theme_bw(base_size = 12) +
  scale_colour_okabe_ito() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title   = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 12))

ggsave(filename = file.path(outdir, "Figure2B.svg"),
       plot = p_factor,
       device = svglite,
       width = 6,height = 4)

p <- plot_factor(model, factors = c(1,2,3,4,5,6),
                 color_by = "cell_line",
                 dot_size = 3, dodge = T,
                 legend = F, add_violin = T,
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

for (v in views) {
  for (f in factors) {
    
    message("Processing: ", v, " | Factor ", f)
    
    #Weights plot
    p_weights <- plot_weights(model,view = v,factor = f,nfeatures = 10,scale = TRUE,
                              abs = FALSE,text_size = 2)
    
    ggsave(file.path(outdir, paste0(v, "_F", f, "_weights.svg")),
           p_weights, width = 6, height = 4, device = svglite)
    
    #Top weights
    p_top <- plot_top_weights(model,view = v,factor = f,nfeatures = 20)
    
    ggsave(file.path(outdir, paste0(v, "_F", f, "_topweights.svg")),
      p_top, width = 6, height = 4, device = svglite)
    
    #Scatter
    p_scatter <- plot_data_scatter(model,view = v,factor = f,features = 12,
                                   add_lm = TRUE,color_by = "cell_line")
    
    ggsave(file.path(outdir, paste0(v, "_F", f, "_scatter.svg")),
      p_scatter, width = 10, height = 8, device = svglite)
    
    #Heatmap
    png(file.path(outdir, paste0(v, "_F", f, "_heatmap.png")),width = 1800,
        height = 1200,res = 300)
    
    plot_data_heatmap(model,view = v,factor = f,features = 20,cluster_rows = TRUE,
                      cluster_cols = TRUE,show_rownames = TRUE, show_colnames = TRUE,
                      fontsize_col = 4)
    dev.off()
  }
}

#Run UMAP
set.seed(42)
model <- run_umap(model)#, factors = c(1,2)) - Adjust the factors to include

p_umap <- plot_dimred(model,method = "UMAP",color_by = "cell_line",
                      color_name = "Cell Line",shape_by = NULL,label = FALSE,
                      dot_size = 2,legend = TRUE)

p_umap <- p_umap +
  scale_colour_okabe_ito() +
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),legend.text = element_text(size = 16))

ggsave(
  filename = file.path(outdir, "UMAP_cellline.svg"),
  plot = p_umap,
  device = svglite,
  width = 8,
  height = 6
)