library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(MOFAdata)
library(dplyr)
library(MOFA2)

#Parameters
top_genes_n  <- 3
max_pathways <- 15
max_genes    <- 8

outdir <- "PATH_TO_OUTPUT_DIRECTORY"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#Load model
setwd("PATH_TO_WORKING_DIRECTORY")
model <- load_model("PATH/TO/SCOTplus_trained_model.hdf5")

model@samples_metadata <- model@samples_metadata %>%
  mutate(cell_line = sub("^[^_]+_([^_]+).*", "\\1", sample))

data("MSigDB_v6.0_C5_mouse")

#Settings
views  <- c("Transcriptomics","Proteomics")
signs  <- c("positive","negative")
factors_all <- seq_len(get_dimensions(model)$K)

#Normalize feature names once only
for(v in views){
  suffix <- if(v=="Transcriptomics") "_RNA$" else "_PROT$"
  features_names(model)[[v]] <- toupper(features_names(model)[[v]])
  features_names(model)[[v]] <- sub(suffix, "", features_names(model)[[v]])
}

#Helper functions
empty_plot <- function(label){
  ggplot() +
    annotate("text", 0.5, 0.5, label=label, size=5,
             hjust=0.5, vjust=0.5) +
    theme_void()
}

has_sign_weights <- function(model, view, fac, sign, min_genes = 5){
  w <- get_weights(model, views=view, factors=fac, as.data.frame=TRUE)
  if(nrow(w)==0) return(FALSE)
  
  if(sign=="positive"){
    sum(w$value > 0) > min_genes
  } else {
    sum(w$value < 0) > min_genes
  }
}

#Main loop
for(view in views){
  
  message("Processing view: ", view)
  
  for(sign in signs){
    
    message("  Sign: ", sign)
    
    #Filter factors with usable weights
    valid_factors <- factors_all[
      sapply(factors_all, function(f)
        has_sign_weights(model, view, f, sign))
    ]
    
    if(length(valid_factors) == 0){
      message("    No valid factors → skipping")
      next
    }
    
    #Run enrichment for valid factors
    enrichment <- tryCatch({
      run_enrichment(
        model,
        view = view,
        factors = valid_factors,
        feature.sets = MSigDB_v6.0_C5_mouse,
        sign = sign,
        statistical.test = "parametric"
      )
    }, error=function(e){
      message("    Enrichment failed: ", e$message)
      return(NULL)
    })
    
    if(is.null(enrichment)) next
    
    #Heatmap per viwer and sign
    heat_file <- file.path(outdir, paste0(view,"_",sign,"_enrichment_heatmap.png"))
    
    sig_any <- any(sapply(valid_factors,
                          function(f) length(enrichment$sigPathways[[f]]) > 0))
    
    if(sig_any){
      p_heat <- plot_enrichment_heatmap(enrichment, factor = valid_factors)
    } else {
      p_heat <- empty_plot(paste0(view," | ",sign,"\nNo significant pathways"))
    }
    
    ggsave(heat_file, p_heat, width=8, height=6, dpi=300)
    
    #Per factor plots
    for(fac in valid_factors){
      
      message("    Factor ", fac)
      
      fac_dir <- file.path(outdir,
                           paste0(view,"_",sign,"/Factor",fac))
      dir.create(fac_dir, recursive=TRUE, showWarnings=FALSE)
      
      weights <- get_weights(model, views=view, factors=fac, as.data.frame=TRUE)
      sig_pathways <- enrichment$sigPathways[[fac]]
      
      #Top pathways
      p_top <- if(length(sig_pathways)==0){
        empty_plot("No significant pathways")
      } else {
        plot_enrichment(enrichment, factor=fac,
                        max.pathways=max_pathways, text_size = 0.7)
      }
      
      ggsave(file.path(fac_dir,"top_enrichment.png"),
             p_top, width=6, height=5, dpi=300)
      
      #Detailed enrichment
      p_det <- if(length(sig_pathways)==0){
        empty_plot("No significant pathways")
      } else {
        plot_enrichment_detailed(enrichment, factor=fac,
                                 max.genes=max_genes,
                                 max.pathways=5)
      }
      
      ggsave(file.path(fac_dir,"detailed_enrichment.png"),
             p_det, width=7, height=6, dpi=300)
      
    } # factor loop
  } # sign loop
} # view loop

#Mapping specific genes of interest UMAP
genes <- list("WNT5A","RNF2","MKNK1")

genes %>% map(~ plot_factors(model, 
                             factors = c(1,2), 
                             color_by = "cell_line", 
                             scale = T,
                             legend = T
)) %>% cowplot::plot_grid(plotlist=., nrow=1)

genes %>% map(~ plot_factors(model, 
                             factors = c(1,2), 
                             color_by = ., 
                             scale = T,
                             legend = T
)) %>% cowplot::plot_grid(plotlist=., nrow=1)