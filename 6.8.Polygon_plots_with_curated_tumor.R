rm(list=ls())
library(Seurat)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)

#load in the final phenotyped seurat object after running Insitutype (Feb 7, 2023)
manley_data = readRDS("Manley_SMI/results/reports/modeling_tumor/manley_data_tumor_final_fixedFOV.rds")

#identify the files with the polygon information for all FOVs on all slides
polygons_files = list.files("Manley_SMI/data/SMI-0050_BrandonManley_Moffit/5 Raw data/", 
                            recursive = T,pattern = "polygons", full.names = T)
#now read in the polygon slides
polygon_data = lapply(polygons_files, read.csv, check.names = F)

#for merging need to add which slide (tissue) that the polygon data is for
polygon_data[[1]]$tissue = "RCC3"
polygon_data[[2]]$tissue = "RCC4"
polygon_data[[3]]$tissue = "RCC5"

#concatonate the polygon data into a single data frame and convert the fov and cellID to characters for merging with the metadata
polygon_data = polygon_data %>% 
  do.call(bind_rows, .) %>%
  mutate_at(c("cellID"), as.character)

#extract metadata with the phenotype information from the seurat object
meta_data = manley_data@meta.data 
#split the Nanostring "cell_ID" column into cellID that matches with the polygon information data
meta_data$cellID = meta_data$cell_ID %>%
  gsub('.*_', '', .)

#collapse down the phenotypes that were identified with the insitutype function based on brandons recommendations (data/final_celltypes_BJM.xlsx)
full_data = meta_data %>%
  left_join(polygon_data) 

phenotypes = full_data$lasso_final %>% unique() %>% sort()
full_data$lasso_final = factor(full_data$lasso_final,
                                        levels = phenotypes[c(1, 3:4, 8, 11:15, 17:22, 24, 29,
                                                              2, 5:7, 9:10, 16, 23, 25:28, 31, 30)])
immune_color_func = colorRampPalette(c("red", "blue"))
immune_colors = immune_color_func(17)
tissue_color_func = colorRampPalette(c("green", "orange"))
tissue_colors = tissue_color_func(13)
all_colors = c(immune_colors, tissue_colors, "black")
names(all_colors) = levels(full_data$lasso_final)

#now we can use the group_map to output a list of plots for all of the slide/fov combinations
plots = full_data %>%
  group_by(tissue, fov) %>%
  group_map(~{
    .x %>%
      ggplot() + 
      geom_polygon(aes(group = cellID, fill = lasso_final, x = x_local_px, y = y_local_px),
                   color = "black", size = 0.1) +
      theme_bw() +
      labs(title = paste0("Slide ", .y$tissue, " - FOV ", .y$fov)) +
      tune::coord_obs_pred() + 
      guides(fill = guide_legend(title = "InSituType Assignments\n(Manley Collapsed)",ncol = 1)) + 
      scale_fill_manual(values = all_colors)
  })

#save the plots into PDF form into the figures folder
pdf("Manley_SMI/results/figures/InSituType_Phenotypes_polygon_plots_final-phenotypes.pdf", width = 14, height = 10)
tmp = lapply(plots, print)
dev.off()



# Adding jpeg -------------------------------------------------------------

library(imager)

pictures = list.files("Manley_SMI/data/SMI-0050_BrandonManley_Moffit/5 Raw data/", pattern = ".jpg", recursive = T, full.names = T) %>% 
  grep("CellComposite", ., value = T)
pictures = pictures[-c(46, 60)]

plots2 = plots[-4]

pdf("Manley_SMI/results/figures/InSituType_IF_plots.pdf", width = 14, height = 8)
tmp = lapply(seq(pictures), function(p){
  im = grid::rasterGrob(load.image(pictures[p]))
  ggarrange(plots2[[p]], 
            ggarrange(ggplot() +
                        theme_bw() + 
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              panel.background = element_blank()), im, nrow = 2, heights = c(0.3, 1)),
            ncol = 2, widths = c(1.7, 1)) %>% 
    print()
})
dev.off()
