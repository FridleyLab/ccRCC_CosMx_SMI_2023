rm(list=ls())
library(tidyverse)

metadata = readRDS("Manley_SMI/data/final_dataframes/metadata_clinical_spatial.rds") %>%
  filter(fov <= 20, unique_fov != "RCC5_13") %>%
  filter(Site == "Tumor") %>%
  mutate(Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received IO"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received IO"))) %>%
  mutate(Sarcomatoid = gsub(" $", "", Sarcomatoid))
#onnly includes primary tumor enrichment scores
all_fov = readRDS("Manley_SMI/results/GEx/hallmark_pathway_list.rds")

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

test = full_join(polygon_data, metadata %>%
                   mutate(cellID = gsub('.*_', '', cell_ID))
)

all_enrichments = lapply(all_fov, function(x){
  x$Cell_Sets %>% t() %>% data.frame(check.names = F)
}) %>% 
  do.call(bind_rows, .) %>%
  rownames_to_column("id")

## STenrich results
res = readRDS("Manley_SMI/data/pathway_enrichment/STEnrich_output_gsva.rds")
res = lapply(names(res), function(nam){
  res[[nam]] %>% mutate(unique_fov = nam)
}) %>% do.call(bind_rows, .)
##

allo_example1 = test %>%
  filter(unique_fov == "RCC5_1") %>%
  left_join(all_enrichments) %>%
  ggplot() + 
  geom_polygon(aes(group = cellID, fill = HALLMARK_ALLOGRAFT_REJECTION, x = x_local_px, y = y_local_px),
               color = "black", size = 0.1) +
  theme_bw() +
  labs(title = "RCC5-1, Significant Allograft Rejection Clustering") +
  tune::coord_obs_pred() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2()
allo_example2 = test %>%
  filter(unique_fov == "RCC4_5") %>%
  left_join(all_enrichments) %>%
  ggplot() + 
  geom_polygon(aes(group = cellID, fill = HALLMARK_ALLOGRAFT_REJECTION, x = x_local_px, y = y_local_px),
               color = "black", size = 0.1) +
  theme_bw() +
  labs(title = "RCC4_5, Significant Allograft Rejection Clustering") +
  tune::coord_obs_pred() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2()

pdf("Manley_SMI/results/figures/Enrichment/allograft_rejection_example_RCC5-1.pdf", width = 10, height = 7)
allo_example
dev.off()

pdf("Manley_SMI/results/figures/Enrichment/allograft_rejection_example_RCC4_5.pdf", width = 10, height = 7)
allo_example2
dev.off()

hypo_example1 = test %>%
  filter(unique_fov == "RCC5_5") %>%
  left_join(all_enrichments) %>%
  ggplot() + 
  geom_polygon(aes(group = cellID, fill = HALLMARK_HYPOXIA, x = x_local_px, y = y_local_px),
               color = "black", size = 0.1) +
  theme_bw() +
  labs(title = "RCC5_5, Significant Hypoxia Clustering") +
  tune::coord_obs_pred() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2()
hypo_example2 = test %>%
  filter(unique_fov == "RCC5_1") %>%
  left_join(all_enrichments) %>%
  ggplot() + 
  geom_polygon(aes(group = cellID, fill = HALLMARK_HYPOXIA, x = x_local_px, y = y_local_px),
               color = "black", size = 0.1) +
  theme_bw() +
  labs(title = "RCC5_1, Significant Hypoxia Clustering") +
  tune::coord_obs_pred() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2()

pdf("Manley_SMI/results/figures/Enrichment/hypoxia_example_RCC5_5.pdf", width = 10, height = 7)
hypo_example1
dev.off()
pdf("Manley_SMI/results/figures/Enrichment/hypoxia_example_RCC5_1.pdf", width = 10, height = 7)
hypo_example2
dev.off()


# EMT ---------------------------------------------------------------------
#clustered
#RCC4_7
#RCC5_1
#RCC5_19
emt_1 = test %>%
  filter(unique_fov == "RCC5_1") %>%
  left_join(all_enrichments) %>%
  ggplot() + 
  geom_polygon(aes(group = cellID, fill = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, x = x_local_px, y = y_local_px),
               color = "black", size = 0.1) +
  theme_bw() +
  labs(title = "RCC5_1, Significant EMT Clustering") +
  tune::coord_obs_pred() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradientn(colors = c("#191970FF", "#0000FFFF", "#FFFFFFFF", "#FF0000FF", "#8B0000FF"))

pdf("Manley_SMI/results/figures/Enrichment/EMT_example_RCC5_1.pdf", width = 10, height = 7)
emt_1
dev.off()

#RCC3_3
#RCC3_5
#RCC3_9
emt_2 = test %>%
  filter(unique_fov == "RCC3_3") %>%
  left_join(all_enrichments) %>%
  ggplot() + 
  geom_polygon(aes(group = cellID, fill = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, x = x_local_px, y = y_local_px),
               color = "black", size = 0.1) +
  theme_bw() +
  labs(title = "RCC3_3, non-Significant EMT Clustering") +
  tune::coord_obs_pred() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradientn(colors = c("#191970FF", "#0000FFFF", "#FFFFFFFF", "#FF0000FF", "#8B0000FF"))

pdf("Manley_SMI/results/figures/Enrichment/EMT_example_RCC3_3.pdf", width = 10, height = 7)
emt_2
dev.off()
