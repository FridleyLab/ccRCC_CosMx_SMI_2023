
# Creating mIF object for spatial clustering analysis ---------------------
rm(list=ls())
#libraries
library(Seurat)
library(tidyverse)

#import data
manley_data = readRDS("Manley_SMI/results/reports/modeling_tumor/manley_data_tumor_final_fixedFOV.rds")
clinical = read.csv("Manley_SMI/data/manley_files/Clinical.data_acs_02.14.23_final.csv") %>%
  mutate(tissue = paste0("RCC", Slide)) %>%
  rename("fov" = "FOV") %>%
  filter(Histology == "clear cell",
         fov <= 20)
#extract cell level characteristics
meta_data = manley_data@meta.data
#create cell level clinical data
test = left_join(clinical, meta_data) %>%
  filter(!is.na(orig.ident)) %>%
  unite("unique_fov", c(tissue, fov), remove = F) %>%
  filter(grepl("Tumor", Site))

Cell_by_fov = test %>%
  group_by(Slide_name, fov, Histology, Sarcomatoid, IT.Treatment.before.collection) %>%
  summarise(cells = n())

write.csv(Cell_by_fov,
          "Manley_SMI/results/figures/FOV_Cell_counts_PrimaryTumor.csv")

#plotting
waterfall_pl = test %>% group_by(tissue, fov) %>% summarise(cells = n()) %>% mutate(`FOV Name` = paste(tissue, fov)) %>%
  ggplot() +
  geom_bar(aes(x = reorder(`FOV Name`, desc(cells)), y = cells, fill = tissue), stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  labs(x = "FOV Name (Descending Cell Count)", y = "Cells on FOV", title = "Number of Cells on each Field of View")


pdf("Manley_SMI/results/figures/FOV_Cell_counts_PrimaryTumor.pdf", height = 5, width = 8)
waterfall_pl
dev.off()