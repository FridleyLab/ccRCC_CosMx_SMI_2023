
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
  filter(Histology == "Papillary")
#extract cell level characteristics
meta_data = manley_data@meta.data
#create cell level clinical data
test = left_join(clinical, meta_data) %>%
  filter(!is.na(orig.ident)) %>%
  unite("unique_fov", c(tissue, fov), remove = F)

#find and import spatial information for cells
position_files = list.files("Manley_SMI/data/SMI-0050_BrandonManley_Moffit/5 Raw data/", recursive = T, pattern = "metadata_file.csv", full.names = T)
position_data = lapply(position_files, function(x){
  tissue = str_split(x, "/") %>% 
    unlist() %>% 
    grep("metadata", ., value = T) %>%
    str_split(., "_") %>% 
    unlist() %>% 
    grep("RCC", ., value = T)
  tmp = read.csv(x)
  tmp$tissue = tissue
  tmp %>%
    select(fov, Area, CenterX_local_px, CenterY_local_px,
           Width, Height, Mean.CD3, Mean.DAPI, Mean.PanCK) %>%
    return()
}) %>%
  do.call(bind_rows, .)
#join cell-clinical with position information, and collapse down cell assignments so simpler names
test = test %>%
  left_join(position_data)

sp = test %>% select(Specimen.ID.primary, unique_fov, 
                     CenterX_local_px, CenterY_local_px, lasso_final) %>% 
  mutate(pos = 1) %>% 
  pivot_wider(names_from = lasso_final, values_from = pos, values_fill = 0)

#split the cell data into FOV data frames and collapse in a list
spatial_list = lapply(unique(test$unique_fov), function(x){
  dat = sp[sp$unique_fov == x,]
  if(nrow(dat) < 5){
    return(NULL)
  }
  dat
})
#give each the unique ID
names(spatial_list) = unique(test$unique_fov)
#remove those that had less than 5 cells 
spatial_list = spatial_list[!sapply(spatial_list, is.null)]
saveRDS(test, "Manley_SMI/data/final_dataframes/metadata_clinical_spatial_papillary.rds")
scaled_expression = manley_data@assays$SCT@scale.data[,colnames(manley_data@assays$SCT@scale.data) %in% test$id] %>% t()
saveRDS(scaled_expression, "Manley_SMI/data/final_dataframes/sct_scaled_data_papillary.rds")
counts_expression = manley_data@assays$SCT@counts[,colnames(manley_data@assays$SCT@counts) %in% test$id] %>% SeuratDisk::Transpose()
saveRDS(counts_expression, "Manley_SMI/data/final_dataframes/sct_counts_data_papillary.rds")

# Create mIF Object -------------------------------------------------------

#1.3.3.3 feature-alex branch
library(spatialTIME)
#pull out the collapsed cell types
markers = test$lasso_final %>% unique() %>% sort()

#calculate fov cell counts and abundance
sample_summary = lapply(spatial_list, function(x){
  x %>% 
    group_by(Specimen.ID.primary, unique_fov) %>%
    summarise(across(any_of(markers %>% sort()), list(Cells = ~ sum(.x)), .names = "{.col} {.fn}"),
              `Total Cells` = n()) %>%
    mutate(across(contains(" Cells"), list(`%` = ~ .x / `Total Cells` * 100), .names = "{.col} {.fn}"))
}) %>%
  do.call(bind_rows, .)

#create manleys mif
manley_mif = create_mif(clinical = clinical %>%
                          unite("unique_fov", c(tissue, fov), remove = F),
                        sample_data = sample_summary,
                        spatial_list = spatial_list,
                        sample_id = "unique_fov",
                        patient_id = "Specimen.ID.primary")
#save raw mif
saveRDS(manley_mif, "Manley_SMI/data/manual_insitutype_raw_mif/manley_mif_papillary.rds")
