
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(spatialTIME)
library(data.table)
library(openxlsx)

# data --------------------------------------------------------------------

spatial_files = list.files("Manley_SMI/data/MCC20148_3 RCC Panel 2 July 2023/Manley_P2_MCC20148 3 RCC TMA_ObjectData/", pattern = "csv", full.names = T)
names(spatial_files) = spatial_files %>%
  gsub("\\/", "\\\\", .) %>%
  gsub(".*\\\\", "", .) %>%
  gsub("tif.*", "tif", .)
spatial_data = lapply(spatial_files, function(x){
  fread(x) %>%
    mutate(`Image Tag` = gsub(".*\\\\", "", `Image Location`))
})

summary_data = read.xlsx("Manley_SMI/data/MCC20148_3 RCC Panel 2 July 2023/Manley_P2_MCC20148 3 RCC Combined TMA.xlsx") %>%
  rename_with(~ gsub("\\.", " ", .x)) %>%
  mutate(`Dummy Variable` = `Image Tag`, .before = 1)

clinical_data = data.frame(`Image Tag` = names(spatial_files),
                           `Dummy Variable` = names(spatial_files),
                           check.names = FALSE)

# making mif --------------------------------------------------------------

mif = create_mif(clinical = clinical_data,
                 sample_data = summary_data,
                 spatial_list = spatial_data,
                 patient_id = "Dummy Variable",
                 sample_id = "Image Tag")

# Ripley's K --------------------------------------------------------------

markers = spatial_data[[1]] %>%
  colnames() %>%
  grep("Pos|\\+", ., value = T) %>%
  grep("Cyto|Nucl", ., value = T, invert = T)

mif2 = ripleys_k(mif = mif,
                 mnames = markers,
                 r_range = 0:500,
                 edge_correction = "translation",
                 permute = F,
                 workers = 5)

#save
saveRDS(mif2, "Manley_SMI/data/UnivariateK_mif/manley_mif_validation.rds")
