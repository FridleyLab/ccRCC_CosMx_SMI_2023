rm(list=ls())
# Univariate RipK ---------------------------------------------------------

#libraries
library(spatialTIME)
library(tidyverse)

#read in datas
manley_mif = readRDS("Manley_SMI/data/manual_insitutype_raw_mif/manley_mif.rds")

#grab all of the markers
markers = colnames(manley_mif$sample) %>%
  grep("\\%", ., value = T) %>%
  grep("Total", ., value = T, invert = T) %>%
  gsub(" Cells \\%", "", .)

system.time({
  manley_mif2 =  ripleys_k(mif = manley_mif,
                          mnames = markers,
                          r_range = seq(0, 500, 5),
                          num_permutations = 10,
                          edge_correction = "translation",
                          method = "K",
                          permute = F,
                          keep_permutation_distribution = F,
                          workers = 32,
                          overwrite = T,
                          xloc = "CenterX_local_px",
                          yloc = "CenterY_local_px")
})

saveRDS(manley_mif2, "Manley_SMI/data/UnivariateK_mif/manley_mif.rds")
