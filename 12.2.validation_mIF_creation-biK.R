# Libraries ---------------------------------------------------------------

library(tidyverse)
library(spatialTIME)
library(data.table)
library(openxlsx)


# Univar mIF --------------------------------------------------------------

mif2 = readRDS("Manley_SMI/data/UnivariateK_mif/manley_mif_validation.rds")
#markers
markers = mif2$spatial[[1]] %>%
  colnames() %>%
  grep("Pos|\\+", ., value = T) %>%
  grep("Cyto|Nucl", ., value = T, invert = T)


# bivar ripley's k  -------------------------------------------------------

mif3 = bi_ripleys_k(mif = mif2, mnames = markers, r_range = 0:500, edge_correction = "translation", permute = F, workers = 32, big = 10000, nlarge = 10000)


# save to bivariate_mif folder --------------------------------------------

saveRDS(mif3, "Manley_SMI/data/bivarK_mif/manley_mif_validation.rds")


# Dixon's -----------------------------------------------------------------

mif4 = dixons_s(mif = mif3, mnames = markers, num_permutations = 1000, type = "Z", workers = 32)


saveRDS(mif4, "Manley_SMI/data/bivarK_mif/manley_mif_validation_Dixons.rds")
