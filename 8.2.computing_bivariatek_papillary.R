rm(list=ls())
# Univariate RipK ---------------------------------------------------------

#libraries
library(spatialTIME)
library(tidyverse)

#read in datas
manley_mif = readRDS("Manley_SMI/data/UnivariateK_mif/manley_mif_papillary.rds")

#grab all of the markers
markers = colnames(manley_mif$sample) %>%
  grep("\\%", ., value = T) %>%
  grep("Total", ., value = T, invert = T) %>%
  gsub(" Cells \\%", "", .)

# system.time({
#   manley_mif2 =  bi_ripleys_k(mif = manley_mif,
#                               mnames = markers,
#                               r_range = seq(0, 500, 5),
#                               num_permutations = 10,
#                               edge_correction = "translation",
#                               permute = FALSE,
#                               keep_permutation_distribution = FALSE,
#                               workers = 6,
#                               overwrite = T,
#                               xloc = "CenterX_local_px",
#                               yloc = "CenterY_local_px")
# })
# 
# saveRDS(manley_mif2, "Manley_SMI/data/bivarK_mif/manley_mif_papillary.rds")


system.time({
  manley_mif2 =  bi_ripleys_k(mif = manley_mif,
                              mnames = markers[c(1,3,4,7,8,10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 29)],
                              r_range = seq(0, 500, 25),
                              num_permutations = 1000,
                              edge_correction = "translation",
                              permute = TRUE,
                              keep_permutation_distribution = FALSE,
                              workers = 32,
                              overwrite = T,
                              xloc = "CenterX_local_px",
                              yloc = "CenterY_local_px",
                              big = 5000,
                              nlarge = 5000)
})
saveRDS(manley_mif2, "Manley_SMI/data/bivarK_mif/manley_mif_papillary_1000perms.rds")

