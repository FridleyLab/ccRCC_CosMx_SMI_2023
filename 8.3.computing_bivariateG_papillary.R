library(spatialTIME)
library(tidyverse)
mif = readRDS("Manley_SMI/data/bivarK_mif/manley_mif_papillary_1000perms.rds")

#look at bik first before digging too far
bik = mif$derived$bivariate_Count
bik %>%
  filter(Anchor == "CD8 T cell", Counted == "M2 macrophage (CD163)") %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Permutation`, group = unique_fov)) +
  geom_line(data = . %>% filter(unique_fov %in% c("RCC3_19", "RCC3_17")),
            aes(x = r, y = `Degree of Clustering Permutation`, color = unique_fov))

mif$spatial = lapply(mif$spatial, function(dat){
  dat$xloc = dat$CenterX_local_px
  dat$yloc = dat$CenterY_local_px
  return(dat)
})

mif = bi_NN_G(mif = mif, 
              mnames = c("CD8 T cell", "M2 macrophage (CD163)"),
              r_range = 0:1000,
              num_permutations = 1000,
              edge_correction = "rs",
              keep_perm_dis = FALSE,
              workers = 6,
              xloc = "CenterX_local_px",
              yloc = "CenterY_local_px",
              overwrite = TRUE)


big = mif$derived$bivariate_NN
big %>%
  filter(Anchor == "CD8 T cell", Counted == "M2 macrophage (CD163)") %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Permutation`, group = unique_fov)) +
  geom_line(data = . %>% filter(!(unique_fov %in% c("RCC4_19", "RCC4_20"))),
            aes(x = r, y = `Degree of Clustering Permutation`, color = unique_fov))
