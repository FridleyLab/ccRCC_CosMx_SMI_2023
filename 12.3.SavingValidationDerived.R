
# clean slate -------------------------------------------------------------

rm(list=ls())

# libraries ---------------------------------------------------------------

library(tidyverse)
library(openxlsx)

# import data -------------------------------------------------------------

mif = readRDS('Manley_SMI/data/bivarK_mif/manley_mif_validation.rds')
clinical = read.csv("Manley_SMI/data/manley_files/Clinical.data_acs_02.14.23_final.csv", row.names = 1) %>%
  mutate(tissue = paste0("RCC", Slide)) %>%
  rename("fov" = "FOV") %>%
  filter(Histology == "clear cell")
key = read.xlsx("Manley_SMI/data/manley_files/FOV_TMA_key.xlsx", check.names = F)
summary_data = mif$sample

# get TMA names from image names ------------------------------------------

image_names = data.frame(`Image Tag` = mif$sample$`Image Tag`,
                         check.names = F)
image_names$tmp_key = image_names$`Image Tag` %>% gsub("\\].*", "", .) %>% gsub('.*\\-', "", .)
image_names$Slide = image_names$tmp_key %>% gsub("\\ .*", "", .) %>% paste0("RCC", .)
image_names$TMA.Row = image_names$tmp_key %>% gsub(".*\\[", "", .) %>% gsub("\\,.*", "", .) %>% as.numeric()
image_names$TMA.Col = image_names$tmp_key %>% gsub(".*\\,", "", .)


# merge key and tma sample names ------------------------------------------

key = key %>% 
  full_join(image_names) %>%
  mutate(tissue = Slide,
         fov = FOV) %>%
  full_join(clinical, by = c("tissue" = "tissue",
                             "fov" = "fov")) %>%
  full_join(summary_data)


# univar and bivar outputs ------------------------------------------------

univar = mif$derived$univariate_Count
bivar = mif$derived$bivariate_Count


# merge uni and bivar to clinical -----------------------------------------

uni_clin = key %>%
  full_join(univar)
bi_clin = key %>%
  full_join(bivar)


# SMI r=150 meaning -------------------------------------------------------
150 * 0.180
#27 um from r = 150px and 0.180um/px

#subset mIF to r=27um or r = 54px
uni_clin_27 = uni_clin %>%
  filter(r == 54)
bi_clin_27 = bi_clin %>%
  filter(r == 54)

# save validation data ----------------------------------------------------

saveRDS(uni_clin, "Manley_SMI/data/validation_dataframes/full_validation_clinical_univariate.rds")
saveRDS(bi_clin, "Manley_SMI/data/validation_dataframes/full_validation_clinical_bivariate.rds")
saveRDS(uni_clin_27, "Manley_SMI/data/validation_dataframes/validation_clinical_univariate_27um.rds")
saveRDS(bi_clin_27, "Manley_SMI/data/validation_dataframes/validation_clinical_bivariate_27um.rds")
