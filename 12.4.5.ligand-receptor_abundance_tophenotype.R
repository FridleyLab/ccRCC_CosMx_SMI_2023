
# clean slate -------------------------------------------------------------

rm(list=ls())


# libraries ---------------------------------------------------------------

library(tidyverse)
library(VGAM)


# read in data ------------------------------------------------------------

uni_k_dat = readRDS("Manley_SMI/data/validation_dataframes/validation_clinical_univariate_27um.rds") %>%
  filter(!is.na(Source)) %>%
  mutate(Site = gsub("\\ $", "", Site),
         Sarcomatoid = gsub("\\ $", "", Sarcomatoid),
         Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received IO"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received IO"))) %>%
  mutate(unique_fov = paste0(tissue, "_FOV", fov)) %>%
  filter(Site == "Tumor",
         FOV <= 20,
         !(unique_fov %in% c("RCC5_FOV3", "RCC5_FOV13", "RCC5_FOV17", "RCC5_FOV19", "RCC5_FOV20")))
mif = readRDS("Manley_SMI/data/UnivariateK_mif/manley_mif_validation.rds")

# notes from brandon ------------------------------------------------------
#VIM and SMA in stroma = fibroblasts/epithelial
#VIM in tumor = fibroblast/malignant

markers = uni_k_dat$Marker %>%
  unique()


# calculating proportions -------------------------------------------------

sma_itgav_props = lapply(mif$spatial, function(spat){
  spat %>% 
    mutate(`Image Tag` = basename(gsub("\\\\", "\\/",`Image Location`))) %>%
    group_by(`Image Tag`, 
             `SMA (Opal 570) Positive Classification`,
             `ITGAV (Opal 520) Positive Classification`) %>%
    summarise(n())
}) %>%
  do.call(bind_rows, .) 
sma_col4a1_props = lapply(mif$spatial, function(spat){
  spat %>% 
    mutate(`Image Tag` = basename(gsub("\\\\", "\\/",`Image Location`))) %>%
    group_by(`Image Tag`, 
             `SMA (Opal 570) Positive Classification`,
             `COL4 (Opal 540) Positive Classification`) %>%
    summarise(n())
}) %>%
  do.call(bind_rows, .)
sma_col4a1_props = lapply(mif$spatial, function(spat){
  spat %>% 
    mutate(`Image Tag` = basename(gsub("\\\\", "\\/",`Image Location`))) %>%
    group_by(`Image Tag`, 
             `SMA (Opal 570) Positive Classification`,
             `COL4 (Opal 540) Positive Classification`) %>%
    summarise(n())
}) %>%
  do.call(bind_rows, .)
pck_col4a1_props = lapply(mif$spatial, function(spat){
  spat %>% 
    mutate(`Image Tag` = basename(gsub("\\\\", "\\/",`Image Location`))) %>%
    group_by(`Image Tag`, 
             `PCK (Opal 780) Positive Classification`,
             `COL4 (Opal 540) Positive Classification`) %>%
    summarise(n())
}) %>%
  do.call(bind_rows, .)

# calculating frequencies -------------------------------------------------

sma_itgav_props2 = sma_itgav_props %>%
  group_by(`Image Tag`) %>%
  mutate(freq = `n()` / sum(`n()`),
         total = sum(`n()`)) %>%
  unite("unique_phenotype", `SMA (Opal 570) Positive Classification`:`ITGAV (Opal 520) Positive Classification`, sep = "", remove = F)
sma_col4a1_props2 = sma_col4a1_props %>%
  group_by(`Image Tag`) %>%
  mutate(freq = `n()` / sum(`n()`),
         total = sum(`n()`)) %>%
  unite("unique_phenotype", `SMA (Opal 570) Positive Classification`:`COL4 (Opal 540) Positive Classification`, sep = "", remove = F)
pck_col4a1_props2 = pck_col4a1_props %>%
  group_by(`Image Tag`) %>%
  mutate(freq = `n()` / sum(`n()`),
         total = sum(`n()`)) %>%
  unite("unique_phenotype", `PCK (Opal 780) Positive Classification`:`COL4 (Opal 540) Positive Classification`, sep = "", remove = F)

# cohorts -----------------------------------------------------------------

pre_IO = uni_k_dat %>%
  filter(Sarcomatoid == "No", IT.Treatment.before.collection == "None") %>%
  select(unique_fov, Slide.x, FOV, Sarcomatoid, IT.Treatment.before.collection, Source, Pretreatment.IO, `Image Tag`) %>% 
  distinct()
post_IO = uni_k_dat %>%
  filter(IT.Treatment.before.collection != "None") %>%
  select(unique_fov, Slide.x, FOV, Sarcomatoid, IT.Treatment.before.collection, Source, Pretreatment.IO, `Image Tag`) %>% 
  distinct()
sarc = uni_k_dat %>%
  filter(Sarcomatoid == "Yes", IT.Treatment.before.collection == "None") %>%
  select(unique_fov, Slide.x, FOV, Sarcomatoid, IT.Treatment.before.collection, Source, `Image Tag`) %>% 
  distinct()


# beta_binomial -----------------------------------------------------------
#ITGAV
#pre-post stroma
wkn_dat = sma_itgav_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "11") %>% 
  filter(Source == "Stroma")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = coefficients(summary(model_fit)) %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column("Coef") %>% 
  filter(grepl("Pre", Coef)) %>%
  mutate(Cohort = "Pre-/Post-", Tissue = "Stroma", `Ligand/Receptor` = "ITGAV", `Cell Class` = "SMA")
#pre-post tumor
wkn_dat = sma_itgav_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "11") %>% 
  filter(Source == "Tumor")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Pre", Coef)) %>%
                      mutate(Cohort = "Pre-/Post-", Tissue = "Tumor", `Ligand/Receptor` = "ITGAV", `Cell Class` = "SMA"))
#pre-sarc stroma
wkn_dat = sma_itgav_props2 %>%
  right_join(bind_rows(pre_IO, sarc)) %>%
  filter(unique_phenotype == "11") %>% 
  filter(Source == "Stroma")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Sarcomatoid")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Sarc", Coef)) %>%
                      mutate(Cohort = "Pre-/Sarcomatoid", Tissue = "Stroma", `Ligand/Receptor` = "ITGAV", `Cell Class` = "SMA"))
#pre-post tumor
wkn_dat = sma_itgav_props2 %>%
  right_join(bind_rows(pre_IO, sarc)) %>%
  filter(unique_phenotype == "11") %>% 
  filter(Source == "Tumor")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Sarcomatoid")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Sarc", Coef)) %>%
                      mutate(Cohort = "Pre-/Sarcomatoid", Tissue = "Tumor", `Ligand/Receptor` = "ITGAV", `Cell Class` = "SMA"))

#COL4A1
#pre-post stroma
wkn_dat = sma_col4a1_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "11") %>% 
  filter(Source == "Stroma")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Pre", Coef)) %>%
                      mutate(Cohort = "Pre-/Post-", Tissue = "Stroma", `Ligand/Receptor` = "COL4A1", `Cell Class` = "SMA"))
#pre-post tumor
wkn_dat = sma_col4a1_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "11") %>% 
  filter(Source == "Tumor")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Pre", Coef)) %>%
                      mutate(Cohort = "Pre-/Post-", Tissue = "Tumor", `Ligand/Receptor` = "COL4A1", `Cell Class` = "SMA"))
#pre-sarc stroma
wkn_dat = sma_col4a1_props2 %>%
  right_join(bind_rows(pre_IO, sarc)) %>%
  filter(unique_phenotype == "11") %>% 
  filter(Source == "Stroma")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Sarcomatoid")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Sarc", Coef)) %>%
                      mutate(Cohort = "Pre-/Sarcomatoid", Tissue = "Stroma", `Ligand/Receptor` = "COL4A1", `Cell Class` = "SMA"))
#pre-post tumor
wkn_dat = sma_col4a1_props2 %>%
  right_join(bind_rows(pre_IO, sarc)) %>%
  filter(unique_phenotype == "11") %>% 
  filter(Source == "Tumor")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Sarcomatoid")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Sarc", Coef)) %>%
                      mutate(Cohort = "Pre-/Sarcomatoid", Tissue = "Tumor", `Ligand/Receptor` = "COL4A1", `Cell Class` = "SMA"))

#trying tumor cells
#pre-post stroma
wkn_dat = pck_col4a1_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "11") %>% 
  filter(Source == "Stroma")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Pre", Coef)) %>%
                      mutate(Cohort = "Pre-/Post-", Tissue = "Stroma", `Ligand/Receptor` = "COL4A1", `Cell Class` = "PCK"))
#pre-post tumor
wkn_dat = pck_col4a1_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "11") %>% 
  filter(Source == "Tumor")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Pre", Coef)) %>%
                      mutate(Cohort = "Pre-/Post-", Tissue = "Tumor", `Ligand/Receptor` = "COL4A1", `Cell Class` = "PCK"))

# save the results! -------------------------------------------------------

write.csv(results, "Manley_SMI/results/mIF_validation/beta-binomial_fibroblast-tumor_LR_results.csv")


# Fibroblast abundances ---------------------------------------------------

sma_props = lapply(mif$spatial, function(spat){
  spat %>% 
    mutate(`Image Tag` = basename(gsub("\\\\", "\\/",`Image Location`))) %>%
    group_by(`Image Tag`, 
             `SMA (Opal 570) Positive Classification`) %>%
    summarise(n())
}) %>%
  do.call(bind_rows, .) 
sma_props2 = sma_props %>%
  group_by(`Image Tag`) %>%
  mutate(freq = `n()` / sum(`n()`),
         total = sum(`n()`)) %>%
  unite("unique_phenotype", `SMA (Opal 570) Positive Classification`, sep = "", remove = F)
# beta_binomial -----------------------------------------------------------
#SMA
#pre-post stroma
wkn_dat = sma_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Stroma")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = coefficients(summary(model_fit)) %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column("Coef") %>% 
  filter(grepl("Pre", Coef)) %>%
  mutate(Cohort = "Pre-/Post-", Tissue = "Stroma", `Ligand/Receptor` = "None", `Cell Class` = "SMA")
#pre-post tumor
wkn_dat = sma_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Tumor")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Pre", Coef)) %>%
                      mutate(Cohort = "Pre-/Post-", Tissue = "Tumor", `Ligand/Receptor` = "None", `Cell Class` = "SMA"))

#pre-sarc stroma
wkn_dat = sma_props2 %>%
  right_join(bind_rows(pre_IO, sarc)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Stroma")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Sarcomatoid")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Sarc", Coef)) %>%
                      mutate(Cohort = "Pre-/Sarcomatoid", Tissue = "Stroma", `Ligand/Receptor` = "None", `Cell Class` = "SMA"))
#pre-sarc tumor
wkn_dat = sma_props2 %>%
  right_join(bind_rows(pre_IO, sarc)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Tumor")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Sarcomatoid")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Sarc", Coef)) %>%
                      mutate(Cohort = "Pre-/Sarcomatoid", Tissue = "Tumor", `Ligand/Receptor` = "None", `Cell Class` = "SMA"))


write.csv(results, "Manley_SMI/results/mIF_validation/beta-binomial_SMA_results.csv")
