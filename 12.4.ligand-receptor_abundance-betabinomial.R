
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
         !(unique_fov %in% c("RCC5_FOV3", "RCC5_FOV13", "RCC4_FOV17", "RCC5_FOV17", "RCC5_FOV19", "RCC5_FOV20"))) #17 and 3 are Pre-IO
mif = readRDS("Manley_SMI/data/UnivariateK_mif/manley_mif_validation.rds")

# Markers -----------------------------------------------------------------

markers = uni_k_dat$Marker %>%
  unique()

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


# ITGAV and COL4A1 Abundances ---------------------------------------------

itgav_props = lapply(mif$spatial, function(spat){
  spat %>% 
    mutate(`Image Tag` = basename(gsub("\\\\", "\\/",`Image Location`))) %>%
    group_by(`Image Tag`, 
             `ITGAV (Opal 520) Positive Classification`) %>%
    summarise(n())
}) %>%
  do.call(bind_rows, .) 
itgav_props2 = itgav_props %>%
  group_by(`Image Tag`) %>%
  mutate(freq = `n()` / sum(`n()`),
         total = sum(`n()`)) %>%
  unite("unique_phenotype", `ITGAV (Opal 520) Positive Classification`, sep = "", remove = F)


col4a1_props = lapply(mif$spatial, function(spat){
  spat %>% 
    mutate(`Image Tag` = basename(gsub("\\\\", "\\/",`Image Location`))) %>%
    group_by(`Image Tag`, 
             `COL4 (Opal 540) Positive Classification`) %>%
    summarise(n())
}) %>%
  do.call(bind_rows, .) 
col4a1_props2 = col4a1_props %>%
  group_by(`Image Tag`) %>%
  mutate(freq = `n()` / sum(`n()`),
         total = sum(`n()`)) %>%
  unite("unique_phenotype", `COL4 (Opal 540) Positive Classification`, sep = "", remove = F)

# beta_binomial -----------------------------------------------------------
#ITGAV
#pre-post stroma
wkn_dat = itgav_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Stroma")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = coefficients(summary(model_fit)) %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column("Coef") %>% 
  filter(grepl("Pre", Coef)) %>%
  mutate(Cohort = "Pre-/Post-", Tissue = "Stroma", `Ligand/Receptor` = "ITGAV", `Cell Class` = "None")
#pre-post tumor
wkn_dat = itgav_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Tumor")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Pre", Coef)) %>%
                      mutate(Cohort = "Pre-/Post-", Tissue = "Tumor", `Ligand/Receptor` = "ITGAV", `Cell Class` = "None"))
#COL4A1
#pre-post stroma
wkn_dat = col4a1_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Stroma")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Pre", Coef)) %>%
                      mutate(Cohort = "Pre-/Post-", Tissue = "Stroma", `Ligand/Receptor` = "COL4A1", `Cell Class` = "None"))
#pre-post stroma
wkn_dat = col4a1_props2 %>%
  right_join(bind_rows(pre_IO, post_IO)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Tumor")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Pretreatment.IO")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Pre", Coef)) %>%
                      mutate(Cohort = "Pre-/Post-", Tissue = "Tumor", `Ligand/Receptor` = "COL4A1", `Cell Class` = "None"))

#ITGAV
#pre-sarc stroma
wkn_dat = itgav_props2 %>%
  right_join(bind_rows(pre_IO, sarc)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Stroma")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Sarcomatoid")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Sarc", Coef)) %>%
                      mutate(Cohort = "Pre-/Sarcomatoid", Tissue = "Stroma", `Ligand/Receptor` = "ITGAV", `Cell Class` = "None"))
#pre-sarc tumor
wkn_dat = itgav_props2 %>%
  right_join(bind_rows(pre_IO, sarc)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Tumor")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Sarcomatoid")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Sarc", Coef)) %>%
                      mutate(Cohort = "Pre-/Sarcomatoid", Tissue = "Tumor", `Ligand/Receptor` = "ITGAV", `Cell Class` = "None"))
#COL4A1
#pre-sarc stroma
wkn_dat = col4a1_props2 %>%
  right_join(bind_rows(pre_IO, sarc)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Stroma")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Sarcomatoid")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Sarc", Coef)) %>%
                      mutate(Cohort = "Pre-/Sarcomatoid", Tissue = "Stroma", `Ligand/Receptor` = "COL4A1", `Cell Class` = "None"))
#pre-sarc stroma
wkn_dat = col4a1_props2 %>%
  right_join(bind_rows(pre_IO, sarc)) %>%
  filter(unique_phenotype == "1") %>% 
  filter(Source == "Tumor")
fo = as.formula("cbind(wkn_dat[['n()']], wkn_dat[['total']] - wkn_dat[['n()']]) ~ Sarcomatoid")
model_fit = vglm(fo, betabinomial(zero = 2), data = wkn_dat)
results = bind_rows(results, coefficients(summary(model_fit)) %>% 
                      data.frame(check.names = F) %>% 
                      rownames_to_column("Coef") %>% 
                      filter(grepl("Sarc", Coef)) %>%
                      mutate(Cohort = "Pre-/Sarcomatoid", Tissue = "Tumor", `Ligand/Receptor` = "COL4A1", `Cell Class` = "None"))

results = results %>%
  mutate(FDR = p.adjust(`Pr(>|z|)`), .before = Cohort)
write.csv(results, "Manley_SMI/results/mIF_validation/beta-binomial_ITGAV-COL4A1_results.csv")


# plot abundance of primary pre-/post-IO ----------------------------------

summary_dat = mif$sample %>% 
  select(`Dummy Variable`, `ITGAV (Opal 520) Positive Cells`, `COL4 (Opal 540) Positive Cells`, `Total Cells`) %>%
  rename(`ITGAV+` = 2, `COL4+` = 3) %>%
  mutate(across(c(`ITGAV+`, `COL4+`), ~ .x / `Total Cells` * 100, .names = "{col} %")) %>%
  select(`Dummy Variable` , contains("%")) %>%
  gather("Marker", "Abundance", -`Dummy Variable`)

plotting_data = uni_k_dat %>%
  select(unique_fov, Pretreatment.IO, Source, `Dummy Variable`, Sarcomatoid) %>% 
  distinct() %>% 
  left_join(summary_dat)

prepost_pl = plotting_data %>% 
  filter(!(Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive")) %>%
  ggplot() +
  geom_boxplot(aes(x = Pretreatment.IO, y = Abundance, color = Marker)) +
  facet_grid(~Source) +
  theme_bw() +
  labs(title = "Abundance in IO Naive/Exposed") +
  theme(plot.title = element_text(hjust = 0.5))
presarc_pl = plotting_data %>% 
  filter(Pretreatment.IO != "Received IO") %>%
  ggplot() +
  geom_boxplot(aes(x = Sarcomatoid, y = Abundance, color = Marker)) +
  facet_grid(~Source) +
  theme_bw() +
  labs(title = "Abundance in IO Naive/Sarcomatoid") +
  theme(plot.title = element_text(hjust = 0.5))

pdf("Manley_SMI/results/mIF_validation/COL4A1_ITGAV_abundance_boxplots.pdf", height = 5, width = 7)
prepost_pl
presarc_pl
dev.off()