rm(list=ls())
#libraries
library(tidyverse)
library(pbmcapply)

`%notin%` = Negate(`%in%`)

#data
scaled_expression = readRDS("Manley_SMI/data/final_dataframes/sct_scaled_data.rds")
metadata = readRDS("Manley_SMI/data/final_dataframes/metadata_clinical_spatial.rds") %>%
  filter(fov <= 20, unique_fov != "RCC5_13") %>%
  filter(unique_fov %notin% c("RCC5_19", "RCC5_20")) %>%
  mutate(Sarcomatoid = gsub(" ", "", Sarcomatoid),
         Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received IO"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received IO")))
tumor_pre_post = metadata %>% 
  filter(Source == "Tumor", Site == "Tumor") %>%
  filter(!(Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive"))

GEX_out = lapply(colnames(scaled_expression), function(gene){
  gex_dat = scaled_expression[tumor_pre_post$id,gene] %>% data.frame(check.names = F) %>% rename(gene := 1) %>% rownames_to_column("id") %>%
    inner_join(tumor_pre_post %>% select(Pretreatment.IO, unique_fov, id), by = "id") %>%
    group_by(Pretreatment.IO, unique_fov) %>%
    summarise(gene = mean((gene)))
  res = t.test(gene ~ Pretreatment.IO, data = gex_dat)
  return(data.frame(Gene = gene,
                    P.value = res$p.value) %>%
           bind_cols(res$estimate %>% t() %>% data.frame(check.names = F)))
}) %>% 
  do.call(bind_rows, .) %>%
  mutate(P.adj = p.adjust(P.value, method = "fdr"), .after = P.value)

#saveRDS(GEX_out, "Manley_SMI/results/GEx/bulk_tumor_pre-postIO.rds")

#Stroma DEGs
stroma_pre_post = metadata %>% 
  filter(Source == "Stroma", Site == "Tumor") %>%
  filter(!(Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive"))

GEX_out = pbmclapply(colnames(scaled_expression), function(gene){
  gex_dat = scaled_expression[stroma_pre_post$id,gene] %>% data.frame(check.names = F) %>% rename(gene := 1) %>% rownames_to_column("id") %>%
    inner_join(stroma_pre_post %>% select(Pretreatment.IO, unique_fov, id), by = "id") %>%
    group_by(Pretreatment.IO, unique_fov) %>%
    summarise(gene = mean((gene)))
  res = t.test(gene ~ Pretreatment.IO, data = gex_dat)
  return(data.frame(Gene = gene,
                    P.value = res$p.value) %>%
           bind_cols(res$estimate %>% t() %>% data.frame(check.names = F)))
}, mc.cores = 5) %>% 
  do.call(bind_rows, .) %>%
  mutate(P.adj = p.adjust(P.value, method = "fdr"), .after = P.value)

#saveRDS(GEX_out, "Manley_SMI/results/GEx/bulk_stroma_pre-postIO.rds")


# Saving CSV --------------------------------------------------------------


primary_tumor = readRDS("Manley_SMI/results/GEx/bulk_tumor_pre-postIO.rds") %>% mutate(Source = "Tumor")
primary_stroma = readRDS("Manley_SMI/results/GEx/bulk_stroma_pre-postIO.rds") %>% mutate(Source = "Stroma")

bulk_primary = primary_tumor %>%
  bind_rows(primary_stroma)
write.csv(bulk_primary, "Manley_SMI/results/GEx/bulk_primary_DEGs.csv")
