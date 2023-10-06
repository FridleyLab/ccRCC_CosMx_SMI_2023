rm(list=ls())
#libraries
library(tidyverse)
library(pbmcapply)

#data
scaled_expression = readRDS("Manley_SMI/data/final_dataframes/sct_scaled_data.rds")
metadata = readRDS("Manley_SMI/data/final_dataframes/metadata_clinical_spatial.rds") %>%
  filter(fov <= 20, unique_fov != "RCC5_13") %>%
  mutate(Sarcomatoid = gsub(" ", "", Sarcomatoid),
         Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received IO"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received IO")),
         Sarcomatoid.Group = ifelse(Sarcomatoid == "No", "non-Sarcomatoid", "Sarcomatoid"),
         Sarcomatoid.Group = factor(Sarcomatoid.Group, levels = c("non-Sarcomatoid", "Sarcomatoid"))) #sarcomatoid *compared to* non-sarcomatoid
non_sarc_tumor = metadata %>% 
  filter(Source == "Tumor", Site == "Tumor", Pretreatment.IO == "Treatment Naive")

GEX_out = lapply(colnames(scaled_expression), function(gene){
  gex_dat = scaled_expression[non_sarc_tumor$id,gene] %>% data.frame(check.names = F) %>% rename(gene := 1) %>% rownames_to_column("id") %>%
    inner_join(non_sarc_tumor %>% select(Sarcomatoid.Group, unique_fov, id), by = "id") %>%
    group_by(Sarcomatoid.Group, unique_fov) %>%
    summarise(gene = mean((gene)))
  res = t.test(gene ~ Sarcomatoid.Group, data = gex_dat)
  return(data.frame(Gene = gene,
                    P.value = res$p.value) %>%
           bind_cols(res$estimate %>% t() %>% data.frame(check.names = F)))
}) %>% 
  do.call(bind_rows, .) %>%
  mutate(P.adj = p.adjust(P.value, method = "fdr"), .after = P.value)

#saveRDS(GEX_out, "Manley_SMI/results/GEx/bulk_tumor_sarcomatoid.rds")

#Stroma DEGs
non_sarc_stroma = metadata %>% 
  filter(Source == "Stroma", Site == "Tumor", Pretreatment.IO == "Treatment Naive")

GEX_out = pbmclapply(colnames(scaled_expression), function(gene){
  gex_dat = scaled_expression[non_sarc_stroma$id,gene] %>% data.frame(check.names = F) %>% rename(gene := 1) %>% rownames_to_column("id") %>%
    inner_join(non_sarc_stroma %>% select(Sarcomatoid.Group, unique_fov, id), by = "id") %>%
    group_by(Sarcomatoid.Group, unique_fov) %>%
    summarise(gene = mean((gene)))
  res = t.test(gene ~ Sarcomatoid.Group, data = gex_dat)
  return(data.frame(Gene = gene,
                    P.value = res$p.value) %>%
           bind_cols(res$estimate %>% t() %>% data.frame(check.names = F)))
}, mc.cores = 5) %>% 
  do.call(bind_rows, .) %>%
  mutate(P.adj = p.adjust(P.value, method = "fdr"), .after = P.value)

#saveRDS(GEX_out, "Manley_SMI/results/GEx/bulk_stroma_sarcomatoid.rds")


# Saving CSV --------------------------------------------------------------


primary_tumor = readRDS("Manley_SMI/results/GEx/bulk_tumor_sarcomatoid.rds") %>% mutate(Source = "Tumor")
primary_stroma = readRDS("Manley_SMI/results/GEx/bulk_stroma_sarcomatoid.rds") %>% mutate(Source = "Stroma")

bulk_primary = primary_tumor %>%
  bind_rows(primary_stroma)
write.csv(bulk_primary, "Manley_SMI/results/GEx/bulk_sarcomatoid_DEGs.csv")
