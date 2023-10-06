rm(list=ls())
#libraries
library(tidyverse)
library(spatstat)
library(pbmcapply)
library(spdep)
library(sfdep)
library(msigdbr)

#getting ligand pair
ligand_pairs = readRDS("Manley_SMI/data/ligand_receptor_celltalkdb/human_lr_pair.rds")
rm(list=ls())

#data
scaled_expression = readRDS("Manley_SMI/data/final_dataframes/sct_scaled_data.rds")
metadata = readRDS("Manley_SMI/data/final_dataframes/metadata_clinical_spatial.rds") %>%
  filter(fov <= 20, unique_fov != "RCC5_13") %>%
  mutate(Sarcomatoid = gsub(" ", "", Sarcomatoid),
         Site = gsub("\ .*", "", Site),
         Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received IO"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received IO")))

#filter ligand list to those in data
genes = colnames(scaled_expression)
ligand_pairs = readRDS("Manley_SMI/data/ligand_receptor_celltalkdb/human_lr_pair.rds")
ligand_pairs_available = ligand_pairs %>%
  filter(ligand_gene_symbol %in% genes, receptor_gene_symbol %in% genes)

#gene sets
all_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
msigdbr_list = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)

#which ligand-receptors are in EMT
EMT_genes = intersect(msigdbr_list$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, colnames(scaled_expression))
EMT_lr_pairs = ligand_pairs_available %>%
  filter(ligand_gene_symbol %in% EMT_genes, receptor_gene_symbol %in% EMT_genes)

#loop to run
fovs = unique(metadata$unique_fov)

morans_list = pbmclapply(setNames(fovs, fovs), function(fov){
  fov_meta = metadata %>% filter(unique_fov == !!fov)
  fov_expression = scaled_expression[fov_meta$id,]
  spdf = st_as_sf(fov_meta %>%
                    left_join(as.data.frame(fov_expression, check.names = F) %>%
                                rownames_to_column("id")),
                  coords = c("CenterX_local_px", "CenterY_local_px"))
  knn = knearneigh(st_coordinates(spdf), k = 3)
  knn.nb = knn2nb(knn)
  knn.nb.listw = nb2listw(knn.nb, style = "B")
  
  df = lapply(EMT_lr_pairs$lr_pair, function(pair){
    global_moran = moran_bv(spdf[[EMT_lr_pairs$ligand_gene_symbol[EMT_lr_pairs$lr_pair == pair]]], 
                            spdf[[EMT_lr_pairs$receptor_gene_symbol[EMT_lr_pairs$lr_pair == pair]]], knn.nb.listw, nsim = 100)
    data.frame(original = global_moran$t0,
               lr_pair = pair)
  }) %>% 
    do.call(bind_rows, .) %>%
    mutate(unique_fov = fov)
}, mc.cores = 5)
#saveRDS(morans_list, "Manley_SMI/results/ligand_receptor/Manley_ligand_receptor_CellTalkDB-knn3_EMT.rds")
morans_list = readRDS("Manley_SMI/results/ligand_receptor/Manley_ligand_receptor_CellTalkDB-knn3_EMT.rds")

#Making data frame of data
morans_df = do.call(bind_rows, morans_list) %>%
  spread("lr_pair", "original")
lr_pairs=EMT_lr_pairs$lr_pair

#Associations with primary tumor treatment tumor
primary_prepost_tumor_meta = metadata %>%
  filter(Site == "Tumor", Source == "Tumor") %>%
  filter(!(Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive")) %>%
  filter(!(unique_fov %in% c("RCC5_19", "RCC5_20"))) %>%
  select(unique_fov, Pretreatment.IO, Sarcomatoid) %>% distinct()

primary_prepost_tumor_res = lapply(setNames(lr_pairs, lr_pairs), function(lr_pair){
  df = primary_prepost_tumor_meta %>%
    left_join(morans_df %>% select(unique_fov, !!lr_pair), by = join_by(unique_fov))
  t.test(get(lr_pair) ~ Pretreatment.IO, data = df) %>%unlist() %>% t() %>% data.frame(check.names = F) %>% mutate(lr_pair = !!lr_pair)
}) %>%
  do.call(bind_rows, .) %>%
  mutate(p.adj = p.adjust(p.value, "fdr"), .after = p.value)

#Associations with primary tumor treatment stroma
primary_prepost_stroma_meta = metadata %>%
  filter(Site == "Tumor", Source == "Stroma") %>%
  filter(!(Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive")) %>%
  filter(!(unique_fov %in% c("RCC5_19", "RCC5_20"))) %>%
  select(unique_fov, Pretreatment.IO) %>% distinct()

primary_prepost_stroma_res = lapply(setNames(lr_pairs, lr_pairs), function(lr_pair){
  df = primary_prepost_stroma_meta %>%
    left_join(morans_df %>% select(unique_fov, !!lr_pair), by = join_by(unique_fov))
  t.test(get(lr_pair) ~ Pretreatment.IO, data = df) %>%unlist() %>% t() %>% data.frame(check.names = F) %>% mutate(lr_pair = !!lr_pair)
}) %>%
  do.call(bind_rows, .) %>%
  mutate(p.adj = p.adjust(p.value, "fdr"), .after = p.value)

#primary pre-io to sarcomatoid post-io tumor
sarcomatoid_primary_tumor_meta = metadata %>%
  filter(Site == "Tumor", Source == "Tumor") %>%
  select(unique_fov, Pretreatment.IO, Sarcomatoid) %>% distinct() %>%
  filter(Pretreatment.IO != "Received IO") %>% 
  mutate(Sarcomatoid.Class = ifelse(Pretreatment.IO == "Treatment Naive" & Sarcomatoid == "No",
                                    "Treatment Naive non-Sarcomatoid",
                                    "Sarcomatoid"))

sarcomatoid_primary_tumor_res = lapply(setNames(lr_pairs, lr_pairs), function(lr_pair){
  df = sarcomatoid_primary_tumor_meta %>%
    left_join(morans_df %>% select(unique_fov, !!lr_pair), by = join_by(unique_fov))
  t.test(get(lr_pair) ~ Sarcomatoid.Class, data = df) %>% unlist() %>% t() %>% data.frame(check.names = F) %>% mutate(lr_pair = !!lr_pair)
}) %>%
  do.call(bind_rows, .) %>%
  mutate(p.adj = p.adjust(p.value, "fdr"), .after = p.value)

#primary pre-io to sarcomatoid post-io stroma
sarcomatoid_primary_stroma_meta = metadata %>%
  filter(Site == "Tumor", Source == "Stroma") %>%
  select(unique_fov, Pretreatment.IO, Sarcomatoid) %>% distinct() %>%
  filter(Pretreatment.IO != "Received IO") %>% 
  mutate(Sarcomatoid.Class = ifelse(Pretreatment.IO == "Treatment Naive" & Sarcomatoid == "No",
                                    "Treatment Naive non-Sarcomatoid",
                                    "Sarcomatoid"))

sarcomatoid_primary_stroma_res = lapply(setNames(lr_pairs, lr_pairs), function(lr_pair){
  df = sarcomatoid_primary_stroma_meta %>%
    left_join(morans_df %>% select(unique_fov, !!lr_pair), by = join_by(unique_fov))
  t.test(get(lr_pair) ~ Sarcomatoid.Class, data = df) %>% unlist() %>% t() %>% data.frame(check.names = F) %>% mutate(lr_pair = !!lr_pair)
}) %>%
  do.call(bind_rows, .) %>%
  mutate(p.adj = p.adjust(p.value, "fdr"), .after = p.value)


# Save results of ligand - receptor to excel ------------------------------

library(openxlsx)

#create workbook
wb = createWorkbook()
addWorksheet(wb, "Primary Tumor Pre Post IO")
writeData(wb, sheet = 1, primary_prepost_tumor_res %>% arrange((p.value)))
addWorksheet(wb, "Primary Stroma Pre Post IO")
writeData(wb, sheet = 2, primary_prepost_stroma_res %>% arrange((p.value)))
addWorksheet(wb, "Sarco Tumor Pre Prim Sarco Post")
writeData(wb, sheet = 3, sarcomatoid_primary_tumor_res %>% arrange((p.value)))
addWorksheet(wb, "Sarco Strom Pre Prim Sarco Post")
writeData(wb, sheet = 4, sarcomatoid_primary_stroma_res %>% arrange((p.value)))
saveWorkbook(wb, "Manley_SMI/results/ligand_receptor/Primary_Sarcomatoid_ligand-receptor_bivarMoransI_knn3_EMT.xlsx", overwrite = TRUE)
