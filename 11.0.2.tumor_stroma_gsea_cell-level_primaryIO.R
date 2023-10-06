rm(list=ls())
# libraries ---------------------------------------------------------------

library(pbmcapply)
library(msigdbr)
library(GSVA)
library(ComplexHeatmap)
library(tidyverse)

`%notin%` = Negate(`%in%`)

#data
scaled_expression = readRDS("Manley_SMI/data/final_dataframes/sct_scaled_data.rds")
metadata = readRDS("Manley_SMI/data/final_dataframes/metadata_clinical_spatial.rds") %>%
  filter(fov <= 20, unique_fov != "RCC5_13") %>%
  mutate(Site = gsub(" $", "", Site)) %>%
  filter(Site == "Tumor") %>%
  mutate(Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received IO"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received IO"))) 
#cell types
phenotypes = metadata$lasso_final %>% unique()
#get the stroma and tumor groups
tumor = metadata %>% filter(Source == "Tumor")
stroma = metadata %>% filter(Source == "Stroma")


# Gene Sets ---------------------------------------------------------------

all_gene_sets = msigdbr(species = "Homo sapiens") %>%
  filter(gs_cat == "H") %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
gene_set_list = split(all_gene_sets$gene_symbol, all_gene_sets$gs_name)

fovs = metadata$unique_fov %>% unique()
all_fov = pbmclapply(setNames(fovs, fovs), function(fov_id){
  fov_meta = metadata %>% filter(unique_fov == fov_id)
  fov_genes = scaled_expression[row.names(scaled_expression) %in% fov_meta$id,]
  
  
  cell_gs = gsva(t(fov_genes), gene_set_list)
  phen = unique(fov_meta$lasso_final)
  phen = setNames(phen, phen)
  av_gs = lapply(phen, function(x){
    rowMeans(cell_gs[,fov_meta$id[fov_meta$lasso_final == x]] %>% data.frame())
  }) %>% do.call(bind_cols, .)
  row.names(av_gs) = row.names(cell_gs) %>%
    gsub("HALLMARK_", "", .) %>%
    gsub("_", " ", .)
  return(list(`Cell_Sets` = cell_gs,
              Phenotype_sets = av_gs))
}, mc.cores = 40)
#saveRDS(all_fov, "Manley_SMI/results/GEx/hallmark_pathway_list.rds")
all_fov = readRDS("Manley_SMI/results/GEx/hallmark_pathway_list.rds")

#primary tumor (no sarcomatoid) treatment classes
naive_fovs = tumor %>% filter(Pretreatment.IO == "Treatment Naive", Sarcomatoid == "No") %>% pull(unique_fov) %>% unique()
postIO_fovs = tumor %>% filter(Pretreatment.IO == "Received IO") %>% pull(unique_fov) %>% unique() %>% grep("RCC5_19", ., value = T, invert = T)

naive_phenotype_scores = lapply(setNames(phenotypes, phenotypes), function(cell){
  lapply(naive_fovs, function(fov){
    all_fov[[fov]]$Phenotype_sets[[cell]]
  }) %>%
    do.call(bind_rows, .) %>%
    colMeans(na.rm = T)
}) %>% do.call(bind_rows, .)
rownames(naive_phenotype_scores) = paste0(phenotypes, " - Treatment Naive")

postIO_phenotype_scores = lapply(setNames(phenotypes, phenotypes), function(cell){
  d = lapply(postIO_fovs, function(fov){
    all_fov[[fov]]$Phenotype_sets[[cell]]
  }) %>%
    do.call(bind_rows, .) %>%
    colMeans(na.rm = T)
  if(length(d) == 0){
    d = rep(0, 50)
    names(d) = names(gene_set_list)
    return(d)
  }
  return(d)
}) %>% do.call(bind_rows, .)
rownames(postIO_phenotype_scores) = paste0(phenotypes, " - Received IO")

#heatmap
top_annon = HeatmapAnnotation(df = data.frame(Treatment = c(rep("Treatment Naive", length(phenotypes)),
                                                            rep("Received IO", length(phenotypes)))),
                              col = list(Treatment = c("Treatment Naive" = "black", "Received IO" = "darkgray")))
mat = rbind(naive_phenotype_scores,
            postIO_phenotype_scores) %>% t()
rownames(mat) = rownames(mat) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("_", " ", .)
ht = Heatmap(mat, top_annotation = top_annon, name = "GSVA Hallmark\nScore", 
             row_names_max_width = max_text_width(rownames(mat)), column_names_max_height = max_text_width(colnames(mat)),
             column_labels = colnames(mat) %>% gsub("\\ - .*", "", .), column_title = "Cell Level Hallmark GSEA (Averaged) - Non-Sarcomatoid Tumor",
             col = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("blue", "white", "red")))

pdf("Manley_SMI/results/GEx/primary_io/tumor/tumor_primary_pre-postIO_cell-level-GSEA-averaged.pdf", height = 10, width = 15)
ht
dev.off()


#primary stroma (no sarcomatoid) treatment classes
naive_fovs = stroma %>% filter(Pretreatment.IO == "Treatment Naive", Sarcomatoid == "No") %>% pull(unique_fov) %>% unique()
postIO_fovs = stroma %>% filter(Pretreatment.IO == "Received IO") %>% pull(unique_fov) %>% unique() %>% grep("RCC5_20", ., value = T, invert = T)

naive_phenotype_scores = lapply(setNames(phenotypes, phenotypes), function(cell){
  lapply(naive_fovs, function(fov){
    all_fov[[fov]]$Phenotype_sets[[cell]]
  }) %>%
    do.call(bind_rows, .) %>%
    colMeans(na.rm = T)
}) %>% do.call(bind_rows, .)
rownames(naive_phenotype_scores) = paste0(phenotypes, " - Treatment Naive")

postIO_phenotype_scores = lapply(setNames(phenotypes, phenotypes), function(cell){
  d = lapply(postIO_fovs, function(fov){
    all_fov[[fov]]$Phenotype_sets[[cell]]
  }) %>%
    do.call(bind_rows, .) %>%
    colMeans(na.rm = T)
  if(length(d) == 0){
    d = rep(0, 50)
    names(d) = names(gene_set_list)
    return(d)
  }
  return(d)
}) %>% do.call(bind_rows, .)
rownames(postIO_phenotype_scores) = paste0(phenotypes, " - Received IO")

#heatmap
top_annon = HeatmapAnnotation(df = data.frame(Treatment = c(rep("Treatment Naive", length(phenotypes)),
                                                            rep("Received IO", length(phenotypes)))),
                              col = list(Treatment = c("Treatment Naive" = "black", "Received IO" = "darkgray")))
mat = rbind(naive_phenotype_scores,
            postIO_phenotype_scores) %>% t()
rownames(mat) = rownames(mat) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("_", " ", .)
ht = Heatmap(mat, top_annotation = top_annon, name = "GSVA Hallmark\nScore", 
             row_names_max_width = max_text_width(rownames(mat)), column_names_max_height = max_text_width(colnames(mat)),
             column_labels = colnames(mat) %>% gsub("\\ - .*", "", .), column_title = "Cell Level Hallmark GSEA (Averaged) - Non-Sarcomatoid Stroma",
             col = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("blue", "white", "red")))
ht

pdf("Manley_SMI/results/GEx/primary_io/stroma/stroma_primary_pre-postIO_cell-level-GSEA-averaged.pdf", height = 10, width = 15)
ht
dev.off()