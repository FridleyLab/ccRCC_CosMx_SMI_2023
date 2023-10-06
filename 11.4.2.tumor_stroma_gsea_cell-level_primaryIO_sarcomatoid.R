rm(list=ls())
# libraries ---------------------------------------------------------------

library(pbmcapply)
library(msigdbr)
library(GSVA)
library(ComplexHeatmap)
library(tidyverse)
#data
metadata = readRDS("Manley_SMI/data/final_dataframes/metadata_clinical_spatial.rds") %>%
  filter(fov <= 20, unique_fov != "RCC5_13") %>%
  mutate(Site = gsub(" $", "", Site)) %>%
  filter(Site == "Tumor") %>%
  mutate(Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received IO"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received IO")),
         Sarcomatoid = gsub(" ", "", Sarcomatoid),
         Sarcomatoid.Group = case_when(IT.Treatment.before.collection == "None" & Sarcomatoid == "No" ~ "Treatment Naive Non-Sarcomatoid",
                                       IT.Treatment.before.collection != "None" & Sarcomatoid == "No" ~ "Post IO Non-Sarcomatoid",
                                       IT.Treatment.before.collection != "None" & Sarcomatoid == "Yes" ~ "Post IO Sarcomatoid",
                                       T ~ "Pre IO Sarcomatoid"), 
         Sarcomatoid.Group = factor(Sarcomatoid.Group, levels = c("Treatment Naive Non-Sarcomatoid", "Post IO Non-Sarcomatoid",
                                                                  "Post IO Sarcomatoid", "Pre IO Sarcomatoid")))
#cell types
phenotypes = metadata$lasso_final %>% unique()
#get the stroma and tumor groups
tumor = metadata %>% filter(Source == "Tumor", Sarcomatoid.Group %in% c("Treatment Naive Non-Sarcomatoid", "Pre IO Sarcomatoid"))
stroma = metadata %>% filter(Source == "Stroma", Sarcomatoid.Group %in% c("Treatment Naive Non-Sarcomatoid", "Pre IO Sarcomatoid"))


# Gene Sets ---------------------------------------------------------------

all_gene_sets = msigdbr(species = "Homo sapiens") %>%
  filter(gs_cat == "H") %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
gene_set_list = split(all_gene_sets$gene_symbol, all_gene_sets$gs_name)

fovs = metadata$unique_fov %>% unique()

all_fov = readRDS("Manley_SMI/results/GEx/hallmark_pathway_list.rds")

#primary tumor (no sarcomatoid) treatment classes
naive_fovs = tumor %>% filter(Sarcomatoid.Group == "Treatment Naive Non-Sarcomatoid") %>% pull(unique_fov) %>% unique()
sarc_fovs = tumor %>% filter(Sarcomatoid.Group == "Pre IO Sarcomatoid") %>% pull(unique_fov) %>% unique()

naive_phenotype_scores = lapply(setNames(phenotypes, phenotypes), function(cell){
  d = lapply(naive_fovs, function(fov){
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
rownames(naive_phenotype_scores) = paste0(phenotypes, " - Treatment Naive Non-Sarcomatoid")

sarc_phenotype_scores = lapply(setNames(phenotypes, phenotypes), function(cell){
  d = lapply(sarc_fovs, function(fov){
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
rownames(sarc_phenotype_scores) = paste0(phenotypes, " - Pre IO Sarcomatoid")

#heatmap
top_annon = HeatmapAnnotation(df = data.frame(Sample = c(rep("Treatment Naive Non-Sarcomatoid", length(phenotypes)),
                                                            rep("Pre IO Sarcomatoid", length(phenotypes)))),
                              col = list(Sample = c("Treatment Naive Non-Sarcomatoid" = "black", "Pre IO Sarcomatoid" = "darkgray")))
mat = rbind(naive_phenotype_scores,
            sarc_phenotype_scores) %>% t()
rownames(mat) = rownames(mat) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("_", " ", .)
ht = Heatmap(mat, top_annotation = top_annon, name = "GSVA Hallmark\nScore", 
        row_names_max_width = max_text_width(rownames(mat)), column_names_max_height = max_text_width(colnames(mat)),
        column_labels = colnames(mat) %>% gsub("\\ - .*", "", .), column_title = "Cell Level Hallmark GSEA (Averaged) - Sarcomatoid Tumor",
        col = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("blue", "white", "red")))

pdf("Manley_SMI/results/GEx/primary_io/sarcomatoid_tumor/sarcomatoid_tumor_primary_pre-non-sarc-postIO_cell-level-GSEA-averaged.pdf", height = 10, width = 15)
ht
dev.off()


#primary stroma (no sarcomatoid) treatment classes
naive_fovs = stroma %>% filter(Sarcomatoid.Group == "Treatment Naive Non-Sarcomatoid") %>% pull(unique_fov) %>% unique()
sarc_fovs = stroma %>% filter(Sarcomatoid.Group == "Pre IO Sarcomatoid") %>% pull(unique_fov) %>% unique()

naive_phenotype_scores = lapply(setNames(phenotypes, phenotypes), function(cell){
  d = lapply(naive_fovs, function(fov){
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
rownames(naive_phenotype_scores) = paste0(phenotypes, " - Treatment Naive Non-Sarcomatoid")

sarc_phenotype_scores = lapply(setNames(phenotypes, phenotypes), function(cell){
  d = lapply(sarc_fovs, function(fov){
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
rownames(sarc_phenotype_scores) = paste0(phenotypes, " - Pre IO Sarcomatoid")

#heatmap
top_annon = HeatmapAnnotation(df = data.frame(Sample = c(rep("Treatment Naive Non-Sarcomatoid", length(phenotypes)),
                                                            rep("Pre IO Sarcomatoid", length(phenotypes)))),
                              col = list(Sample = c("Treatment Naive Non-Sarcomatoid" = "black", "Pre IO Sarcomatoid" = "darkgray")))
mat = rbind(naive_phenotype_scores,
            sarc_phenotype_scores) %>% t()
rownames(mat) = rownames(mat) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("_", " ", .)
ht = Heatmap(mat, top_annotation = top_annon, name = "GSVA Hallmark\nScore", 
             row_names_max_width = max_text_width(rownames(mat)), column_names_max_height = max_text_width(colnames(mat)),
             column_labels = colnames(mat) %>% gsub("\\ - .*", "", .), column_title = "Cell Level Hallmark GSEA (Averaged) - Sarcomatoid Stroma",
             col = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("blue", "white", "red")))
ht

pdf("Manley_SMI/results/GEx/primary_io/sarcomatoid_stroma/sarcomatoid_stroma_primary_pre-non-sarc-postIO_cell-level-GSEA-averaged.pdf", height = 10, width = 15)
ht
dev.off()
