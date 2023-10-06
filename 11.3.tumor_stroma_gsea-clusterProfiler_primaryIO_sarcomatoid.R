rm(list=ls())
#libraries
library(tidyverse)
library(spatstat)
library(lmerTest)
library(ggpubr)
library(pbmcapply)
library(msigdbr)
library(clusterProfiler)
library(openxlsx)
library(ComplexHeatmap)

#data
scaled_expression = readRDS("Manley_SMI/data/final_dataframes/sct_scaled_data.rds")
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
                                                                  "Post IO SaSarcomatoidrcoma", "Pre IO Sarcomatoid")))
#cell types
phenotypes = metadata$lasso_final %>% unique()
#get the stroma and tumor groups
tumor = metadata %>% filter(Source == "Tumor", Sarcomatoid.Group %in% c("Treatment Naive Non-Sarcomatoid", "Pre IO Sarcomatoid"))
stroma = metadata %>% filter(Source == "Stroma", Sarcomatoid.Group %in% c("Treatment Naive Non-Sarcomatoid", "Pre IO Sarcomatoid"))

# Tumor Cells -------------------------------------------------------------
#just quick test
tumor_output = pbmclapply(phenotypes, function(cell){
  cat(cell, "\n")
  dat = tumor %>%
    filter(lasso_final == cell)
  all_expression = scaled_expression[row.names(scaled_expression) %in% dat$id,]
  if(is(all_expression, "numeric")) return(NULL) #catches those that only have a low cell counts and fail subsetting or modeling
  if(nrow(all_expression)<2) return(NULL) 
  variances = apply(all_expression, 2, var)
  variable_genes = names(sort(variances, decreasing = T))[1:20]
  dat = dat %>% 
    left_join(all_expression[,variable_genes] %>%
                data.frame(check.names = F) %>%
                rownames_to_column("id"), by = "id")
  gene_mods = lapply(variable_genes, function(gene){
    form = formula(paste0("`", gene, "` ~ Sarcomatoid.Group + (1|unique_fov)"))
    t = try(lmer(form, data = dat) %>% 
      summary() %>% 
      coef() %>% 
      data.frame(check.names = F) %>%
      rownames_to_column("Level") %>% 
      filter(grepl("Sarco", Level)) %>%
      mutate(Gene = gene,
             Phenotype = cell))
    if(is(t, "try-error")){
      return(NULL)
    }
    return(t)
  }) %>% 
    do.call(bind_rows, .)
  
  if(nrow(gene_mods) > 0) {
    gene_mods = gene_mods %>%
    mutate(padj = p.adjust(`Pr(>|t|)`, method = "fdr"), .before = "Gene") %>%
    arrange(padj)
  } else {
    return(NULL)
  }
  
  gene_plots = lapply(gene_mods$Gene, function(gene){
    dat %>%
      ggplot(aes(x = Sarcomatoid.Group, y = get(gene))) +
      geom_violin() +
      geom_jitter(alpha = 0.2, width = 0.3) +
      stat_summary(fun.y=mean, geom="line", color="red", group = 'Sarcomatoid.Group') +
      labs(title = paste0(gene),
           subtitle = paste0("(pdj = ", signif(gene_mods$padj[gene_mods$Gene == gene], digits = 3), ")"),
           y = gene) + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  })
  
  pl = ggarrange(plotlist = gene_plots, nrow = 4, ncol = 5) %>%
    annotate_figure(., top = text_grob(cell, face = "bold", size = 14))
  
  return(list(`Model Stats` = gene_mods,
              Plots = pl))
}, mc.cores = 5)

names(tumor_output) = phenotypes
tumor_output = tumor_output[!sapply(tumor_output, is.null)]

tmp = lapply(seq(tumor_output), function(num){
  pdf(paste0("Manley_SMI/results/GEx/primary_io/sarcomatoid_tumor/tumor_primary_pre-postIO_", names(tumor_output)[num],".pdf"), height = 10, width = 15)
  print(tumor_output[[num]]$Plots)
  dev.off()
})


# Tumor All Genes Differential --------------------------------------------

#getting p vals for all genes for all phenotypes
tumor_output = pbmclapply(phenotypes, function(cell){
  cat(cell, "\n")
  dat = tumor %>%
    filter(lasso_final == cell)
  all_expression = scaled_expression[row.names(scaled_expression) %in% dat$id,]
  if(is(all_expression, "numeric")) return(NULL) #catches those that only have a low cell counts and fail subsetting or modeling
  if(nrow(all_expression)<2) return(NULL) 
  variable_genes = colnames(all_expression)
  dat = dat %>% 
    left_join(all_expression[] %>%
                data.frame(check.names = F) %>%
                rownames_to_column("id"), by = "id")
  gene_mods = lapply(variable_genes, function(gene){
    form = formula(paste0("`", gene, "` ~ Sarcomatoid.Group + (1|unique_fov)"))
    t = try(lmer(form, data = dat) %>% 
              summary() %>% 
              coef() %>% 
              data.frame(check.names = F) %>%
              rownames_to_column("Level") %>% 
              filter(grepl("Sarcoma", Level)) %>%
              mutate(Gene = gene,
                     Phenotype = cell))
    if(is(t, "try-error")){
      return(NULL)
    }
    return(t)
  }) %>% 
    do.call(bind_rows, .)
  
  if(nrow(gene_mods) > 0) {
    gene_mods = gene_mods %>%
      mutate(padj = p.adjust(`Pr(>|t|)`, method = "fdr"), .before = "Gene") %>%
      arrange(padj)
  } else {
    return(NULL)
  }
  
  return(list(`Model Stats` = gene_mods))
}, mc.cores = 5)

names(tumor_output) = phenotypes
#saveRDS(tumor_output, "Manley_SMI/results/GEx/primary_io/sarcomatoid_tumor/tumor_cell-type_gex.rds")

#enrichment
tumor_output = tumor_output[!sapply(tumor_output, is.null)]
all_gene_sets = msigdbr(species = "Homo sapiens") %>%
  filter(gs_cat == "H") %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
test = lapply(tumor_output, function(cell){
  up = try(clusterProfiler::enricher(cell$`Model Stats` %>% filter(padj < 0.1, Estimate > 0) %>% pull(Gene), TERM2GENE = all_gene_sets)@result)
  down = try(clusterProfiler::enricher(cell$`Model Stats` %>% filter(padj < 0.1, Estimate < 0) %>% pull(Gene), TERM2GENE = all_gene_sets)@result)
  if(is(up, "try-error") & is(down, "try-error")){
    return(NULL)
  } else if(is(up, "try-error")){
    return(down %>%
             mutate(Direction = "Up in Treatment Naive Non-Sarcomatoid") %>%
             mutate(Phenotype = unique(cell$`Model Stats`$Phenotype)))
  } else if(is(down, "try-error")){
    return(up %>%
             mutate(Direction = "Up in Pre IO Sarcomatoid") %>%
             mutate(Phenotype = unique(cell$`Model Stats`$Phenotype)))
  }
  out = bind_rows(up %>%
    mutate(Direction = "Up in Pre IO Sarcomatoid"), down %>%
    mutate(Direction = "Up in Treatment Naive Non-Sarcomatoid")) %>%
    mutate(Phenotype = unique(cell$`Model Stats`$Phenotype))
  return(out)
})

#save tumor outputs
wb = createWorkbook()
addWorksheet(wb, "Differential Genes")
addWorksheet(wb, "Enrichment of p<0.1 Genes")
writeData(wb, 1, lapply(tumor_output, function(x) x$`Model Stats`) %>% do.call(bind_rows, .) %>% 
            mutate(`Higher In` = ifelse(Estimate >= 0, "Pre IO Sarcomatoid", "Treatment Naive Non-Sarcomatoid"), .before = Estimate))
writeData(wb, 2, do.call(bind_rows, test))
saveWorkbook(wb, "Manley_SMI/results/GEx/primary_io/sarcomatoid_tumor/sarcomatoid_tumor_Differential_gene_expression_no-spatial.xlsx", overwrite = T)

#visualize

path_prep = do.call(bind_rows, test) %>% select(Description, p.adjust, Direction, Phenotype)
mat_rio = path_prep %>%
  filter(Direction == "Up in Pre IO Sarcomatoid") %>%
  select(Description, p.adjust, Phenotype) %>%
  spread("Phenotype", "p.adjust", fill = 1) %>%
  mutate(Description = gsub("HALLMARK_", "", Description),
         Description = gsub("_", " ", Description))
mat_rio[phenotypes[!phenotypes %in% colnames(mat_rio)]] = 1
mat_rio = mat_rio %>%
  column_to_rownames("Description") %>% as.matrix() %>% t()
ht_rio = Heatmap(-log10(mat_rio), name = "-log10 FDR", column_title = "Gene Set Enrichment in Tumor FOVs of Pre-IO  Sarcomatoid",
                 column_names_max_height = max_text_width(colnames(mat_rio)))

mat_tn = path_prep %>%
  filter(Direction == "Up in Treatment Naive Non-Sarcomatoid") %>%
  select(Description, p.adjust, Phenotype) %>%
  spread("Phenotype", "p.adjust", fill = 1) %>%
  mutate(Description = gsub("HALLMARK_", "", Description),
         Description = gsub("_", " ", Description))
mat_tn[phenotypes[!phenotypes %in% colnames(mat_tn)]] = 1
mat_tn = mat_tn %>%
  column_to_rownames("Description") %>% as.matrix() %>% t()
ht_tn = Heatmap(-log10(mat_tn), name = "-log10 FDR", column_title = "Gene Set Enrichment in Treatment Naive Non-Sarcoma Tumor FOVs", 
                column_names_max_height = max_text_width(colnames(mat_tn)))

pdf("Manley_SMI/results/GEx/primary_io/sarcomatoid_tumor/sarcomatoid_tumor_pre-post-io_geneset-pvalue-heatmaps.pdf", width = 13, height = 10)
ht_rio
ht_tn
dev.off()



# Stroma ------------------------------------------------------------------
#just quick test
stroma_output = pbmclapply(phenotypes, function(cell){
  cat(cell, "\n")
  dat = stroma %>%
    filter(lasso_final == cell)
  all_expression = scaled_expression[row.names(scaled_expression) %in% dat$id,]
  if(is(all_expression, "numeric")) return(NULL) #catches those that only have a low cell counts and fail subsetting or modeling
  if(nrow(all_expression)<2) return(NULL) 
  variances = apply(all_expression, 2, var)
  variable_genes = names(sort(variances, decreasing = T))[1:20]
  dat = dat %>% 
    left_join(all_expression[,variable_genes] %>%
                data.frame(check.names = F) %>%
                rownames_to_column("id"), by = "id")
  gene_mods = lapply(variable_genes, function(gene){
    form = formula(paste0("`", gene, "` ~ Sarcomatoid.Group + (1|unique_fov)"))
    t = try(lmer(form, data = dat) %>% 
              summary() %>% 
              coef() %>% 
              data.frame(check.names = F) %>%
              rownames_to_column("Level") %>% 
              filter(grepl("Sarcoma", Level)) %>%
              mutate(Gene = gene,
                     Phenotype = cell))
    if(is(t, "try-error")){
      return(NULL)
    }
    return(t)
  }) %>% 
    do.call(bind_rows, .)
  
  if(nrow(gene_mods) > 0) {
    gene_mods = gene_mods %>%
      mutate(padj = p.adjust(`Pr(>|t|)`, method = "fdr"), .before = "Gene") %>%
      arrange(padj)
  } else {
    return(NULL)
  }
  
  gene_plots = lapply(gene_mods$Gene, function(gene){
    dat %>%
      ggplot(aes(x = Sarcomatoid.Group, y = get(gene))) +
      geom_violin() +
      geom_jitter(alpha = 0.2, width = 0.3) +
      stat_summary(fun.y=mean, geom="line", color="red", group = 'Sarcomatoid.Group') +
      labs(title = paste0(gene),
           subtitle = paste0("(pdj = ", signif(gene_mods$padj[gene_mods$Gene == gene], digits = 3), ")"),
           y = gene) + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  })
  
  pl = ggarrange(plotlist = gene_plots, nrow = 4, ncol = 5) %>%
    annotate_figure(., top = text_grob(cell, face = "bold", size = 14))
  
  return(list(`Model Stats` = gene_mods,
              Plots = pl))
}, mc.cores = 5)

names(stroma_output) = phenotypes
stroma_output = stroma_output[!sapply(stroma_output, is.null)]

tmp = lapply(seq(stroma_output), function(num){
  pdf(paste0("Manley_SMI/results/GEx/primary_io/sarcomatoid_stroma/stroma_primary_pre-postIO_", names(stroma_output)[num],".pdf"), height = 10, width = 15)
  print(stroma_output[[num]]$Plots)
  dev.off()
})

#getting p vals for all genes for all phenotypes
stroma_output = pbmclapply(phenotypes, function(cell){
  cat(cell, "\n")
  dat = stroma %>%
    filter(lasso_final == cell)
  all_expression = scaled_expression[row.names(scaled_expression) %in% dat$id,]
  variable_genes = colnames(all_expression)
  dat = dat %>% 
    left_join(all_expression[] %>%
                data.frame(check.names = F) %>%
                rownames_to_column("id"), by = "id")
  gene_mods = lapply(variable_genes, function(gene){
    form = formula(paste0("`", gene, "` ~ Sarcomatoid.Group + (1|unique_fov)"))
    t = try(lmer(form, data = dat) %>% 
              summary() %>% 
              coef() %>% 
              data.frame(check.names = F) %>%
              rownames_to_column("Level") %>% 
              filter(grepl("Sarco", Level)) %>%
              mutate(Gene = gene,
                     Phenotype = cell))
    if(is(t, "try-error")){
      return(NULL)
    }
    return(t)
  }) %>% 
    do.call(bind_rows, .)
  
  if(nrow(gene_mods) > 0) {
    gene_mods = gene_mods %>%
      mutate(padj = p.adjust(`Pr(>|t|)`, method = "fdr"), .before = "Gene") %>%
      arrange(padj)
  } else {
    return(NULL)
  }
  
  return(list(`Model Stats` = gene_mods))
}, mc.cores = 5)

names(stroma_output) = phenotypes
#saveRDS(stroma_output, "Manley_SMI/results/GEx/primary_io/sarcomatoid_stroma/stroma_cell-type_gex.rds")

#enrichment
stroma_output = stroma_output[!sapply(stroma_output, is.null)]
all_gene_sets = msigdbr(species = "Homo sapiens") %>%
  filter(gs_cat == "H") %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
stroma_pathays = lapply(stroma_output, function(cell){
  up = try(clusterProfiler::enricher(cell$`Model Stats` %>% filter(padj < 0.1, Estimate > 0) %>% pull(Gene), TERM2GENE = all_gene_sets)@result)
  down = try(clusterProfiler::enricher(cell$`Model Stats` %>% filter(padj < 0.1, Estimate < 0) %>% pull(Gene), TERM2GENE = all_gene_sets)@result)
  if(is(up, "try-error") & is(down, "try-error")){
    return(NULL)
  } else if(is(up, "try-error")){
    return(down %>%
             mutate(Direction = "Up in Treatment Naive Non-Sarcomatoid") %>%
             mutate(Phenotype = unique(cell$`Model Stats`$Phenotype)))
  } else if(is(down, "try-error")){
    return(up %>%
             mutate(Direction = "Up in Pre IO Sarcomatoid") %>%
             mutate(Phenotype = unique(cell$`Model Stats`$Phenotype)))
  }
  out = bind_rows(up %>%
                    mutate(Direction = "Up in Pre IO Sarcomatoid"), down %>%
                    mutate(Direction = "Up in Treatment Naive Non-Sarcomatoid")) %>%
    mutate(Phenotype = unique(cell$`Model Stats`$Phenotype))
  return(out)
})

#save stroma outputs
wb = createWorkbook()
addWorksheet(wb, "Differential Genes")
addWorksheet(wb, "Enrichment of p<0.1 Genes")
writeData(wb, 1, lapply(stroma_output, function(x) x$`Model Stats`) %>% do.call(bind_rows, .) %>% 
            mutate(`Higher In` = ifelse(Estimate >= 0, "Pre IO Sarcomatoid", "Treatment Naive Non-Sarcoma"), .before = Estimate))
writeData(wb, 2, do.call(bind_rows, stroma_pathays))
saveWorkbook(wb, "Manley_SMI/results/GEx/primary_io/sarcomatoid_stroma/sarcomatoid_stroma_Differential_gene_expression_no-spatial.xlsx", overwrite = T)

#visualize

path_prep = do.call(bind_rows, stroma_pathays) %>% select(Description, p.adjust, Direction, Phenotype)
mat_rio = path_prep %>%
  filter(Direction == "Up in Pre IO Sarcomatoid") %>%
  select(Description, p.adjust, Phenotype) %>%
  spread("Phenotype", "p.adjust", fill = 1) %>%
  mutate(Description = gsub("HALLMARK_", "", Description),
         Description = gsub("_", " ", Description))
mat_rio[phenotypes[!phenotypes %in% colnames(mat_rio)]] = 1
mat_rio = mat_rio %>%
  column_to_rownames("Description") %>% as.matrix() %>% t()
ht_rio = Heatmap(-log10(mat_rio), name = "-log10 FDR", column_title = "Gene Set Enrichment in Stroma FOVs of Pre IO Sarcomatoid",
                 column_names_max_height = max_text_width(colnames(mat_rio)))

mat_tn = path_prep %>%
  filter(Direction == "Up in Treatment Naive Non-Sarcomatoid") %>%
  select(Description, p.adjust, Phenotype) %>%
  spread("Phenotype", "p.adjust", fill = 1) %>%
  mutate(Description = gsub("HALLMARK_", "", Description),
         Description = gsub("_", " ", Description))
mat_tn[phenotypes[!phenotypes %in% colnames(mat_tn)]] = 1
mat_tn = mat_tn %>%
  column_to_rownames("Description") %>% as.matrix() %>% t()
ht_tn = Heatmap(-log10(mat_tn), name = "-log10 FDR", column_title = "Gene Set Enrichment in Treatment Naive Non-Sarcoma Stroma FOVs", 
                column_names_max_height = max_text_width(colnames(mat_tn)))

pdf("Manley_SMI/results/GEx/primary_io/sarcomatoid_stroma/sarcomatoid_stroma_pre-post-io_geneset-pvalue-heatmaps.pdf", width = 13, height = 10)
ht_rio
ht_tn
dev.off()
