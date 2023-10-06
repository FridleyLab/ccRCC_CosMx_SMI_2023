rm(list=ls())
#libraries
library(tidyverse)
#library(spatialGE)
library(msigdbr) #7.5.1
#library(pbmcapply)
library(parallel)
library(ComplexHeatmap)

#data
#scaled_expression = readRDS("Manley_SMI/data/final_dataframes/sct_scaled_data.rds")
metadata = readRDS("Manley_SMI/data/final_dataframes/metadata_clinical_spatial.rds") %>%
  filter(fov <= 20, unique_fov != "RCC5_13") %>%
  mutate(Site = gsub(" $", "", Site)) %>%
  filter(Site == "Tumor") %>%
  mutate(Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received IO"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received IO"))) %>%
  mutate(Sarcomatoid = gsub(" $", "", Sarcomatoid))
#onnly includes primary tumor enrichment scores
all_fov = readRDS("Manley_SMI/results/GEx/hallmark_pathway_list.rds")

cores = 8
min_genes = 5
num_sds = 1
min_units = 20
reps = 1000

all_cell_enrichment = lapply(all_fov, function(x){
  x$Cell_Sets %>%
    t() %>% data.frame(check.names = F)
}) %>%
  do.call(bind_rows, .)

# #gene sets
# all_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
# msigdbr_list = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)

#loop each FOV
slide_fovs = unique(metadata$unique_fov)

res = mclapply(slide_fovs, function(sf){
  cat(sf, "\n")
  tissue_spots = metadata$id[metadata$unique_fov == sf]
  exp = all_cell_enrichment[tissue_spots,] %>% t()
  
  coords_df = metadata %>%
    filter(unique_fov == sf) %>%
    select(id, CenterX_local_px, CenterY_local_px)
  coords_df = coords_df[match(colnames(exp), coords_df$id),]
  
  pathways = row.names(exp)
  
  pval_ls = mclapply(setNames(pathways, pathways), function(pw){
    #need to get pathway average
    cat("\t", pw, "\n")
    pw_avg_exp = exp[pw,]
    avg = mean(pw_avg_exp)
    std = sd(pw_avg_exp)
    exp_thresh = avg + (num_sds * std)
    
    high_spots_bc = names(which(pw_avg_exp >= exp_thresh))
    
    # Are there at least 'min_units' number of cells/spots?
    if (length(high_spots_bc) >= min_units) {
      # Compute the distance between these spots' XY coordinates with high expression
      coords_high_spots = coords_df[coords_df$id %in% high_spots_bc,]
      distances_high_spots = as.matrix(dist(coords_high_spots[, c('CenterX_local_px', 'CenterY_local_px')], method =
                                                      'euclidean')) # Distance computation
      distances_high_spots[upper.tri(distances_high_spots)] = 0 # Make upper half zero to avoid sum of distances twice
      sum_high_distances = sum(distances_high_spots)
      
      # Compute random distance permutations
      
      #collect size of true spots above SD in data
      size = length(high_spots_bc)
      #get range of rows to sample
      rs = 1:nrow(coords_df)
      #calculate master distance matrix of all FOV cells
      tmp = as.matrix(dist(coords_df[,c("CenterX_local_px", "CenterY_local_px")], method = "euclidean"))
      #for number of permutations, sample true positive number of cells from FOV cells, then sum the distances dividing by 2
      sum_rand_distances2 = mclapply(1:reps, function(rep){
        vs = sample(x = rs, size = size, replace = F)
        sum(tmp[vs,vs]/2)
      }) %>% unlist()
      
      
      # Compute p-value
      # Ho: sum of observed distances is larger than sum of random distances
      # Ha: sum of observed distances is smaller than sum of random distances
      count_test = sum(sum_rand_distances2 < sum_high_distances) # Times observed dist was higher than random dists
      p_val = count_test / reps
      
    } else{
      p_val = NA
    } # Close test minimum spots # Close test minimum genes
    
    # Get results in data frame form
    pval_tmp = tibble::tibble(
      gene_set = pw,
      p_value = p_val
    )
    
    return(pval_tmp)
  }) %>%
    do.call(bind_rows, .)
  
  
  return(pval_ls)
}, mc.cores = 32, mc.allow.recursive = T)
names(res) = slide_fovs

# saveRDS(res, "Manley_SMI/data/pathway_enrichment/STEnrich_output_gsva.rds")
res = readRDS("Manley_SMI/data/pathway_enrichment/STEnrich_output_gsva.rds")
res = res[!(names(res) %in% c("RCC5_19", "RCC5_20"))]

joined_significance = lapply(names(res), function(nam){
  d = res[[nam]]
  colnames(d)[2] = nam
  return(data.table::data.table(d))
}) %>% Reduce(merge, .) %>%
  data.frame(check.names = F) %>%
  column_to_rownames("gene_set") %>%
  t() %>%
  data.frame(check.names = F)

joined_significance2 = -log10(joined_significance+0.0001)
joined_significance2_tmp = joined_significance2 %>%
  rownames_to_column("unique_fov") %>%
  right_join(metadata %>%
               select(Source, Sarcomatoid, unique_fov, Pretreatment.IO) %>%
               distinct(), .) %>%
  right_join(metadata %>%
               group_by(unique_fov) %>%
               summarise(Tumor_Percent = sum(lasso_final == "Tumor") / n() * 100), .)

mat = joined_significance2_tmp %>%
  select(-Source, -Sarcomatoid, -Tumor_Percent, -Pretreatment.IO) %>%
  column_to_rownames("unique_fov") %>% 
  as.matrix()
colnames(mat) = colnames(mat) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("_", " ", .)
top_anno = HeatmapAnnotation(df = joined_significance2_tmp %>% select(1, 3, 4, 5) %>% column_to_rownames("unique_fov"),
                             Tumor_Percent = anno_barplot(joined_significance2_tmp %>% pull(Tumor_Percent, unique_fov)),
                             col = list(Source = c("Tumor" = "blue", "Stroma" = "red"),
                                        Sarcomatoid = c("Yes" = "orange", "No" = "green"),
                                        Pretreatment.IO = c("Treatment Naive" = "skyblue", "Received IO" = "magenta")))
col_fun = circlize::colorRamp2(c(0, 1, 5), c("blue", "white", "red"))
ht = Heatmap(t(mat), top_annotation = top_anno, cluster_column_slices = TRUE, column_split = joined_significance2_tmp$Source,
        col = col_fun, name = "-log10(p)", row_names_max_width = max_text_width(colnames(mat)), 
        column_title = "STenrich - Clustering of Enrichment Scores"); ht

pdf("Manley_SMI/results/figures/Enrichment/STenrich_msigdb_gsva-scores.pdf", width = 15, height = 12)
ht
dev.off()


# non-sarcomatoid ---------------------------------------------------------
res = readRDS("Manley_SMI/data/pathway_enrichment/STEnrich_output_gsva.rds")
res = res[!(names(res) %in% c("RCC5_19", "RCC5_20"))]

joined_significance = lapply(names(res), function(nam){
  d = res[[nam]]
  colnames(d)[2] = nam
  return(data.table::data.table(d))
}) %>% Reduce(merge, .) %>%
  data.frame(check.names = F) %>%
  column_to_rownames("gene_set") %>%
  t() %>%
  data.frame(check.names = F)

joined_significance2 = -log10(joined_significance+0.0001)
joined_significance2_tmp = joined_significance2 %>%
  rownames_to_column("unique_fov") %>%
  right_join(metadata %>%
               select(Source, Sarcomatoid, unique_fov, Pretreatment.IO) %>%
               distinct(), .) %>%
  right_join(metadata %>%
               group_by(unique_fov) %>%
               summarise(Tumor_Percent = sum(lasso_final == "Tumor") / n() * 100), .) %>%
  filter(Sarcomatoid == "No")

mat = joined_significance2_tmp %>%
  select(-Source, -Sarcomatoid, -Tumor_Percent, -Pretreatment.IO) %>%
  column_to_rownames("unique_fov") %>% 
  as.matrix()
colnames(mat) = colnames(mat) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("_", " ", .)
top_anno = HeatmapAnnotation(df = joined_significance2_tmp %>% select(1, 3, 4, 5) %>% column_to_rownames("unique_fov"),
                             Tumor_Percent = anno_barplot(joined_significance2_tmp %>% pull(Tumor_Percent, unique_fov)),
                             col = list(Source = c("Tumor" = "blue", "Stroma" = "red"),
                                        Sarcomatoid = c("Yes" = "orange", "No" = "green"),
                                        Pretreatment.IO = c("Treatment Naive" = "skyblue", "Received IO" = "magenta")))
col_fun = circlize::colorRamp2(c(0, 1, 5), c("blue", "white", "red"))
ht = Heatmap(t(mat), top_annotation = top_anno, cluster_column_slices = TRUE, column_split = joined_significance2_tmp$Source,
             col = col_fun, name = "-log10(p)", row_names_max_width = max_text_width(colnames(mat)), 
             column_title = "STenrich - Clustering of Enrichment Scores"); ht

pdf("Manley_SMI/results/figures/Enrichment/STenrich_msigdb_gsva-scores_nonsarcomatoid.pdf", width = 12, height = 12)
ht
dev.off()


# Pre-Post Tumor ----------------------------------------------------------

res = readRDS("Manley_SMI/data/pathway_enrichment/STEnrich_output_gsva.rds")

joined_significance = lapply(names(res), function(nam){
  d = res[[nam]]
  colnames(d)[2] = nam
  return(data.table::data.table(d))
}) %>% Reduce(merge, .) %>%
  data.frame(check.names = F) %>%
  column_to_rownames("gene_set") %>%
  t() %>%
  data.frame(check.names = F)

joined_significance2 = -log10(joined_significance+0.0001)
joined_significance2_tmp = joined_significance2 %>%
  rownames_to_column("unique_fov") %>%
  right_join(metadata %>%
               select(Source, Sarcomatoid, unique_fov, Pretreatment.IO) %>%
               distinct(), .) %>%
  right_join(metadata %>%
               group_by(unique_fov) %>%
               summarise(Tumor_Percent = sum(lasso_final == "Tumor") / n() * 100), .) %>%
  filter(Source == "Tumor") %>%
  filter(!(Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive"))


mat = joined_significance2_tmp %>%
  select(-Source, -Sarcomatoid, -Tumor_Percent, -Pretreatment.IO) %>%
  column_to_rownames("unique_fov") %>% 
  as.matrix()
colnames(mat) = colnames(mat) %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("_", " ", .)
top_anno = HeatmapAnnotation(df = joined_significance2_tmp %>% select(1, 3, 4, 5) %>% column_to_rownames("unique_fov"),
                             Tumor_Percent = anno_barplot(joined_significance2_tmp %>% pull(Tumor_Percent, unique_fov)),
                             col = list(Source = c("Tumor" = "blue", "Stroma" = "red"),
                                        Sarcomatoid = c("Yes" = "orange", "No" = "green"),
                                        Pretreatment.IO = c("Treatment Naive" = "skyblue", "Received IO" = "magenta")))
col_fun = circlize::colorRamp2(c(0, 1, 5), c("blue", "white", "red"))
ht = Heatmap(t(mat), top_annotation = top_anno, cluster_column_slices = TRUE, column_split = joined_significance2_tmp$Source,
             col = col_fun, name = "-log10(p)", row_names_max_width = max_text_width(colnames(mat)), 
             column_title = "STenrich - Clustering of Enrichment Scores"); ht

pdf("Manley_SMI/results/figures/Enrichment/STenrich_msigdb_gsva-scores_pre-postIO_Tumor.pdf", width = 10, height = 12)
ht
dev.off()


# Saving dataframe --------------------------------------------------------

clinical = metadata %>%
  select(Site:tissue, slide_fov) %>%
  distinct() %>%
  full_join(joined_significance %>% rownames_to_column("unique_fov"))
write.csv(clinical, "Manley_SMI/results/GEx/STenrich_hallmark_pathways.csv")











