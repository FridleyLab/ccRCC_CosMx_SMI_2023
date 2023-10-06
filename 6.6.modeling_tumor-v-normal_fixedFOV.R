library(tidyverse)
library(Seurat)
library(plotly)
library(RColorBrewer)
library(ggpubr)
# library(patchwork)
#library(future)
#plan("multisession", workers = 4)

manley_data = readRDS("Manley_SMI/results/rds/InSituType/manley_data_manual_TandMonocytes.rds")
manley_data@meta.data = manley_data@meta.data %>%
  mutate(final_insitutype_simplified = case_when(final_insitutype %in% c("Ascending vasa recta endothelium",
                                                                         "Descending vasa recta endothelium",
                                                                         "Thick ascending limb of Loop of Henle")
                                                 ~ "Vasa recta endothelium",
                                                 final_insitutype %in% c("Connecting tubule",
                                                                         "Distinct proximal tubule 1",
                                                                         "Distinct proximal tubule 2",
                                                                         "Proliferating Proximal Tubule")
                                                 ~ "Proximal tubule",
                                                 final_insitutype %in% c("Peritubular capillary endothelium 1",
                                                                         "Peritubular capillary endothelium 2")
                                                 ~ "Capillary endothelium",
                                                 final_insitutype %in% c("Indistinct intercalated cell",
                                                                         "Type A intercalated cell",
                                                                         "Type B intercalated cell")
                                                 ~ "Intercalated Cell",
                                                 T ~ final_insitutype))

meta_data = manley_data@meta.data
meta_data = meta_data %>%
  mutate(cell_class = case_when(final_insitutype_simplified %in% c("Fibroblast", "Myofibroblast") ~ "tissue forming",
                                final_insitutype_simplified %in% c("M2 macrophage (CD163)", "B cell", "M1 Macrophage (STAT1)", "Neutrophil",
                                                                   "CD8 T cell", "gdT cell", "Regulatory T cell", "Naive CD4 T cell", "Intermediate monocyte",
                                                                   "Plasmacytoid dendritic cell", "Myeloid DC", "Naive T cell", "NK cell",
                                                                   "Naive B cell", "Mast cell", "Non-classical monocyte", "Classical monocyte") ~
                                  "immune",
                                final_insitutype_simplified %in% c("Possibly Mid-Rep (Misc. Cells)", "Pelvic epithelium") ~ "other",
                                T ~ "tissue"),
         tumor_baseline = case_when(cell_class == "tissue" & fov == "5" & Slide_name == "RCC3" ~ "Tumor",
                                    cell_class == "tissue" & fov == "11" & Slide_name == "RCC3" ~ "Tumor",
                                    cell_class == "tissue" & fov == "15" & Slide_name == "RCC3" ~ "Tumor",
                                    cell_class == "tissue" & fov == "3" & Slide_name == "RCC4" ~ "Tumor",
                                    cell_class == "tissue" & fov == "11" & Slide_name == "RCC4" ~ "Tumor",
                                    T ~ final_insitutype_simplified))
manley_data@meta.data = meta_data

#FOVs with distinct proximal tubule clusters
#tumor FOVs
#3-1, 3-5, 3-7, 3-11, 3-13, 3-21, 4-1, 4-3, 4-7, 4-13, 5-1, 5-5, 5-7, 5-9, 5-19, 5-21
#stroma FOVs
#3-2, 3-20, 3-22, 4-20, 5-6, 5-10, 5-20, 5-22, 

#making new metadata column
manley_data@meta.data = manley_data@meta.data %>%
  mutate(tumor_class_for_DEGS = case_when((Slide_name == "RCC3" &
                                             fov %in% c(1, 5, 7, 11, 13, 21) |
                                             Slide_name == "RCC4" &
                                             fov %in% c(1, 3, 7, 13) | 
                                             Slide_name == "RCC5" &
                                             fov %in% c(1, 5, 7, 9, 19, 21)) &
                                            final_insitutype_simplified == "Proximal tubule" ~ "tumor_pt",
                                          (Slide_name == "RCC3" &
                                             fov %in% c(2, 20, 22) |
                                             Slide_name == "RCC4" &
                                             fov %in% c(20) | 
                                             Slide_name == "RCC5" &
                                             fov %in% c(6, 10, 20, 22)) &
                                            final_insitutype_simplified == "Proximal tubule" ~ "stroma_pt",
                                          T ~ "other"))
Idents(manley_data) = "tumor_class_for_DEGS"
pt_diff_genes = FindMarkers(manley_data, ident.1 = "tumor_pt", ident.2 = "stroma_pt")
pt_genes_up = pt_diff_genes %>%
  filter(avg_log2FC > 0.5) %>% 
  row.names()
pt_genes_down = pt_diff_genes %>%
  filter(avg_log2FC < -0.5) %>% 
  row.names()
pt_genes = c(pt_genes_up, pt_genes_down)

FeaturePlot(manley_data, features = pt_genes) #looks like VIMis everywhere, not just proximal tubule. apparently LFC >1.3 in tumor

#create data to plot with and determine fitness of gene
plt_dat = left_join(manley_data@meta.data %>% 
                      filter(cell_class == "tissue") %>%
                      rownames_to_column("slide_fov_cell"),
                    manley_data@assays$SCT@scale.data[c(pt_genes_up,pt_genes_down),] %>%
                      as.matrix() %>% t() %>% data.frame(check.names=F) %>% 
                      rownames_to_column("slide_fov_cell")) %>% 
  column_to_rownames("slide_fov_cell")
#plots definitely show that even proximal tubule cells and vasa recta endothelium on stroma FOVs have greater than zero expression
plt_dat %>% 
  ggplot() +
  geom_density(aes(x = VIM, group = annotation, color = annotation)) +
  facet_grid(final_insitutype_simplified~.)
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = VIM, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~.)

#looking at other genes up
#pdf("Manley_SMI/results/figures/classifying_tumor/top_10_DEGs_T_S_proximal_tubule-tissue_applied.pdf", height = 12, width = 10)
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = CCND1, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y") # 0 seems like okay threshold - above peak with mostly tumor above 0
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = CD24, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y") # 0 seems fair here, proximal tubule in stroma does show few cells above
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = DUSP1, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y") # main peak falls off before zero but proximal tubule peak is highest after - i think 0 again here
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = FGG, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y")# most cells are below 0, with small peak close 10 on proximal tubule
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = HILPDA, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y") # most cells below zero except tumor podocytes and proximal tubule
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = IGFBP3, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y") #similarly mainly tumor above 0, capillary, podo, proximal, and vasa recta in the tumor
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = SPINK1, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y") # almost everything below 0 except very very few podocytes and then a peak for proximal tubule at 10 0, about 1000 cells
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = VEGFA, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y") # 0 again for mostly tumor FOVs, stroma podo does have a buildup of 100 cells at 10 and stragglers in proximal tubule
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = VIM, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y") # dont like this one
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = IGFBP5, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y")
#dev.off()
#genes down
#pdf("Manley_SMI/results/figures/classifying_tumor/bottom_4_DEGs_T_S_proximal_tubule-tissue_applied.pdf", height = 12, width = 10)
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = IGFBP7, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y") # 0 seems like okay threshold - above peak with mostly tumor above 0
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = LTF, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y") # 0 seems fair here, proximal tubule in stroma does show few cells above
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = PIGR, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y") # main peak falls off before zero but proximal tubule peak is highest after - i think 0 again here
plt_dat %>% 
  ggplot() +
  geom_histogram(aes(x = S100A6, group = annotation, fill = annotation), alpha = 0.5, bins = 100) +
  facet_grid(final_insitutype_simplified~annotation, scales="free_y")# most cells are below 0, with small peak close 10 on proximal tubule
dev.off()
##try just usin if any of these genes have an expression greater than 0 ~ Tumor, else normal


# Adding annotation to metadata -------------------------------------------

# Polygon Plots with Tumor Definitions ------------------------------------
mat = manley_data@assays$SCT@scale.data %>% t() %>% data.frame(check.names = F) %>% dplyr::select(!!pt_genes)
wkn_tbl = full_join(mat %>% rownames_to_column("rowname"),
                    manley_data@meta.data %>% rownames_to_column("rowname"))
up = apply(wkn_tbl, 1, function(x) (TRUE %in% (as.numeric(x[pt_genes_up[-9]]) > 0)))
manley_data$pt10genes_0thresh = ifelse(wkn_tbl$cell_class == "tissue" &
                                         up == TRUE, "Tumor", manley_data$final_insitutype_simplified)
wkn_tbl$tubule_tumor_ptgenes0 = ifelse(wkn_tbl$cell_class == "tissue" &
                                         up == TRUE, "Tumor", wkn_tbl$final_insitutype_simplified)

#identify the files with the polygon information for all FOVs on all slides
polygons_files = list.files("Manley_SMI/data/SMI-0050_BrandonManley_Moffit/5 Raw data/", 
                            recursive = T,pattern = "polygons", full.names = T)
#now read in the polygon slides
polygon_data = lapply(polygons_files, read.csv, check.names = F)

#for merging need to add which slide (tissue) that the polygon data is for
polygon_data[[1]]$tissue = "RCC3"
polygon_data[[2]]$tissue = "RCC4"
polygon_data[[3]]$tissue = "RCC5"

#concatonate the polygon data into a single data frame and convert the fov and cellID to characters for merging with the metadata
polygon_data = polygon_data %>% 
  do.call(bind_rows, .) %>%
  mutate_at(c("cellID"), as.character)

wkn_tbl$cellID = wkn_tbl$cell_ID %>%
  gsub('.*_', '', .)

#prepping colors for phenotypes
phenotypes = wkn_tbl$tubule_tumor_ptgenes0 %>% unique() %>% sort()
wkn_tbl$tubule_tumor_ptgenes0 = factor(wkn_tbl$tubule_tumor_ptgenes0,
                                       levels = phenotypes[c(1, 3:4, 8, 11:15, 17:22, 24, 29, #immune
                                                             2, 5:7, 9:10, 16, 23, 25:28, 31, #tissue
                                                             30)]) #tumor
#create plots
immune_color_func = colorRampPalette(c("red", "blue"))
immune_colors = immune_color_func(17)
tissue_color_func = colorRampPalette(c("green", "orange"))
tissue_colors = tissue_color_func(13)
all_colors = c(immune_colors, tissue_colors, 'black')
names(all_colors) = levels(wkn_tbl$tubule_tumor_ptgenes0)

full_data = wkn_tbl %>%
  left_join(polygon_data) 
#now we can use the group_map to output a list of plots for all of the slide/fov combinations
plots = full_data %>%
  group_by(tissue, fov) %>%
  group_map(~{
    .x %>%
      ggplot() + 
      geom_polygon(aes(group = cellID, fill = tubule_tumor_ptgenes0, x = x_local_px, y = y_local_px),
                   color = "black", linewidth = 0.1) +
      theme_bw() +
      labs(title = paste0("Slide ", .y$tissue, " - FOV ", .y$fov)) +
      tune::coord_obs_pred() + 
      guides(fill = guide_legend(title = "InSituType Assignments\n(Manley Collapsed)",ncol = 1)) + 
      scale_fill_manual(values = all_colors)
  })


# #save the plots into PDF form into the figures folder
# pdf("Manley_SMI/results/figures/InSituType_Phenotypes_polygon_plots_ptgenes0-tumor.pdf", width = 14, height = 10)
# tmp = lapply(plots, print)
# dev.off()



# Modeling to predict tumor v normal --------------------------------------

model_dat = plt_dat %>% filter(tumor_class_for_DEGS != "other") %>%
  rownames_to_column("rowname") %>% 
  left_join(manley_data@assays$SCT@scale.data[row.names(pt_diff_genes),] %>%
              as.matrix() %>% t() %>% data.frame(check.names = F) %>% rownames_to_column("rowname"))
model_dat$model_class = ifelse(model_dat$annotation == "Tumor ",0,1)
tt_split = sample(c(1, 2), nrow(model_dat), prob = c(0.2, 0.8), replace = T)
train_dat = model_dat[tt_split == 2,c("model_class", row.names(pt_diff_genes))]
train_dat = train_dat[c(sample(which(train_dat$model_class == 0),
                               min(table(train_dat$model_class)), replace = F),
                        +                   which(train_dat$model_class == 1)),]
test_dat = model_dat[tt_split == 1,c("model_class",row.names(pt_diff_genes))]
#form = formula(paste("model_class ~", paste0(row.names(pt_diff_genes), collapse = "+")))
logit_mod = glm(model_class~., data = train_dat, family = "binomial")
#summary(logit_mod)
logit_mod_2 = glm(model_class~., data = subset(train_dat, select=-c(ADGRG1, HIF1A, NDRG1, COL18A1, CAV1,
                                                                    S100A10, TNFRSF12A, CLDN4, YBX3, CLU,
                                                                    SAA1, MMP7, RARRES2, NEAT1, PSAP, SPP1)), family = "binomial") #, SOD2, CXCL14, SPP1
summary(logit_mod_2)
sapply(list(`Full model` = logit_mod,
            `Reduced Model` = logit_mod_2), AIC)

data.frame(Coeff = logit_mod_2$coefficients) %>% 
  rownames_to_column("Gene") %>%
  arrange(desc(Coeff)) %>%
  pull(Gene) %>%
  paste0(., collapse = ", ")
# reduced model is actually better with lower AIC and less coefficients
model_genes_all = intersect(row.names(pt_diff_genes), colnames(logit_mod_2$data))
model_dat_all = model_dat[,c("model_class", model_genes_all)]
logit_model_all = glm(model_class ~ ., data = model_dat_all, family = "binomial")

gex = manley_data@assays$SCT@scale.data[model_genes_all,] %>%
  t() %>% data.frame(check.names = F)
#predict all cells from RCC3 1&2 tumor/stroma model
preds_all = predict(logit_model_all, gex, type = "response")

#Add to manley_data, metadata
predicted_df = data.frame("rowname" = names(preds_all),
                          Predicted = preds_all,
                          Thresholded = ifelse(manley_data$cell_class == 'tissue' &
                                                 preds_all < 0.5,
                                               "Tumor", manley_data$final_insitutype_simplified)) %>% 
  full_join(full_data)
manley_data$pt_manualAIC_optim = ifelse(manley_data$cell_class == 'tissue' &
                                          preds_all < 0.5,
                                        "Tumor", manley_data$final_insitutype_simplified)
#prepping colors for phenotypes
phenotypes = predicted_df$Thresholded %>% unique() %>% sort()
predicted_df$Thresholded = factor(predicted_df$Thresholded,
                                  levels = phenotypes[c(1, 3:4, 8, 11:15, 17:22, 24, 29, #immune
                                                        2, 5:7, 9:10, 16, 23, 25:28, 31, #tissue
                                                        30)]) #tumor


plots = predicted_df %>%
  group_by(tissue, fov) %>%
  group_map(~{
    .x %>%
      ggplot() + 
      geom_polygon(aes(group = cellID, fill = Thresholded, x = x_local_px, y = y_local_px),
                   color = "black", linewidth = 0.1) +
      theme_bw() +
      labs(title = paste0("Slide ", .y$tissue, " - FOV ", .y$fov)) +
      tune::coord_obs_pred() + 
      guides(fill = guide_legend(title = "InSituType Assignments\n(Manley Collapsed) Predicted",ncol = 1)) + 
      scale_fill_manual(values = all_colors)
  })


# #save the plots into PDF form into the figures folder
# pdf("Manley_SMI/results/figures/InSituType_Phenotypes_polygon_plots_predicted-tumor.pdf", width = 14, height = 10)
# tmp = lapply(plots, print)
# dev.off()


# LASSO for gene selection ------------------------------------------------

library(glmnet)
model_dat = plt_dat %>% filter(tumor_class_for_DEGS != "other") %>%
  rownames_to_column("rowname") %>% 
  left_join(manley_data@assays$SCT@scale.data[row.names(pt_diff_genes),] %>%
              as.matrix() %>% t() %>% data.frame(check.names = F) %>% rownames_to_column("rowname"))
model_dat$model_class = ifelse(model_dat$annotation == "Tumor ",0,1)

lasso_dat = model_dat
lasso_predictors = model_dat[,row.names(pt_diff_genes)] %>% as.matrix()
lasso_response = factor(model_dat$tumor_class_for_DEGS, levels = c("tumor_pt", "stroma_pt"))
set.seed(333)
lasso_fit = cv.glmnet(lasso_predictors, lasso_response, family = "binomial", type.measure = "class")
lasso_genes = coef(lasso_fit) %>% as.matrix() %>% data.frame() %>% filter(s1 != 0) %>% row.names() %>% grep("\\(", ., value = T, invert = T)
lasso_dat_all = lasso_dat[,c("model_class", lasso_genes)]
lasso_model_all = glm(model_class ~ ., data = lasso_dat_all, family = "binomial") #0 is tumor, 1 is stroma

lasso_gex = manley_data@assays$SCT@scale.data[lasso_genes,] %>%
  t() %>% data.frame(check.names = F)
#predict all cells from RCC3 1&2 tumor/stroma model
lasso_preds_all = predict(lasso_model_all, lasso_gex, type = "response")

#Add to manley_data, metadata
# lasso_predicted_df = data.frame("rowname" = names(lasso_preds_all),
#                                 Predicted = lasso_preds_all,
#                                 `lasso_genes` = ifelse((manley_data$cell_class == 'tissue' |
#                                                           manley_data$final_insitutype_simplified == "Pelvic epithelium") &
#                                                          lasso_preds_all < 0.5,
#                                                        "Tumor", manley_data$final_insitutype_simplified)) %>%
#   full_join(manley_data@meta.data %>% rownames_to_column("rowname"))
manley_data$lasso_optim = ifelse((manley_data$cell_class == 'tissue' |
                                   manley_data$final_insitutype_simplified == "Pelvic epithelium") &
                                   lasso_preds_all < 0.5,
                                 "Tumor", manley_data$final_insitutype_simplified)
manley_data$lasso_prob = lasso_preds_all
manley_data$lasso_prob025 = ifelse((manley_data$cell_class == 'tissue' |
                                      manley_data$final_insitutype_simplified == "Pelvic epithelium") &
                                     manley_data$lasso_prob < 0.25,
                                   "Tumor", manley_data$final_insitutype_simplified)
#identify the files with the polygon information for all FOVs on all slides
polygons_files = list.files("Manley_SMI/data/SMI-0050_BrandonManley_Moffit/5 Raw data/", 
                            recursive = T,pattern = "polygons", full.names = T)
#now read in the polygon slides
polygon_data = lapply(polygons_files, read.csv, check.names = F)

#for merging need to add which slide (tissue) that the polygon data is for
polygon_data[[1]]$tissue = "RCC3"
polygon_data[[2]]$tissue = "RCC4"
polygon_data[[3]]$tissue = "RCC5"

#concatonate the polygon data into a single data frame and convert the fov and cellID to characters for merging with the metadata
polygon_data = polygon_data %>% 
  do.call(bind_rows, .) %>%
  mutate_at(c("cellID"), as.character)
wkn_tbl = manley_data@meta.data %>% rownames_to_column("rowname")
wkn_tbl$cellID = wkn_tbl$cell_ID %>%
  gsub('.*_', '', .)

#prepping colors for phenotypes
phenotypes = wkn_tbl$lasso_optim %>% unique() %>% sort()
wkn_tbl$lasso_optim = factor(wkn_tbl$lasso_optim,
                             levels = phenotypes[c(1, 3:4, 8, 11:15, 17:22, 24, 29, #immune
                                                   2, 5:7, 9:10, 16, 23, 25:28, 31, #tissue
                                                   30)]) #tumor
wkn_tbl$lasso_prob025 = factor(wkn_tbl$lasso_prob025,
                             levels = phenotypes[c(1, 3:4, 8, 11:15, 17:22, 24, 29, #immune
                                                   2, 5:7, 9:10, 16, 23, 25:28, 31, #tissue
                                                   30)]) #tumor
#create plots
immune_color_func = colorRampPalette(c("red", "blue"))
immune_colors = immune_color_func(17)
tissue_color_func = colorRampPalette(c("green", "orange"))
tissue_colors = tissue_color_func(13)
all_colors = c(immune_colors, tissue_colors, 'black')
names(all_colors) = levels(wkn_tbl$lasso_optim)

full_data = wkn_tbl %>%
  left_join(polygon_data) 
#now we can use the group_map to output a list of plots for all of the slide/fov combinations
plots = full_data %>%
  group_by(tissue, fov) %>%
  group_map(~{
    .x %>%
      ggplot() + 
      geom_polygon(aes(group = cellID, fill = lasso_optim, x = x_local_px, y = y_local_px),
                   color = "black", linewidth = 0.1) +
      theme_bw() +
      labs(title = paste0("Slide ", .y$tissue, " - FOV ", .y$fov),
           caption = "Kidney tissue cells with 'pelvic epithelium'.") +
      tune::coord_obs_pred() + 
      guides(fill = guide_legend(title = "InSituType Assignments\n(Manley Collapsed)",ncol = 1)) + 
      scale_fill_manual(values = all_colors)
  })

plots2 = full_data %>%
  group_by(tissue, fov) %>%
  group_map(~{
    .x %>%
      ggplot() + 
      geom_polygon(aes(group = cellID, fill = lasso_prob025, x = x_local_px, y = y_local_px),
                   color = "black", linewidth = 0.1) +
      theme_bw() +
      labs(title = paste0("Slide ", .y$tissue, " - FOV ", .y$fov, "(LASSO Threshold = 0.25)"),
           caption = "Kidney tissue cells with 'pelvic epithelium'.") +
      tune::coord_obs_pred() + 
      guides(fill = guide_legend(title = "InSituType Assignments\n(Manley Collapsed)",ncol = 1)) + 
      scale_fill_manual(values = all_colors)
  })

plots3 = full_data %>%
  group_by(tissue, fov) %>%
  group_map(~{
    .x %>%
      ggplot() + 
      geom_polygon(aes(group = cellID, fill = lasso_prob, x = x_local_px, y = y_local_px),
                   color = "black", linewidth = 0.1) +
      theme_bw() +
      labs(title = paste0("Slide ", .y$tissue, " - FOV ", .y$fov, "(LASSO Stroma Probabilities)"),
           caption = "Kidney tissue cells with 'pelvic epithelium'.") +
      tune::coord_obs_pred() + 
      guides(fill = guide_legend(title = "LASSO Stroma\nProbability",ncol = 1))
  })

# #save the plots into PDF form into the figures folder
# pdf("Manley_SMI/results/figures/InSituType_Phenotypes_polygon_plots_lasso_predicted-tumor_fixedFOV.pdf", width = 14, height = 10)
# tmp = lapply(plots, print)
# dev.off()
# pdf("Manley_SMI/results/figures/InSituType_Phenotypes_polygon_plots_lasso025_predicted-tumor_fixedFOV.pdf", width = 14, height = 10)
# tmp = lapply(plots2, print)
# dev.off()
# pdf("Manley_SMI/results/figures/InSituType_Phenotypes_polygon_plots_lasso_probabilities_fixedFOV.pdf", width = 14, height = 10)
# tmp = lapply(plots3, print)
# dev.off()


# Increasing Threshold ----------------------------------------------------
gex2 = manley_data@assays$SCT@scale.data %>% 
  as.matrix() %>% t() %>%
  data.frame(check.names = F) %>%
  select(!!pt_genes_up, -VIM)
up2 = apply(gex2 > 0.5, 1, function(x){
  any(x)
})
manley_data$pt10genes_0.5thresh = ifelse(wkn_tbl$cell_class == "tissue" &
                                           up2 == TRUE, "Tumor", manley_data$final_insitutype_simplified)

wkn_tbl = manley_data@meta.data
wkn_tbl$cellID = wkn_tbl$cell_ID %>%
  gsub('.*_', '', .)

#prepping colors for phenotypes
phenotypes = wkn_tbl$pt10genes_0.5thresh %>% unique() %>% sort()
wkn_tbl$pt10genes_0.5thresh = factor(wkn_tbl$pt10genes_0.5thresh,
                             levels = phenotypes[c(1, 3:4, 8, 11:15, 17:22, 24, 29, #immune
                                                   2, 5:7, 9:10, 16, 23, 25:28, 31, #tissue
                                                   30)]) #tumor
#create plots
immune_color_func = colorRampPalette(c("red", "blue"))
immune_colors = immune_color_func(17)
tissue_color_func = colorRampPalette(c("green", "orange"))
tissue_colors = tissue_color_func(13)
all_colors = c(immune_colors, tissue_colors, 'black')
names(all_colors) = levels(wkn_tbl$pt10genes_0.5thresh)

full_data = wkn_tbl %>%
  left_join(polygon_data) 
#now we can use the group_map to output a list of plots for all of the slide/fov combinations
plots = full_data %>%
  group_by(tissue, fov) %>%
  group_map(~{
    .x %>%
      ggplot() + 
      geom_polygon(aes(group = cellID, fill = pt10genes_0.5thresh, x = x_local_px, y = y_local_px),
                   color = "black", linewidth = 0.1) +
      theme_bw() +
      labs(title = paste0("Slide ", .y$tissue, " - FOV ", .y$fov)) +
      tune::coord_obs_pred() + 
      guides(fill = guide_legend(title = "InSituType Assignments\n(Manley Collapsed)",ncol = 1)) + 
      scale_fill_manual(values = all_colors)
  })

# #save the plots into PDF form into the figures folder
pdf("Manley_SMI/results/figures/InSituType_Phenotypes_polygon_plots_ptgenes0.5-tumor_fixedFOV.pdf", width = 14, height = 10)
tmp = lapply(plots, print)
dev.off()

saveRDS(manley_data, "Manley_SMI/results/reports/modeling_tumor/manley_data_tumor_fixedFOV.rds")

# Scatter of Tumor Percent By Sample Group --------------------------------
#hand model vs threshold at 0
#cleaning environment of big objects
#rm(polygon_data, predicted_df, meta_data, full_data, wkn_tbl)
#import clinical and create column that matches that of the metadata
clinical = read.csv("Manley_SMI/data/manley_files/Clinical.data_acs_02.14.23.csv") %>%
  mutate(slide_fov = paste0("slide", Slide, "_fov", FOV),
         Histology = gsub(" $", "", Histology),
         Sarcomatoid = gsub(" $", "", Sarcomatoid),
         Site = gsub(" $", "", Site))
#extract metadata to object
meta_data_updated = manley_data@meta.data %>%
  mutate(slide_fov = paste0("slide", gsub("RCC", "", tissue), "_fov", fov))
#merge clinical and metadata
meta_clinical = right_join(clinical, 
                           meta_data_updated)
treament_naive = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(`Treatment Naive`, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = `Treatment Naive`, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold vs Modeling Pre-Treatment") +
  scale_size_manual(values = c(3,6))

sarcomatoid = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = Sarcomatoid, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold vs Modeling Sarcomatoid") +
  scale_size_manual(values = c(3,6))

metastatic = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold vs Modeling Metastatic") +
  scale_size_manual(values = c(3,6))

fov_callouts = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site, slide_fov) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  ggrepel::geom_text_repel(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`,
                               label = slide_fov), box.padding = 0.5, max.overlaps = 20) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold vs Modeling Names") +
  scale_size_manual(values = c(3,6)) 

pdf("Manley_SMI/results/figures/classifying_tumor/Thresholding_v_Modeling_TumorPercent.pdf", width = 12, height = 10)
treament_naive
sarcomatoid
metastatic
fov_callouts
dev.off()

#hand model vs lasso model
#extract metadata to object
meta_data_updated = manley_data@meta.data %>%
  mutate(slide_fov = paste0("slide", gsub("RCC", "", tissue), "_fov", fov))
#merge clinical and metadata
meta_clinical = right_join(clinical, 
                           meta_data_updated)
treament_naive_lm = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(`Treatment Naive`, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - LASSO` = sum(lasso_optim == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - LASSO` = `Tumor Cells - LASSO` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - LASSO`, y = `Tumor Percent - Model`, 
                 color = `Treatment Naive`, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "LASSO vs Modeling Pre-Treatment") +
  scale_size_manual(values = c(3,6))

sarcomatoid_lm = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - LASSO` = sum(lasso_optim == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - LASSO` = `Tumor Cells - LASSO` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - LASSO`, y = `Tumor Percent - Model`, 
                 color = Sarcomatoid, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "LASSO vs Modeling Sarcomatoid") +
  scale_size_manual(values = c(3,6))

metastatic_lm = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - LASSO` = sum(lasso_optim == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - LASSO` = `Tumor Cells - LASSO` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - LASSO`, y = `Tumor Percent - Model`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "LASSO vs Modeling Metastatic") +
  scale_size_manual(values = c(3,6))

fov_callouts_lm = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site, slide_fov) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - LASSO` = sum(lasso_optim == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - LASSO` = `Tumor Cells - LASSO` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - LASSO`, y = `Tumor Percent - Model`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  ggrepel::geom_text_repel(aes(x = `Tumor Percent - LASSO`, y = `Tumor Percent - Model`,
                               label = slide_fov), box.padding = 0.5, max.overlaps = 20) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "LASSO vs Modeling Names") +
  scale_size_manual(values = c(3,6)) 

pdf("Manley_SMI/results/figures/classifying_tumor/LASSO_v_Modeling_TumorPercent.pdf", width = 12, height = 10)
treament_naive_lm
sarcomatoid_lm
metastatic_lm
fov_callouts_lm
dev.off()

#extract metadata to object
meta_data_updated = manley_data@meta.data %>%
  mutate(slide_fov = paste0("slide", gsub("RCC", "", tissue), "_fov", fov))
#merge clinical and metadata
meta_clinical = right_join(clinical, 
                           meta_data_updated)
treament_naive_mdls = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(`Treatment Naive`, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Cells - LASSO` = sum(lasso_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - LASSO` = `Tumor Cells - LASSO` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - LASSO`, 
                 color = `Treatment Naive`, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold vs LASSO Pre-Treatment") +
  scale_size_manual(values = c(3,6))

sarcomatoid_mdls = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Cells - LASSO` = sum(lasso_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - LASSO` = `Tumor Cells - LASSO` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - LASSO`, 
                 color = Sarcomatoid, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold vs LASSO Sarcomatoid") +
  scale_size_manual(values = c(3,6))

metastatic_mdls = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Cells - LASSO` = sum(lasso_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - LASSO` = `Tumor Cells - LASSO` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - LASSO`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold vs LASSO Metastatic") +
  scale_size_manual(values = c(3,6))

fov_callouts_mdls = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site, slide_fov) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Cells - LASSO` = sum(lasso_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - LASSO` = `Tumor Cells - LASSO` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - LASSO`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  ggrepel::geom_text_repel(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - LASSO`,
                               label = slide_fov), box.padding = 0.5, max.overlaps = 20) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold vs LASSO Names") +
  scale_size_manual(values = c(3,6)) 

pdf("Manley_SMI/results/figures/classifying_tumor/Thresholding_v_LASSO_TumorPercent.pdf", width = 12, height = 10)
treament_naive_mdls
sarcomatoid_mdls
metastatic_mdls
fov_callouts_mdls
dev.off()


# Scatter with 0.5 Threshold ----------------------------------------------
#extract metadata to object
meta_data_updated = manley_data@meta.data %>%
  mutate(slide_fov = paste0("slide", gsub("RCC", "", tissue), "_fov", fov))
#merge clinical and metadata
meta_clinical = right_join(clinical, 
                           meta_data_updated)
#threshold 0.5 with manual model
treament_naive05_ml = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(`Treatment Naive`, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = `Treatment Naive`, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs Modeling Pre-Treatment") +
  scale_size_manual(values = c(3,6))

sarcomatoid05_ml = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = Sarcomatoid, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs Modeling Sarcomatoid") +
  scale_size_manual(values = c(3,6))

metastatic05_ml = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs Modeling Metastatic") +
  scale_size_manual(values = c(3,6))

fov_callouts05_ml = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site, slide_fov) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Model` = sum(pt_manualAIC_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  ggrepel::geom_text_repel(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`,
                               label = slide_fov), box.padding = 0.5, max.overlaps = 20) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs Modeling Names") +
  scale_size_manual(values = c(3,6)) 

pdf("Manley_SMI/results/figures/classifying_tumor/Thresholding0.5_v_Modeling_TumorPercent.pdf", width = 12, height = 10)
treament_naive05_ml
sarcomatoid05_ml
metastatic05_ml
fov_callouts05_ml
dev.off()

#threshold 0.5 with manual model
treament_naive05_lasso = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(`Treatment Naive`, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Model` = sum(lasso_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = `Treatment Naive`, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs LASSO Pre-Treatment") +
  scale_size_manual(values = c(3,6))

sarcomatoid05_lasso = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Model` = sum(lasso_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = Sarcomatoid, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs LASSO Sarcomatoid") +
  scale_size_manual(values = c(3,6))

metastatic05_lasso = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Model` = sum(lasso_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs LASSO Metastatic") +
  scale_size_manual(values = c(3,6))

fov_callouts05_lasso = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site, slide_fov) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Model` = sum(lasso_optim == "Tumor"),
            `Tumor Percent - Threshold` = `Tumor Cells - Threshold` / `Total Cells` * 100,
            `Tumor Percent - Model` = `Tumor Cells - Model` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  ggrepel::geom_text_repel(aes(x = `Tumor Percent - Threshold`, y = `Tumor Percent - Model`,
                               label = slide_fov), box.padding = 0.5, max.overlaps = 20) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs LASSO Names") +
  scale_size_manual(values = c(3,6)) 

pdf("Manley_SMI/results/figures/classifying_tumor/Thresholding0.5_v_LASSO_TumorPercent.pdf", width = 12, height = 10)
treament_naive05_lasso
sarcomatoid05_lasso
metastatic05_lasso
fov_callouts05_lasso
dev.off()

#threshold 0.5 with threshold 0
treament_naive05_0 = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(`Treatment Naive`, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold 0.5` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Threshold 0` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Percent - Threshold 0.5` = `Tumor Cells - Threshold 0.5` / `Total Cells` * 100,
            `Tumor Percent - Threshold 0` = `Tumor Cells - Threshold 0` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold 0.5`, y = `Tumor Percent - Threshold 0`, 
                 color = `Treatment Naive`, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs Threshold 0 Pre-Treatment") +
  scale_size_manual(values = c(3,6))

sarcomatoid05_0 = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold 0.5` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Threshold 0` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Percent - Threshold 0.5` = `Tumor Cells - Threshold 0.5` / `Total Cells` * 100,
            `Tumor Percent - Threshold 0` = `Tumor Cells - Threshold 0` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold 0.5`, y = `Tumor Percent - Threshold 0`, 
                 color = Sarcomatoid, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs Threshold 0 Sarcomatoid") +
  scale_size_manual(values = c(3,6))

metastatic05_0 = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold 0.5` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Threshold 0` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Percent - Threshold 0.5` = `Tumor Cells - Threshold 0.5` / `Total Cells` * 100,
            `Tumor Percent - Threshold 0` = `Tumor Cells - Threshold 0` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold 0.5`, y = `Tumor Percent - Threshold 0`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs Threshold 0 Metastatic") +
  scale_size_manual(values = c(3,6))

fov_callouts05_0 = meta_clinical %>%
  mutate(`Treatment Naive` = ifelse(IT.Treatment.before.collection == "None", "Yes", "No")) %>%
  group_by(Sarcomatoid, Slide, FOV, annotation, Site, slide_fov) %>%
  summarise(Histology = unique(Histology),
            `Total Cells` = n(),
            `Tumor Cells - Threshold 0.5` = sum(pt10genes_0.5thresh == "Tumor"),
            `Tumor Cells - Threshold 0` = sum(pt10genes_0thresh == "Tumor"),
            `Tumor Percent - Threshold 0.5` = `Tumor Cells - Threshold 0.5` / `Total Cells` * 100,
            `Tumor Percent - Threshold 0` = `Tumor Cells - Threshold 0` / `Total Cells` * 100) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = `Tumor Percent - Threshold 0.5`, y = `Tumor Percent - Threshold 0`, 
                 color = Site, shape = Histology, size = annotation), alpha = 0.5) +
  ggrepel::geom_text_repel(aes(x = `Tumor Percent - Threshold 0.5`, y = `Tumor Percent - Threshold 0`,
                               label = slide_fov), box.padding = 0.5, max.overlaps = 20) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Threshold 0.5 vs Threshold 0 Names") +
  scale_size_manual(values = c(3,6)) 

pdf("Manley_SMI/results/figures/classifying_tumor/Thresholding0.5_v_thresholding0_TumorPercent.pdf", width = 12, height = 10)
treament_naive05_0
sarcomatoid05_0
metastatic05_0
fov_callouts05_0
dev.off()

# Testing methods ---------------------------------------------------------

# preds_RCC3_1_2 = predict(logit_mod_RCC3_1_2, test_dat_RCC3_1_2, type = "response")
# preds_scaled_RCC3_1_2 = ifelse(preds_RCC3_1_2 < 0.5, 0, 1)
# con_tab_RCC3_1_2 = table(test_dat_RCC3_1_2$fov, preds_scaled_RCC3_1_2)
# con_tab_RCC3_1_2
# sum(con_tab_RCC3_1_2[c(1,4)])/sum(con_tab_RCC3_1_2)
# 
# 
# # RCC3 FOV 21 and 22 ------------------------------------------------------
# 
# proximal_tubules_RCC3_21_22 = subset(x = slide3, cells = row.names(slide3@meta.data %>% filter(fov %in% c(21, 22), 
#                                                                                                Slide_name == "RCC3",
#                                                                                                final_insitutype_simplified == "Proximal tubule")))
# 
# proximal_tubules_RCC3_21_22 = RunPCA(object = proximal_tubules_RCC3_21_22, assay = "SCT", npcs = 20, reduction.name = "pca_focus")
# proximal_tubules_RCC3_21_22 = RunUMAP(object = proximal_tubules_RCC3_21_22, reduction = "pca_focus", dims = 1:20, assay = "SCT",
#                                       n.components = 3, reduction.name = "SCT_umap_focus", repulsion.strength = 10)
# 
# DimPlot(proximal_tubules_RCC3_21_22, reduction = "SCT_umap_focus", group.by = "fov",
#         #cells.highlight = list(WhichCells(tumors, idents = c("Tumor"))),
#         repel = T,
#         raster.dpi = c(512, 512)) +
#   #theme(legend.position = "none") +
#   labs(title = "RCC3-21 and RCC3-22 Proximal Tubule Cells")
# 
# Idents(proximal_tubules_RCC3_21_22) = "fov"
# degs_RCC3_21_22 = FindAllMarkers(proximal_tubules_RCC3_21_22)
# model_genes_RCC3_21_22 = degs_RCC3_21_22 %>% filter(avg_log2FC > 1) %>% pull(gene)
# 
# #DEGs common between 1-2 and 21-22
# intersect(model_genes_RCC3_1_2, model_genes_RCC3_21_22)
# 
# model_dat_RCC3_21_22 = full_join(proximal_tubules_RCC3_21_22@meta.data %>% 
#                                    select(fov) %>%
#                                    rownames_to_column("slide_fov_cell"),
#                                  proximal_tubules_RCC3_21_22@assays$SCT@data[model_genes_RCC3_21_22,] %>%
#                                    as.matrix() %>% t() %>% data.frame(check.names=F) %>% 
#                                    rownames_to_column("slide_fov_cell")) %>% 
#   column_to_rownames("slide_fov_cell")
# 
# model_dat_RCC3_21_22 %>%
#   mutate(fov = as.factor(fov)) %>%
#   gather("gene", "expression", -fov) %>%
#   ggplot() + 
#   geom_density(aes(x = expression, color = fov, group = fov)) +
#   facet_wrap(~gene)
# 
# model_dat_RCC3_21_22$fov = model_dat_RCC3_21_22$fov - 21 #remember 2 is now 1 and 1 is now 0
# tt_split_RCC3_21_22 = sample(c(1, 2), nrow(model_dat_RCC3_21_22), prob = c(0.2, 0.8), replace = T)
# train_dat_RCC3_21_22 = model_dat_RCC3_21_22[tt_split_RCC3_21_22 == 2,]
# test_dat_RCC3_21_22 = model_dat_RCC3_21_22[tt_split_RCC3_21_22 == 1,]
# logit_mod_RCC3_21_22 = glm(fov ~ ., data = train_dat_RCC3_21_22, family = "binomial")
# summary(logit_mod_RCC3_21_22)
# preds_RCC3_21_22 = predict(logit_mod_RCC3_21_22, test_dat_RCC3_21_22, type = "response")
# preds_scaled_RCC3_21_22 = ifelse(preds_RCC3_21_22 < 0.5, 0, 1)
# con_tab_RCC3_21_22 = table(test_dat_RCC3_21_22$fov, preds_scaled_RCC3_21_22)
# con_tab_RCC3_21_22
# sum(con_tab_RCC3_21_22[c(1,4)])/sum(con_tab_RCC3_21_22)
# 
# 
# # Exploring prediction on all tissue --------------------------------------
# 
# 
# tissue_RCC3_21_22 = subset(x = slide3, cells = row.names(slide3@meta.data %>% filter(fov %in% c(21, 22), 
#                                                                                      Slide_name == "RCC3",
#                                                                                      cell_class == "tissue")))
# 
# model_dat_RCC3_21_22_tissue = full_join(tissue_RCC3_21_22@meta.data %>% 
#                                           select(fov, final_insitutype_simplified) %>%
#                                           rownames_to_column("slide_fov_cell"),
#                                         tissue_RCC3_21_22@assays$SCT@data[model_genes_RCC3_21_22,] %>%
#                                           as.matrix() %>% t() %>% data.frame(check.names=F) %>% 
#                                           rownames_to_column("slide_fov_cell")) %>% 
#   column_to_rownames("slide_fov_cell")
# sample_21 = model_dat_RCC3_21_22_tissue %>% filter(fov == 21)
# 
# pred_tissue = predict(logit_mod_RCC3_21_22, sample_21 %>% select(-final_insitutype_simplified), type = "response")
# pred_tissue = ifelse(pred_tissue < 0.5, 0 ,1)
# table(pred_tissue, sample_21$final_insitutype_simplified)
# 
# 
# tissue_RCC3_1 = subset(x = slide3, cells = row.names(slide3@meta.data %>% filter(fov %in% c(8), 
#                                                                                  Slide_name == "RCC3",
#                                                                                  cell_class == "tissue")))
# 
# model_dat_RCC3_1_tissue = full_join(tissue_RCC3_1@meta.data %>% 
#                                       select(fov, final_insitutype_simplified) %>%
#                                       rownames_to_column("slide_fov_cell"),
#                                     tissue_RCC3_1@assays$SCT@data[model_genes_RCC3_21_22,] %>%
#                                       as.matrix() %>% t() %>% data.frame(check.names=F) %>% 
#                                       rownames_to_column("slide_fov_cell")) %>% 
#   column_to_rownames("slide_fov_cell")
# sample_1 = model_dat_RCC3_1_tissue %>% filter(fov == 8)
# 
# pred_tissue = predict(logit_mod_RCC3_21_22, sample_1 %>% select(-final_insitutype_simplified), type = "response")
# pred_tissue = ifelse(pred_tissue < 0.5, 0 ,1)
# table(pred_tissue, sample_1$final_insitutype_simplified)
# 
# 
# # Fitting Model on Proximal Tubules ---------------------------------------
# 
# clinical = read.csv("Manley_SMI/data/manley_files/Clinical.data_acs_02.14.23.csv") %>%
#   mutate(slide_fov = paste0("slide", Slide, "_fov", FOV))
# modeling_data = subset(manley_data, subset = "slide_fov" %in% clinical %>% 
#                          filter(IT.Treatment.before.collection == "None") %>% 
#                          pull(slide_fov))
# picked_genes = intersect(model_genes_RCC3_1_2, model_genes_RCC3_21_22)
# modeling_data = left_join(manley_data@meta.data %>%
#                             filter(slide_fov %in% unique(clinical %>% 
#                                                            filter(IT.Treatment.before.collection == "None") %>% 
#                                                            pull(slide_fov)),
#                                    final_insitutype_simplified == "Proximal tubule") %>% 
#                             select(fov, final_insitutype_simplified, annotation) %>%
#                             rownames_to_column("slide_fov_cell"),
#                           manley_data@assays$SCT@data[picked_genes,] %>%
#                             as.matrix() %>% t() %>% data.frame(check.names=F) %>% 
#                             rownames_to_column("slide_fov_cell")) %>% 
#   column_to_rownames("slide_fov_cell")
# 
# modeling_data$annotation = factor(modeling_data$annotation)
# logit_mod_tumor = glm(annotation ~ ., data = modeling_data %>% 
#                         select(-fov, -final_insitutype_simplified), family = "binomial")
# summary(logit_mod_tumor)
# all_cells_pred = predict(logit_mod_tumor, manley_data@assays$SCT@data %>% 
#                            as.matrix() %>% t() %>% data.frame(check.names = F) %>%
#                            select(!!picked_genes), type = "response")
# all_cells_pred = predict(logit_mod_RCC3_21_22, manley_data@assays$SCT@data %>% 
#                            as.matrix() %>% t() %>% data.frame(check.names = F) %>%
#                            select(!!model_genes_RCC3_21_22), type = "response")
# all_cells_pred = predict(logit_mod_RCC3_21_22, manley_data@assays$SCT@data %>% 
#                            as.matrix() %>% t() %>% data.frame(check.names = F) %>%
#                            select(!!model_genes_RCC3_21_22) %>%
#                            rownames_to_column("slide_fov_cell") %>%
#                            filter(slide_fov_cell %in% (manley_data@meta.data %>%
#                                                          rownames_to_column("slide_fov_cell") %>%
#                                                          filter(cell_class == "tissue") %>% pull(slide_fov_cell))), type = "response")
# ggplot() +
#   geom_histogram(aes(x = all_cells_pred), binwidth = 0.01)
