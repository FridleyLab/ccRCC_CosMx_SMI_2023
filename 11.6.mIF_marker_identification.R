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
  filter(!(unique_fov %in% c("RCC5_19", "RCC5_20"))) %>%
  mutate(Sarcomatoid = gsub(" ", "", Sarcomatoid),
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

morans_list = readRDS("Manley_SMI/results/ligand_receptor/Manley_ligand_receptor_CellTalkDB-knn3_EMT.rds")

#Making data frame of data
morans_df = do.call(bind_rows, morans_list) %>%
  spread("lr_pair", "original")
lr_pairs=EMT_lr_pairs$lr_pair

#Associations with primary tumor treatment stroma
primary_prepost_stroma_meta = metadata %>%
  filter(Site == "Tumor", Source == "Stroma") %>%
  filter(!(Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive")) %>%
  select(unique_fov, Pretreatment.IO, Sarcomatoid) %>% distinct()

lr_pair = "COL4A1_ITGAV"
df = primary_prepost_stroma_meta %>%
  left_join(morans_df %>% select(unique_fov, !!lr_pair), by = join_by(unique_fov))


# Spatstat Smooth ---------------------------------------------------------

meta_4_8 = metadata %>%
  filter(unique_fov == "RCC4_8")
gene_4_8 = scaled_expression[,(lr_pair %>% str_split("_") %>% unlist())] %>%
  data.frame(check.names=F) %>%
  rownames_to_column("id")
data_4_8 = left_join(meta_4_8, gene_4_8)

win = convexhull.xy(data_4_8$CenterX_local_px,
                    data_4_8$CenterY_local_px)
COL4A1.ppp = ppp(data_4_8$CenterX_local_px,
                 data_4_8$CenterY_local_px,
                 window = win, 
                 marks = data_4_8$COL4A1)
COL4A1.smooth = Smooth(COL4A1.ppp, sigma = 50)
ITGAV.ppp = ppp(data_4_8$CenterX_local_px,
                data_4_8$CenterY_local_px,
                window = win, 
                marks = data_4_8$ITGAV)
ITGAV.smooth = Smooth(ITGAV.ppp, sigma = 50)

diff.smooth = COL4A1.smooth
diff.smooth$v = ((COL4A1.smooth$v + min(COL4A1.smooth, na.rm = T)) / max(COL4A1.smooth$v, na.rm = T)) -
  ((ITGAV.smooth$v + min(ITGAV.smooth, na.rm = T)) / max(ITGAV.smooth$v, na.rm = T))

# diff.smooth$v[diff.smooth$v < -1] = NA
# diff.smooth$v[diff.smooth$v > 1] = NA

pdf("Manley_SMI/results/figures/ligand_receptor/RCC4-8_COL4A1-ITGAV_smooth_MI0.177.pdf", height = 15, width = 5)
par(mfrow = c(3, 1))
plot(COL4A1.smooth, main = "COL4A1 Expression")
plot(ITGAV.smooth, main = "ITGAV Expression")
plot(diff.smooth, main = "COL4A1 - ITGAV Expression")
dev.off()

#negative

meta_3_6 = metadata %>%
  filter(unique_fov == "RCC3_6")
gene_3_6 = scaled_expression[,(lr_pair %>% str_split("_") %>% unlist())] %>%
  data.frame(check.names=F) %>%
  rownames_to_column("id")
data_3_6 = left_join(meta_3_6, gene_3_6)

win = convexhull.xy(data_3_6$CenterX_local_px,
                    data_3_6$CenterY_local_px)
COL4A1.ppp = ppp(data_3_6$CenterX_local_px,
                 data_3_6$CenterY_local_px,
                 window = win, 
                 marks = data_3_6$COL4A1)
COL4A1.smooth = Smooth(COL4A1.ppp, sigma = 50)
ITGAV.ppp = ppp(data_3_6$CenterX_local_px,
                data_3_6$CenterY_local_px,
                window = win, 
                marks = data_3_6$ITGAV)
ITGAV.smooth = Smooth(ITGAV.ppp, sigma = 50)

diff.smooth = COL4A1.smooth
diff.smooth$v = ((COL4A1.smooth$v + min(COL4A1.smooth, na.rm = T)) / max(COL4A1.smooth$v, na.rm = T)) -
  ((ITGAV.smooth$v + min(ITGAV.smooth, na.rm = T)) / max(ITGAV.smooth$v, na.rm = T))
# diff.smooth$v[diff.smooth$v < -1] = NA
# diff.smooth$v[diff.smooth$v > 1] = NA

pdf("Manley_SMI/results/figures/ligand_receptor/RCC3-6_COL4A1-ITGAV_smooth_MI-0.026.pdf", height = 15, width = 5)
par(mfrow = c(3, 1))
plot(COL4A1.smooth, main = "COL4A1 Expression")
plot(ITGAV.smooth, main = "ITGAV Expression")
plot(diff.smooth, main = "COL4A1 - ITGAV Expression")
dev.off()

# Plotting RCC4_8 - highest Moran's I -------------------------------------


# RCC4_4
# RCC4_2
# RCC4_18

meta_4_8 = metadata %>%
  filter(unique_fov == "RCC4_8")
gene_4_8 = scaled_expression[,(lr_pair %>% str_split("_") %>% unlist())] %>%
  data.frame(check.names=F) %>%
  rownames_to_column("id")
data_4_8 = left_join(meta_4_8, gene_4_8)

gene_boxplot = data_4_8 %>%
  select(lasso_final, COL4A1, ITGAV) %>%
  gather("Gene", "SCT Expression", -lasso_final) %>%
  ggplot() +
  geom_boxplot(aes(x = lasso_final, y = `SCT Expression`, color = Gene)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(title = "RCC4 - FOV8 Expression of EMT Ligand-Receptor")

pdf("Manley_SMI/results/ligand_receptor/RCC4_8_COL4A1-ITGAV_phenotype_expression.pdf", height = 7, width = 10)
gene_boxplot
dev.off()

meta_4_2 = metadata %>%
  filter(unique_fov == "RCC4_2")
gene_4_2 = scaled_expression[,(lr_pair %>% str_split("_") %>% unlist())] %>%
  data.frame(check.names=F) %>%
  rownames_to_column("id")
data_4_2 = left_join(meta_4_2, gene_4_2)

gene_boxplot = data_4_2 %>%
  select(lasso_final, COL4A1, ITGAV) %>%
  gather("Gene", "SCT Expression", -lasso_final) %>%
  ggplot() +
  geom_boxplot(aes(x = lasso_final, y = `SCT Expression`, color = Gene)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(title = "RCC4 - FOV8 Expression of EMT Ligand-Receptor")

pdf("Manley_SMI/results/ligand_receptor/RCC4_2_COL4A1-ITGAV_phenotype_expression.pdf", height = 7, width = 10)
gene_boxplot
dev.off()


# Polygon plots -----------------------------------------------------------

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

full_dat = right_join(polygon_data, data_4_8 %>%
                        mutate(cellID = gsub('.*_', '', cell_ID))
)

COL4A1_Polygon = full_dat %>%
  filter(lasso_final %in% c("Capillary endothelium", "Epithelial progenitor cell", "Fibroblast", "Myofibroblast", "Tumor", "Vasa recta endothelium")) %>%
  ggplot() + 
  # geom_polygon(data = . %>% filter(lasso_final != "Epithelial progenitor cell"),
  #              aes(group = cellID, fill = COL4A1,  x = x_local_px, y = y_local_px),
  #              color = "black", size = 0.1) +
  geom_polygon(data = . %>% filter(lasso_final == "Epithelial progenitor cell"),
               aes(group = cellID, fill = ITGAV,  x = x_local_px, y = y_local_px),
               color = "black", size = 0.1) +
  tune::coord_obs_pred() +
  theme_bw() +
  labs(title = "RCC5-1, Significant Allograft Rejection Clustering") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2()


# Kriging -----------------------------------------------------------------

library(sp)
library(gstat)
#prep data
data_4_8 = left_join(meta_4_8, gene_4_8)
fov_grid = expand.grid(CenterX_local_px = seq(min(data_4_8$CenterX_local_px), max(data_4_8$CenterX_local_px), 100),
                       CenterY_local_px = seq(min(data_4_8$CenterY_local_px), max(data_4_8$CenterY_local_px), 100))
coordinates(fov_grid) = ~ CenterX_local_px + CenterY_local_px

coordinates(data_4_8) = ~ CenterX_local_px + CenterY_local_px

#variogram?
COL4A1.vgm = variogram(COL4A1 ~ 1, data_4_8)
COL4A1.fit = fit.variogram(COL4A1.vgm, model = vgm(1, "Sph", 900, 1))
plot(COL4A1.vgm, COL4A1.fit)

COL4A1.kriged = krige(COL4A1 ~ 1, data_4_8, fov_grid, COL4A1.fit)

COL4A1.kriged %>%
  data.frame(check.names = F) %>%
  ggplot(aes(x = CenterX_local_px, y = CenterY_local_px)) +
  geom_tile(aes(fill = var1.pred)) +
  scale_fill_gradientn(colors = c("lightpink", "red4"), breaks = c(-0.5, 7.5))

#variogram?
ITGAV.vgm = variogram(ITGAV ~ 1, data_4_8)
ITGAV.fit = fit.variogram(ITGAV.vgm, model = vgm(1, "Sph", 900, 1))
plot(ITGAV.vgm, ITGAV.fit)

ITGAV.kriged = krige(ITGAV ~ 1, data_4_8, fov_grid, ITGAV.fit)

ITGAV.kriged %>%
  data.frame(check.names = F) %>%
  ggplot(aes(x = CenterX_local_px, y = CenterY_local_px)) +
  geom_tile(aes(fill = var1.pred)) +
  scale_fill_gradientn(colors = c("lightyellow", "yellow4"), breaks = c(-0.2, 0.2)) +
  theme(legend.position = "right")

pl = ggpubr::ggarrange(COL4A1.kriged %>%
                         data.frame(check.names = F) %>%
                         ggplot(aes(x = CenterX_local_px, y = CenterY_local_px)) +
                         geom_tile(aes(fill = var1.pred)) +
                         tune::coord_obs_pred() +
                         theme_bw() +
                         scale_fill_gradientn(colors = c("lightpink", "red4"), breaks = c(-0.5, 7.5)) +
                         labs(title = "COL4A1"),
                       ITGAV.kriged %>%
                         data.frame(check.names = F) %>%
                         ggplot(aes(x = CenterX_local_px, y = CenterY_local_px)) +
                         geom_tile(aes(fill = var1.pred)) +
                         tune::coord_obs_pred() +
                         theme_bw() +
                         scale_fill_gradientn(colors = c("lightyellow", "yellow4"), breaks = c(-0.2, 0.2)) +
                         labs(title = "ITGAV"),
                       ncol = 1, nrow = 2, common.legend = FALSE)
pdf("Manley_SMI/results/figures/ligand_receptor/EMT_RCC4-8_COL4A1-ITGAV_kriging.pdf")
pl
dev.off()