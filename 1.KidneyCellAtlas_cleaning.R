library(SCP)
library(reticulate)
library(tidyverse)

#
sc <- import("scanpy")
adata <- sc$read_h5ad("Manley_SMI/data/KidneyCellAtlas/Mature_Full_v3.h5ad")
srt = adata_to_srt(adata)
saveRDS(srt, "Manley_SMI/data/KidneyCellAtlas/Mature_Full_v3_seurat.rds")
#restart r session
srt_sce = as.SingleCellExperiment(srt)

#prepping for creating summarised experiment
KCA_expression = srt@assays$RNA@data
KCA_meta = srt@meta.data

saveRDS(KCA_expression, "Manley_SMI/data/KidneyCellAtlas/Mature_Full_correctedCounts.rds")
saveRDS(KCA_meta, "Manley_SMI/data/KidneyCellAtlas/Mature_Full_v3_metadata.rds")


# Exploring and Cleaning --------------------------------------------------

rm(list=ls())

KCA_expression = readRDS("Manley_SMI/data/KidneyCellAtlas/Mature_Full_correctedCounts.rds")
KCA_meta = readRDS("Manley_SMI/data/KidneyCellAtlas/Mature_Full_v3_metadata.rds")

# Creating seurat object --------------------------------------------------

KCA_expression = readRDS("Manley_SMI/data/KidneyCellAtlas/Mature_Full_correctedCounts.rds")
KCA_meta = readRDS("Manley_SMI/data/KidneyCellAtlas/Mature_Full_v3_metadata.rds")

seu = CreateSeuratObject(counts = KCA_expression, meta.data = KCA_meta)

seu = SCTransform(seu) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

saveRDS(seu, "Manley_SMI/data/KidneyCellAtlas/KidneyCellAtlas_seurat.rds")




