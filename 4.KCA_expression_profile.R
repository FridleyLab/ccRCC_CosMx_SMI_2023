library(tidyverse)
library(Seurat)

KCA = readRDS("Manley_SMI/data/KidneyCellAtlas/KidneyCellAtlas_seurat.rds")
kca_counts = readRDS("Manley_SMI/data/KidneyCellAtlas/KidneyCellAtlas_rawcounts.rds")
kca_meta = KCA@meta.data

libsizes = colSums(kca_counts)
scaling_factor = max(libsizes)/libsizes
kca_counts_scaled = sweep(kca_counts, 2, scaling_factor, '*')

vec = kca_meta$celltype %>% unique()
names(vec) = vec

kca_counts_large = as.matrix(kca_counts_scaled)
kca_means = sapply(vec, function(x){
  rowMeans(kca_counts_large[,which(colnames(kca_counts_large) %in% row.names(kca_meta[kca_meta$celltype == x,]))])
})

saveRDS(kca_means, "Manley_SMI/data/KidneyCellAtlas/KCA_expressionProfile.rds")