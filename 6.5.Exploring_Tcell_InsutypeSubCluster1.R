
# exploring T cell cluster 1 ----------------------------------------------

#library
library(tidyverse)
library(Seurat)

#data
#clustered t cell data from insitutype
t_cell_data = readRDS("Manley_SMI/results/rds/InSituType/T-cell_Louvain_Clustering.rds")
#DEGs from the above clustered T cells
DEGs = readRDS("Manley_SMI/results/rds/InSituType/T-cell_DEGS_inLouvainClusters.rds")

#
meta_data = t_cell_data@meta.data
counts = t_cell_data@assays$Nanostring@counts %>% 
  as.matrix() %>%
  t() %>%
  as.data.frame(check.names = F)

#comparing insitutype cell assignments to seurat_clusters
meta_data %>%
  group_by(seurat_clusters, insitutype) %>%
  summarise(n = n()) %>%
  spread("seurat_clusters", "n") %>% View()
#doesn't apear it is overaly a single cell type
#a little less than half are assigned with CD4 while the other half is CD8, and NKT cells
#468 as NK cells

test_data = full_join(meta_data %>% select(seurat_clusters) %>% rownames_to_column("rowname"),
                      counts %>% rownames_to_column("rowname")) %>%
  mutate(group = ifelse(seurat_clusters == 1, 1, 0))

#doing quick t test to check for differences in unknown v others
formulas = paste0("`",colnames(counts), "` ~ group")
res = sapply(formulas, function(f){
  res = t.test(as.formula(f), data = test_data)
  c(res$estimate, p.value = res$p.value)
}) %>% t()

res = res %>%
  data.frame(check.names = F) %>%
  mutate(p.adj = p.adjust(p.value),
         difference = `mean in group 1` - `mean in group 0`)

#plot FOXP3 count density

test_data %>% 
  ggplot() +
  geom_density(aes(x = FOXP3, color = as.factor(group), y = ..scaled..))
test_data %>% 
  ggplot() +
  geom_density(aes(x = CTLA4, color = as.factor(group), y = ..scaled..))

test_data %>%
  ggplot() + 
  geom_point(aes(x = CTLA4, y = FOXP3, color = as.factor(group)))

