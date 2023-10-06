library(plotly)
manley_data = readRDS("Manley_SMI/results/rds/final_explore.rds")
dat = manley_data@reductions$SCT.umap3D@cell.embeddings %>% data.frame()
clusts = ifelse(grepl("CD|NK", manley_data$SingleR.labels), "T/NK cells", "Other")
clusts = manley_data$SingleR.labels
clusts[!grepl("CD|NK", manley_data$SingleR.labels)] <- "Other"
fig = plot_ly(dat, x = ~sct.umap3d_1, y = ~sct.umap3d_2, z = ~sct.umap3d_3, 
              color = clusts,
              size = 1)
fig
