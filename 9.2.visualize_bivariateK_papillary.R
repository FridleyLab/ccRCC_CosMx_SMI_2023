
# Plotting Univar K -------------------------------------------------------
rm(list=ls())
library(spatialTIME)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

manley_mif = readRDS("Manley_SMI/data/bivarK_mif/manley_mif_papillary.rds")
clinical = manley_mif$clinical %>%
  mutate(Histology = gsub(" $", "", Histology),
         Sarcomatoid = gsub(" $", "", Sarcomatoid),
         Gender = gsub(" $", "", Gender),
         Laterality = gsub(" $", "", Laterality),
         Site = gsub(" $", "", Site),
         Response.Group = case_when(Response.to.IT.treatment.after.collection..1st.line. %in% c("Complete Response", "Patial Response") ~ "Response",
                                    Response.to.IT.treatment.after.collection..1st.line. %in% c("stable disease", "Progression") ~ "Disease",
                                    T ~ "Unknown"),
         Response.Group = as.factor(Response.Group),
         TNM.stage = factor(TNM.stage, levels = c(1,3,4))) #%>% filter(Site == "Tumor") # adding this vastly changes the signifi

biK_dat = manley_mif$derived$bivariate_Count %>%
  filter(r == 275) %>%
  full_join(manley_mif$sample, .) %>%
  right_join(clinical, .) %>%
  mutate(Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received Treatment"),
         Pretreatment.IO = factor(Pretreatment.IO, level = c("Treatment Naive", "Received Treatment")))

# Primary pre v post IO ---------------------------------------------------

biK_dat_simp_naive = biK_dat %>%
  filter(Site == "Tumor", fov <= 20) %>%
  group_by(Source, Pretreatment.IO, Anchor, Counted) %>%
  summarise(`Mean DOCE` = mean(`Degree of Clustering Exact`, na.rm=T))
biK_dat_simp_io = biK_dat %>%
  filter(Site == "Tumor", fov <= 20) %>%
  group_by(Source, Pretreatment.IO, Anchor, Counted) %>%
  summarise(`Mean DOCE` = mean(`Degree of Clustering Exact`, na.rm=T))

get_heatmap = function(mat, title, ord){
  top = HeatmapAnnotation(
    immune = ord,
    col = list(immune = c("Immune" = "blue", "Tissue" = "green")),
    border = T, show_annotation_name = F, show_legend = F
  )
  right = rowAnnotation(
    immune = ord,
    col = list(immune = c("Immune" = "blue", "Tissue" = "green")),
    border = T, show_annotation_name = F
  )
  high = max(abs(mat), na.rm =T)
  col_fun = circlize::colorRamp2(c(-high, 0, high), c("blue", "white", "red"))
  ht = Heatmap(mat, show_column_names = F, 
               bottom_annotation = HeatmapAnnotation(
                 text = anno_text(colnames(mat), rot = 45, location = unit(1, "npc"), just = "right"),
                 annotation_height = max_text_width(colnames(mat))
               ),
               heatmap_legend_param = list(
                 title = "Mean DOCE"
               ),
               column_title = title,
               top_annotation = top,
               right_annotation = right,
               col = col_fun
  ) + rowAnnotation(Counted = anno_text(rownames(mat), just = "left", show_name = F, rot = 45))
  out = draw(ht, padding = unit(c(2, 25, 25, 2), "mm"))
  return(out)
}

naive_stroma = biK_dat_simp_naive %>%
  filter(Source == "Stroma", Pretreatment.IO == "Treatment Naive") %>%
  #filter(Anchor != "Possibly Mid-Rep (Misc. Cells)", Counted != "Possibly Mid-Rep (Misc. Cells)") %>%
  ungroup() %>%
  spread("Counted", "Mean DOCE") %>%
  ungroup() %>%
  select(-1, -2) %>%
  column_to_rownames("Anchor") %>%
  as.matrix()
naive_stroma[is.nan(naive_stroma)] = 0
vec = c(1,2,1,1,2,2,2,1,2,2,1,1,1,1,1,2,1,1,1,1,1,1,2,1,2,2,2,2,1,2,2)
names(vec) = colnames(naive_stroma) 
naive_stroma_ht = get_heatmap(naive_stroma[!as.logical(vec-1), !as.logical(vec-1)],#only immune because tissue is unreliable
                              "Treatment Naive Primary Stroma (r = 275, Stroma)",
                              rep("Immune", 17))

io_stroma = biK_dat_simp_io %>%
  filter(Source == "Stroma", Pretreatment.IO == "Received Treatment") %>%
  ungroup() %>%
  spread("Counted", "Mean DOCE") %>%
  ungroup() %>%
  select(-1, -2) %>%
  column_to_rownames("Anchor") %>%
  as.matrix()
io_stroma[is.nan(io_stroma)] = 0
io_stroma = io_stroma[
  !sapply(1:nrow(io_stroma), function(r) all(is.na(io_stroma[r,]))),
  !sapply(1:ncol(io_stroma), function(c) all(is.na(io_stroma[,c])))
]

io_stroma_ht = get_heatmap(io_stroma[!as.logical(vec-1), !as.logical(vec-1)], 
                           "Post IO Primary Stroma (r = 275, Stroma)",
                           rep("Immune", 17))



# Tumor -------------------------------------------------------------------

naive_tumor = biK_dat_simp_naive %>%
  filter(Source == "Tumor", Pretreatment.IO == "Treatment Naive") %>%
  spread("Counted", "Mean DOCE") %>%
  ungroup() %>%
  select(-1, -2) %>%
  column_to_rownames("Anchor") %>%
  as.matrix()
naive_tumor = naive_tumor[
  !sapply(1:nrow(naive_tumor), function(r) all(is.na(naive_tumor[r,]))),
  !sapply(1:ncol(naive_tumor), function(c) all(is.na(naive_tumor[,c])))
]
keep_cols = ifelse(vec[match(colnames(naive_tumor), names(vec))] == "1", TRUE, FALSE)
naive_tumor_ht = get_heatmap(naive_tumor[keep_cols, keep_cols], "Treatment Naive Primary Tumor (r = 275, Tumor)",
                             rep("Immune", 15))

io_tumor = biK_dat_simp_io %>%
  filter(Source == "Tumor", Pretreatment.IO == "Received Treatment") %>%
  spread("Counted", "Mean DOCE") %>%
  ungroup() %>%
  select(-1, -2) %>%
  column_to_rownames("Anchor") %>%
  as.matrix()
io_tumor = io_tumor[
  !sapply(1:nrow(io_tumor), function(r) all(is.na(io_tumor[r,]))),
  !sapply(1:ncol(io_tumor), function(c) all(is.na(io_tumor[,c])))
]
keep_cols = ifelse(vec[match(colnames(io_tumor), names(vec))] == "1", TRUE, FALSE)
io_tumor_ht = get_heatmap(io_tumor[keep_cols, keep_cols], "Post IO Primary Tumor (r = 275, Tumor)",
                          rep("Immune", 16))

pdf("Manley_SMI/results/papillary/treatment-naive-primaryTumor_stroma_biK-mean-r275_papillary.pdf", height = 10, width = 10)
naive_stroma_ht
dev.off()

pdf("Manley_SMI/results/papillary/received-io-primaryTumor_stroma_biK-mean-r275_papillary.pdf", height = 10, width = 10)
io_stroma_ht
dev.off()

pdf("Manley_SMI/results/papillary/treatment-naive-primaryTumor_tumor_biK-mean-r275_papillary.pdf", height = 10, width = 10)
naive_tumor_ht
dev.off()

pdf("Manley_SMI/results/papillary/received-io-primaryTumor_tumor_biK-mean-r275_papillary.pdf", height = 10, width = 10)
io_tumor_ht
dev.off()
