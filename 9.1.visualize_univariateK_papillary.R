rm(list=ls())
# Plotting Univar K -------------------------------------------------------

library(spatialTIME)
library(tidyverse)

manley_mif = readRDS("Manley_SMI/data/UnivariateK_mif/manley_mif_papillary.rds")
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

uniK_dat = manley_mif$derived$univariate_Count %>%
  full_join(manley_mif$sample, .) %>%
  left_join(clinical, .)

uniK_dat %>%
  filter(Marker == "CD8 T cell") %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Exact`, group = unique_fov, color = Response.Group)) + 
  facet_grid(~Source)

uniK_dat %>%
  filter(Marker == "CD8 T cell") %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Exact`, group = unique_fov, color = Gender)) + 
  facet_grid(~Source)

uniK_dat %>%
  filter(Marker == "CD8 T cell") %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Exact`, group = unique_fov, color = Laterality)) + 
  facet_grid(~Source)

uniK_dat %>%
  filter(Marker == "CD8 T cell") %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Exact`, group = unique_fov, color = tissue)) + 
  facet_grid(Sarcomatoid~Source)

uniK_dat %>%
  filter(Marker == "CD8 T cell",
         Histology == "clear cell") %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Exact`, group = unique_fov, color = IT.Treatment.after.collection)) + 
  facet_grid(~Source)

uniK_dat %>%
  filter(Marker == "CD8 T cell",
         r == 275) %>%
  ggplot() +
  geom_boxplot(aes(x = Source, y = `Degree of Clustering Exact`, color = IT.Treatment.after.collection))

uniK_dat %>%
  filter(Marker == "CD8 T cell",
         r == 275) %>%
  ggplot() +
  geom_boxplot(aes(x = Source, y = `Degree of Clustering Exact`, color = IT.Treatment.before.collection))

uniK_dat %>%
  filter(Marker == "CD8 T cell",
         r == 275) %>%
  ggplot() +
  geom_boxplot(aes(x = Source, y = `Degree of Clustering Exact`, color = Response.Group))

#####
#diff between tumor and stroma
#remove samples without cores (failed)
#remove samples in FOVs > 20
uniK_dat_diff = manley_mif$derived$univariate_Count %>%
  left_join(clinical, .) %>% 
  filter(fov<=20) %>%
  select(Source, Specimen.ID.primary, Marker, r, `Degree of Clustering Exact`) %>% 
  filter(!is.na(r)) %>%
  #since there are now some patients who have 2 tumor or 2 stroma, 
  #have to take the mean in order to not have the spread blowup
  group_by(Source, Specimen.ID.primary, Marker, r) %>%
  summarise(`Degree of Clustering Exact` = mean(`Degree of Clustering Exact`, na.rm = T)) %>%
  pivot_wider(names_from = Source, values_from = `Degree of Clustering Exact`) %>%
  mutate(`Tumor-Stroma` = Tumor - Stroma,
         Missing = case_when(is.na(Stroma) & is.na(Tumor) ~ "Both",
                             is.na(Stroma) & !is.na(Tumor) ~ "Stroma",
                             !is.na(Stroma) & is.na(Tumor) ~ "Tumor",
                             T ~ "Neither")) %>%
  right_join(clinical %>%
               select(-X, -unique_fov, -fov, -Source) %>% 
               distinct(),
             .)
#boxplots
marker_boxplots = uniK_dat_diff %>%
  filter(r == 275) %>%
  group_by(Marker) %>%
  group_map(~{
    .x %>% 
      ggplot() +
      geom_violin(aes(x = TNM.stage, y = `Tumor-Stroma`, color = TNM.stage)) +
      geom_jitter(aes(x = TNM.stage, y = `Tumor-Stroma`, color = TNM.stage, shape = Site)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = .y$Marker, y = "Tumor DOCE - Stroma DOCE") +
      scale_x_discrete(drop = F) + 
      scale_color_discrete(drop = F)
  })

pdf("Manley_SMI/results/papillary/tumor-stroma_tnmstage_r275_boxplot.pdf")
tmp = lapply(marker_boxplots, print)
dev.off()

b_cell_it = uniK_dat %>%
  filter(r == 275, Marker == "B cell") %>%
  ggplot() +
  geom_line(aes(x = Source, y = `Degree of Clustering Exact`, group = Specimen.ID.primary, color = IT.Treatment.before.collection)) +
  labs(title = "B cell clustering in Papillary RCC")
pdf("Manley_SMI/results/papillary/bcell_DOCE_r275_stroma-tumor-line_therapy.pdf")
b_cell_it
dev.off()

cd8_cell_it = uniK_dat %>%
  filter(r == 275, Marker == "CD8 T cell") %>%
  ggplot() +
  geom_line(aes(x = Source, y = `Degree of Clustering Exact`, group = Specimen.ID.primary, color = IT.Treatment.before.collection)) +
  labs(title = "CD8 T cell clustering in Papillary RCC")
pdf("Manley_SMI/results/papillary/cd8t_DOCE_r275_stroma-tumor-line_therapy.pdf")
cd8_cell_it
dev.off()

#scatter plots
marker_scatter = uniK_dat_diff %>%
  filter(r == 150) %>%
  group_by(Marker) %>%
  group_map(~{
    .x %>% ggplot() +
      geom_point(aes(x = Tumor, y = Stroma, color = TNM.stage, shape = Site)) +
      geom_smooth(aes(x = Tumor, y = Stroma, group = TNM.stage, color = TNM.stage), method = lm, se = FALSE) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = .y$Marker) +
      scale_color_discrete(drop = F)
  })

pdf("Manley_SMI/results/papillary/tumor-stroma_tnmstage_r275_scatter.pdf")
tmp = lapply(marker_scatter, print)
dev.off()

#using all tumor_stroma by stage
marker_boxplot_tissue = uniK_dat %>%
  filter(r %in% c(75, 150, 275, 500)) %>%
  group_by(Marker) %>%
  group_map(~{
    .x %>% 
      ggplot() +
      geom_boxplot(aes(x = TNM.stage, y = `Degree of Clustering Exact`, color = Source)) +
      scale_x_discrete(drop = F) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = .y$Marker) +
      facet_grid(r~., scales = "free_y")
  })

pdf("Manley_SMI/results/papillary/source_tnmstage_r75-150-300-500_boxplot.pdf")
tmp = lapply(marker_boxplot_tissue, function(p) try(print(p)))
dev.off()



#not enough samples to even try
# #testing modeling
# tests_dat = uniK_dat_diff %>%
#   filter(r %in% c(50),
#          Marker == "CD8 T cell")
# glm(`Tumor-Stroma` ~ TNM.stage, data = tests_dat) %>% summary()
# glm(`Tumor-Stroma` ~ Sarcomatoid, data = tests_dat ) %>% summary()
# 
# good_markers = uniK_dat %>%
#   filter(!is.na(`Degree of Clustering Exact`),
#          r == 150) %>%
#   group_by(Marker, Source) %>%
#   summarise(Numbs = n()) %>%
#   filter(Numbs > 2) %>%
#   group_by(Marker) %>%
#   summarise(tissues = n()) %>%
#   filter(tissues == 2) %>%
#   pull(Marker)
# 
# paired_tumor_stroma = uniK_dat %>%
#   group_by(Marker) %>%
#   filter(!is.na(`Degree of Clustering Exact`),
#          r == 150,
#          Marker %in% good_markers) %>%
#   group_map(~{
#     mod = lmerTest::lmer(`Degree of Clustering Exact` ~ Source + (1 | Specimen.ID.primary), .x)
#     summary(mod)$coefficients %>% data.frame(check.names = F) %>% rownames_to_column("Level") %>% mutate(Marker = .y$Marker)
#   }) %>%
#   do.call(bind_rows, .)
# 
# #showing the paired T-test results for degree of clustering between stroma and tumor
# paired_tumor_stroma_plot = paired_tumor_stroma %>%
#   filter(grepl("Source", Level)) %>%
#   mutate(`Higher In` = ifelse(Estimate < 0, "Stroma", "Tumor")) %>%
#   ggplot() +
#   geom_hline(yintercept = -log10(0.05), color = "red") +
#   geom_bar(aes(x = reorder(Marker, `Pr(>|t|)`), y = -log10(`Pr(>|t|)`), fill = `Higher In`), stat = "identity", color = "black") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         plot.title = element_text(hjust = 0.5)) +
#   labs(title = "Paired T-test between Tumor and Stroma (r = 150)",
#        x = "Marker", y = bquote('-' ~log[10] (p-value)))
# 
# pdf("Manley_SMI/results/figures/UnivariateRipK/source_marker_r150_paired-ttest.pdf")
# paired_tumor_stroma_plot
# dev.off()
