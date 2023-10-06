rm(list=ls())

# libraries ---------------------------------------------------------------

library(tidyverse)
library(openxlsx)
library(ggpubr)

# data --------------------------------------------------------------------

bi_ripK = readRDS("Manley_SMI/data/validation_dataframes/validation_clinical_bivariate_27um.rds")  %>%
  filter(!is.na(Source)) %>%
  mutate(Site = gsub("\\ $", "", Site),
         Sarcomatoid = gsub("\\ $", "", Sarcomatoid),
         Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received IO"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received IO"))) %>%
  mutate(unique_fov = paste0(tissue, "_FOV", fov)) %>%
  filter(Site == "Tumor",
         FOV <= 20,
         !(unique_fov %in% c("RCC5_FOV3", "RCC5_FOV13", "RCC5_FOV17", "RCC5_FOV19", "RCC5_FOV20")))

#markers
markers = bi_ripK$Anchor %>% 
  unique() %>%
  grep("Pos|\\+", ., value = T) %>%
  grep("Cyto|Nucl", ., value = T, invert = T)


# filter data to be of interest -------------------------------------------

bi_ripK2 = bi_ripK %>%
  filter(Anchor %in% markers[c(22, 23)],
         Counted %in% markers[c(22, 23)])

# cohorts -----------------------------------------------------------------

pre_IO = bi_ripK %>%
  filter(Sarcomatoid == "No", IT.Treatment.before.collection == "None") %>%
  select(unique_fov, Slide.x, FOV, Sarcomatoid, IT.Treatment.before.collection, Source, Pretreatment.IO) %>% 
  distinct()
post_IO = bi_ripK %>%
  filter(IT.Treatment.before.collection != "None") %>%
  select(unique_fov, Slide.x, FOV, Sarcomatoid, IT.Treatment.before.collection, Source, Pretreatment.IO) %>% 
  distinct()
sarc = bi_ripK %>%
  filter(Sarcomatoid == "Yes", IT.Treatment.before.collection == "None") %>%
  select(unique_fov, Slide.x, FOV, Sarcomatoid, IT.Treatment.before.collection, Source) %>% 
  distinct()


# pre v post IO -----------------------------------------------------------

bi_ripK2 %>%
  filter(unique_fov %in% c(pre_IO %>% pull(unique_fov),
                           post_IO %>% pull(unique_fov))) %>%
  ggplot() +
  geom_boxplot(aes(x = Pretreatment.IO, y = `Degree of Clustering Exact`)) +
  facet_grid(Source~Anchor) +
  theme_bw()+
  theme(plot.title = element_text(hjust= 0.5)) +
  labs(title = "Pre- and Post-IO Bivariate Ripley's K Clustering by Source")

pre_post_results = lapply(unique(pre_IO$Source), function(s){
  dat = bi_ripK2 %>%
    filter(grepl("COL4", Anchor),
           Source == s) %>%
    filter(unique_fov %in% pre_IO$unique_fov | unique_fov %in% post_IO$unique_fov)
  res = t.test(dat$`Degree of Clustering Exact` ~ dat$Pretreatment.IO)
  return(data.frame(Anchor = unique(dat$Anchor),
                    Counted = unique(dat$Counted),
                    Statistic = res$statistic,
                    P.value = res$p.value) %>%
           bind_cols(t(data.frame(res$estimate))))
})  %>%
  do.call(bind_rows, .)
row.names(pre_post_results) = unique(pre_IO$Source)


# pre v sarcomatoid -------------------------------------------------------

bi_ripK2 %>%
  filter(unique_fov %in% c(pre_IO %>% pull(unique_fov),
                           sarc %>% pull(unique_fov))) %>%
  ggplot() +
  geom_boxplot(aes(x = Sarcomatoid, y = `Degree of Clustering Exact`)) +
  facet_grid(Source~Anchor) +
  theme_bw()+
  theme(plot.title = element_text(hjust= 0.5)) +
  labs(title = "Pre- and Sarcomatoid Bivariate Ripley's K Clustering by Source")

pre_sarc_results = lapply(unique(pre_IO$Source), function(s){
  dat = bi_ripK2 %>%
    filter(grepl("COL4", Anchor),
           Source == s) %>%
    filter(unique_fov %in% pre_IO$unique_fov | unique_fov %in% sarc$unique_fov)
  res = t.test(dat$`Degree of Clustering Exact` ~ dat$Sarcomatoid)
  return(data.frame(Anchor = unique(dat$Anchor),
                    Counted = unique(dat$Counted),
                    Statistic = res$statistic,
                    P.value = res$p.value) %>%
           bind_cols(t(data.frame(res$estimate))))
})  %>%
  do.call(bind_rows, .)
row.names(pre_sarc_results) = unique(pre_IO$Source)


# save table to excel file ------------------------------------------------

wb = createWorkbook()
addWorksheet(wb, "Pre-Post IO")
addWorksheet(wb, "Pre-Sarc")

writeData(wb, 1, pre_post_results %>% rownames_to_column("Source"))
writeData(wb, 2, pre_sarc_results %>% rownames_to_column("Source"))

saveWorkbook(wb = wb, file = "Manley_SMI/results/mIF_validation/COL4A1_ITGAV_biRipK_tumor-stroma_cohorts.xlsx", overwrite = T)


# Plot samples ------------------------------------------------------------


mif = readRDS("Manley_SMI/data/UnivariateK_mif/manley_mif_validation.rds")

bi_ripK3 = bi_ripK2 %>%
  filter(grepl("COL", Counted)) %>%
  select(`Image Tag`, `Degree of Clustering Exact`, Source, unique_fov, Pretreatment.IO, Sarcomatoid) %>%
  distinct() %>%
  mutate(DOCE = round(`Degree of Clustering Exact`, 3))

pl = lapply(bi_ripK3$`Image Tag`, function(im){
  dat = mif$spatial[[im]] %>%
    mutate(xloc = (XMin + XMax)/2,
           yloc = (YMin + YMax)/2) %>%
    select(xloc, yloc, `ITGAV (Opal 520) Positive Classification`, `COL4 (Opal 540) Positive Classification`, `Classifier Label`) %>%
    rename("ITGAV+" = 3, "COL4A1+" = 4) %>%
    mutate(Background = case_when(`ITGAV+` == 1 & `COL4A1+` == 1 ~ 1,
                                  `ITGAV+` == 0 & `COL4A1+` == 0 ~ 1,
                                  T ~ 0),
           `ITGAV+` = ifelse(Background == 0 & `ITGAV+` == 1, 1, 0),
           `COL4A1+` = ifelse(Background == 0 & `COL4A1+` == 1, 1, 0)) %>%
    gather("Cell Assignment", "Positive", -xloc, -yloc, -`Classifier Label`) %>%
    filter(Positive == 1)
  dat %>%
    ggplot() +
    geom_point(dat = . %>% filter(`Cell Assignment` == "Background"), aes(x = xloc, y = yloc, shape = `Classifier Label`)) +
    geom_point(dat = . %>% filter(`Cell Assignment` != "Background"), aes(x = xloc, y = yloc, color = `Cell Assignment`, shape = `Classifier Label`)) +
    scale_shape_manual(values = c(1, 3)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 10)) + 
    labs(title = im, 
         subtitle = paste0("DOCE: ", bi_ripK3$`DOCE`[bi_ripK3$`Image Tag` == im],
                           "; ", bi_ripK3$unique_fov[bi_ripK3$`Image Tag` == im],
                           ";\n", bi_ripK3$Pretreatment.IO[bi_ripK3$`Image Tag` == im],
                           "; Sarcomaotoid: ", bi_ripK3$Sarcomatoid[bi_ripK3$`Image Tag` == im])) +
    coord_equal()
})
#arrange to make pretty
output_plots = ggarrange(plotlist = pl, ncol = 3, nrow = 2, common.legend = T)

pdf("Manley_SMI/results/mIF_validation/COL4A1_ITGAV_plots.pdf", height = 10, width = 16)
tmp = lapply(output_plots, print) 
dev.off()