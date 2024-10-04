#10 Pre IO Primary tumors vs. Post IO primary tumors and Pre IO Met tumors vs. Post IO primary tumors.
#
# For this tab and the other I removed the replicate FOV for patient 
# (FOV# greater than 20# for slide 3 and 4 and greater than 18# for slide 5; 
# this one had two samples fail so the numbers are off a little)
#
# He were are seeing how the two groups (Pre IO Primary and Pre IO mets) 
# compare to the Post IO Primary tumors. Would be “easier” to compare the 
# tumor FOVs across the cohorts and then the stroma FOVs.
rm(list=ls())
library(spatialTIME)
library(tidyverse)
library(purrr)
library(VGAM)
library(openxlsx)

manley_mif = readRDS("Manley_SMI/data/UnivariateK_mif/manley_mif.rds")
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
  right_join(clinical, .) %>%
  mutate(Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received Treatment"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received Treatment"))) %>%
  filter(Histology == "clear cell") %>%
  filter(fov <= 20,
         !(Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive")) %>% #remove pretreatment sarcomatoid samples
  filter(!(unique_fov %in% c("RCC4_17", "RCC5_19", "RCC5_20")))

# association with abundances ---------------------------------------------

markers = unique(uniK_dat$Marker)

pre_post_io_stroma = lapply(markers, function(marker){
  print(marker)
  tmp = uniK_dat %>%
    filter(Source == "Stroma",
           Site == "Tumor",
           Marker == marker,
           r == 150) %>%
    select(Pretreatment.IO, unique_fov, Site, Source, Histology, MRN, `Degree of Clustering Exact`) %>%
    distinct() %>%
    mutate(Pretreatment.IO = as.factor(Pretreatment.IO))
  fo = as.formula(paste0("`Degree of Clustering Exact` ~ Pretreatment.IO"))
  model_fit_bb = try(
    lm(fo, data = tmp), 
    silent = T)
  if(class(model_fit_bb) == "try-error"){
    return(NULL)
  }
  aov_bb = coefficients(summary(model_fit_bb))
  vec = c(aov_bb[grepl("Pretreatment", row.names(aov_bb)),],
          Marker = marker)
  return(vec)
}) %>%
  do.call(bind_rows, .) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`), .before = Marker)

pre_post_io_tumor = lapply(markers, function(marker){
  print(marker)
  tmp = uniK_dat %>%
    filter(Source == "Tumor",
           Site == "Tumor",
           Marker == marker,
           r == 150) %>%
    select(Pretreatment.IO, unique_fov, Site, Source, Histology, MRN, `Degree of Clustering Exact`) %>%
    distinct() %>%
    mutate(Pretreatment.IO = as.factor(Pretreatment.IO))
  fo = as.formula(paste0("`Degree of Clustering Exact` ~ Pretreatment.IO"))
  model_fit_bb = try(
    lm(fo, data = tmp), 
    silent = T)
  if(class(model_fit_bb) == "try-error"){
    return(NULL)
  }
  aov_bb = coefficients(summary(model_fit_bb))
  vec = c(aov_bb[grepl("Pretreatment", row.names(aov_bb)),],
          Marker = marker)
  return(vec)
}) %>%
  do.call(bind_rows, .) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`), .before = Marker)

#showing the paired T-test results for degree of clustering between stroma and tumor
paired_tumor_stroma_plot = pre_post_io_stroma %>%
  mutate(`Higher In` = ifelse(Estimate < 0, "Treatment Naive", "Received IO\nTreatment Before"),
         `Pr(>|t|)` = as.numeric(`Pr(>|t|)`)) %>%
  ggplot() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_bar(aes(x = reorder(Marker, `Pr(>|t|)`), y = -log10(`Pr(>|t|)`), fill = `Higher In`), stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(color = "black"),
        plot.margin = margin(0, 0, 0, 2, "cm")) +
  labs(title = "Difference in Degree of Clustering (Exact) of \nPre vs Post IO of Primary Tumor (Stroma Compartment, r = 150)",
       x = "Marker", y = bquote('-' ~log[10] (p-value)))

paired_tumor_tumor_plot = pre_post_io_tumor %>%
  mutate(`Higher In` = ifelse(Estimate < 0, "Treatment Naive", "Received IO\nTreatment Before"),
         `Pr(>|t|)` = as.numeric(`Pr(>|t|)`)) %>%
  ggplot() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_bar(aes(x = reorder(Marker, `Pr(>|t|)`), y = -log10(`Pr(>|t|)`), fill = `Higher In`), stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(color = "black"),
        plot.margin = margin(0, 0, 0, 2, "cm")) +
  labs(title = "Difference in Degree of Clustering (Exact) of \nPre vs Post IO of Primary Tumor (Tumor Compartment, r = 150)",
       x = "Marker", y = bquote('-' ~log[10] (p-value)))

pdf("Manley_SMI/results/figures/UnivariateRipK/io_analysis/treatment-naive_vs_post-io-treatment_tumor-and-stroma_primary-DOCE-r150.pdf", 
    height = 7, width = 10)
paired_tumor_stroma_plot
paired_tumor_tumor_plot
dev.off()

wb = createWorkbook()
addWorksheet(wb, "Pre Post IO - Stroma")
addWorksheet(wb, "Pre Post IO - Tumor")
writeData(wb, sheet = 1, x = pre_post_io_stroma)
writeData(wb, sheet = 2, x = pre_post_io_tumor)
saveWorkbook(wb, "Manley_SMI/results/univar_associations/pre-post-io_primary_DOCE-r150.xlsx",
             overwrite = T)


# plots with adjusted p values --------------------------------------------

#showing the paired T-test results for degree of clustering between stroma and tumor
paired_tumor_stroma_plot = pre_post_io_stroma %>%
  mutate(`Higher In` = ifelse(Estimate < 0, "Treatment Naive", "Received IO\nTreatment Before"),
         FDR = as.numeric(FDR)) %>%
  ggplot() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_bar(aes(x = reorder(Marker, FDR), y = -log10(FDR), fill = `Higher In`), stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(color = "black"),
        plot.margin = margin(0, 0, 0, 2, "cm")) +
  labs(title = "Difference in Degree of Clustering (Exact) of \nPre vs Post IO of Primary Tumor (Stroma Compartment, r = 150)",
       x = "Marker", y = bquote('-' ~log[10] (FDR)))

paired_tumor_tumor_plot = pre_post_io_tumor %>%
  mutate(`Higher In` = ifelse(Estimate < 0, "Treatment Naive", "Received IO\nTreatment Before"),
         FDR = as.numeric(FDR)) %>%
  ggplot() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_bar(aes(x = reorder(Marker, FDR), y = -log10(FDR), fill = `Higher In`), stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(color = "black"),
        plot.margin = margin(0, 0, 0, 2, "cm")) +
  labs(title = "Difference in Degree of Clustering (Exact) of \nPre vs Post IO of Primary Tumor (Tumor Compartment, r = 150)",
       x = "Marker", y = bquote('-' ~log[10] (FDR)))

pdf("Manley_SMI/results/figures/UnivariateRipK/io_analysis/treatment-naive_vs_post-io-treatment_tumor-and-stroma_primary-DOCE-r150_FDR.pdf", 
    height = 7, width = 10)
paired_tumor_stroma_plot
paired_tumor_tumor_plot
dev.off()

