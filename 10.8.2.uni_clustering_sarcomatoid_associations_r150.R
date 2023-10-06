#10  Sarcomatoid Primary tumors vs. Non-Sarcomatoid Primary tumors
#
# For this tab and the other I removed the replicate FOV for patient 
# (FOV# greater than 20# for slide 3 and 4 and greater than 18# for slide 5; 
# this one had two samples fail so the numbers are off a little)
#
# Here we are looking for some of the TME differences between two cohorts
# given that sarcomatoid tumors respond really well to IO treatment. Again 
# would be interesting to stratify the tumor FOVs from the stroma FOVs 
# for each patient. By “profiling” these cases we could then see if it 
# predicts response in other tumors that don’t have a formal diagnosis of 
# sarcomatoid histology. A few of the cases in the Non-sarcomatoid primary 
# tumors had excellent responses to IO and would be interesting to see if 
# their TME is more similar to the Sarcomatoid tumors.

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
         Sarcomatoid.Group = case_when(IT.Treatment.before.collection == "None" & Sarcomatoid == "No" ~ "Treatment Naive Non-Sarcoma",
                                       IT.Treatment.before.collection != "None" & Sarcomatoid == "No" ~ "Post IO Non-Sarcoma",
                                       IT.Treatment.before.collection != "None" & Sarcomatoid == "Yes" ~ "Post IO Sarcoma",
                                       T ~ "Pre IO Sarcomatoid"), 
         Sarcomatoid.Group = factor(Sarcomatoid.Group, levels = c("Treatment Naive Non-Sarcoma", "Post IO Non-Sarcoma", "Post IO Sarcoma", "Pre IO Sarcomatoid")),
         TNM.stage = factor(TNM.stage, levels = c(1,3,4))) #%>% filter(Site == "Tumor") # adding this vastly changes the signifi

uniK_dat = manley_mif$derived$univariate_Count %>%
  full_join(manley_mif$sample, .) %>%
  right_join(clinical, .) %>%
  mutate(Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received Treatment"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received Treatment")),
         Tumor.Site = ifelse(Site == "Tumor", "Primary", "Metastatic"),
         Tumor.Site = factor(Tumor.Site, levels = c("Primary", "Metastatic")),
         Sarcomatoid = factor(Sarcomatoid, levels = c("No", "Yes"))) %>%
  filter(Histology == "clear cell") %>%
  filter(fov <= 20)

# association with abundances ---------------------------------------------

markers = unique(uniK_dat$Marker)

pre_post_io_stroma = lapply(markers, function(marker){
  print(marker)
  tmp = uniK_dat %>%
    filter(Source == "Stroma",
           Site == "Tumor",
           Sarcomatoid.Group %in% c("Treatment Naive Non-Sarcoma", "Pre IO Sarcomatoid"),
           Marker == marker,
           r == 150) %>%
    select(`Total Cells`, Sarcomatoid, Pretreatment.IO, unique_fov, Tumor.Site, Source, `Degree of Clustering Exact`, Sarcomatoid.Group) %>%
    distinct() %>%
    mutate(Pretreatment.IO = as.factor(Pretreatment.IO))
  fo = as.formula(paste0("`Degree of Clustering Exact` ~ Sarcomatoid.Group"))
  model_fit_bb = try(
    lm(fo, data = tmp), 
    silent = T)
  if(class(model_fit_bb) == "try-error"){
    return(NULL)
  }
  aov_bb = coefficients(summary(model_fit_bb))
  vec = c(aov_bb[grepl("Sarcomatoid.Group", row.names(aov_bb)),],
          Marker = marker)
  return(vec)
}) %>%
  do.call(bind_rows, .) %>%
  mutate(`Pr(>|t|)` = as.numeric(`Pr(>|t|)`)) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`), .before = Marker)


pre_post_io_tumor = lapply(markers, function(marker){
  print(marker)
  tmp = uniK_dat %>%
    filter(Source == "Tumor",
           Site == "Tumor",
           Sarcomatoid.Group %in% c("Treatment Naive Non-Sarcoma", "Pre IO Sarcomatoid"),
           Marker == marker,
           r == 150) %>%
    select(`Total Cells`, Sarcomatoid, Pretreatment.IO, unique_fov, Tumor.Site, Source, `Degree of Clustering Exact`, Sarcomatoid.Group) %>%
    distinct() %>%
    mutate(Pretreatment.IO = as.factor(Pretreatment.IO))
  fo = as.formula(paste0("`Degree of Clustering Exact` ~ Sarcomatoid.Group"))
  model_fit_bb = try(
    lm(fo, data = tmp), 
    silent = T)
  if(class(model_fit_bb) == "try-error"){
    return(NULL)
  }
  aov_bb = coefficients(summary(model_fit_bb))
  vec = c(aov_bb[grepl("Sarcomatoid.Group", row.names(aov_bb)),],
          Marker = marker)
  return(vec)
}) %>%
  do.call(bind_rows, .) %>%
  mutate(`Pr(>|t|)` = as.numeric(`Pr(>|t|)`)) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`), .before = Marker)

#showing the paired T-test results for degree of clustering between stroma and tumor
paired_tumor_stroma_plot = pre_post_io_stroma %>%
  mutate(`Higher In: Sarcomatoid` = ifelse(Estimate < 0, "No", "Yes"),
         `Pr(>|t|)` = as.numeric(`Pr(>|t|)`)) %>%
  ggplot() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_bar(aes(x = reorder(Marker, `Pr(>|t|)`), y = -log10(`Pr(>|t|)`), fill = `Higher In: Sarcomatoid`), stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(color = "black"),
        plot.margin = margin(0, 0, 0, 2, "cm")) +
  labs(title = "Difference in Degree of Clustering (Exact) of Primary Tumors\n(Pre IO, Stroma Compartment, r = 150)",
       x = "Marker", y = bquote('-' ~log[10] (p-value)))

paired_tumor_tumor_plot = pre_post_io_tumor %>%
  mutate(`Higher In: Sarcomatoid` = ifelse(Estimate < 0, "No", "Yes"),
         `Pr(>|t|)` = as.numeric(`Pr(>|t|)`)) %>%
  ggplot() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_bar(aes(x = reorder(Marker, `Pr(>|t|)`), y = -log10(`Pr(>|t|)`), fill = `Higher In: Sarcomatoid`), stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(color = "black"),
        plot.margin = margin(0, 0, 0, 2, "cm")) +
  labs(title = "Difference in Degree of Clustering (Exact) of \nPrimary Tumors (Pre IO, Tumor Compartment, r = 150)",
       x = "Marker", y = bquote('-' ~log[10] (p-value)))

pdf("Manley_SMI/results/figures/UnivariateRipK/io_analysis/treatment-naive_tumor-and-stroma_sarcomatoid-DOCE-r150.pdf", 
    height = 7, width = 10)
paired_tumor_stroma_plot
paired_tumor_tumor_plot
dev.off()

wb = createWorkbook()
addWorksheet(wb, "Pre Post IO - Stroma")
addWorksheet(wb, "Pre Post IO - Tumor")
writeData(wb, sheet = 1, x = pre_post_io_stroma)
writeData(wb, sheet = 2, x = pre_post_io_tumor)
saveWorkbook(wb, "Manley_SMI/results/univar_associations/pre-post-io_sarcomatoid_DOCE-r150.xlsx",
             overwrite = T)

#showing the paired T-test results for degree of clustering between stroma and tumor
paired_tumor_stroma_plot = pre_post_io_stroma %>%
  mutate(`Higher In: Sarcomatoid` = ifelse(Estimate < 0, "No", "Yes"),
         FDR = as.numeric(FDR)) %>%
  ggplot() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_bar(aes(x = reorder(Marker, FDR), y = -log10(FDR), fill = `Higher In: Sarcomatoid`), stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(color = "black"),
        plot.margin = margin(0, 0, 0, 2, "cm")) +
  labs(title = "Difference in Degree of Clustering (Exact) of Primary Tumors\n(Pre IO, Stroma Compartment, r = 150)",
       x = "Marker", y = bquote('-' ~log[10] (FDR)))

paired_tumor_tumor_plot = pre_post_io_tumor %>%
  mutate(`Higher In: Sarcomatoid` = ifelse(Estimate < 0, "No", "Yes"),
         FDR = as.numeric(FDR)) %>%
  ggplot() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_bar(aes(x = reorder(Marker, FDR), y = -log10(FDR), fill = `Higher In: Sarcomatoid`), stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(color = "black"),
        plot.margin = margin(0, 0, 0, 2, "cm")) +
  labs(title = "Difference in Degree of Clustering (Exact) of \nPrimary Tumors (Pre IO, Tumor Compartment, r = 150)",
       x = "Marker", y = bquote('-' ~log[10] (FDR)))

pdf("Manley_SMI/results/figures/UnivariateRipK/io_analysis/treatment-naive_tumor-and-stroma_sarcomatoid-DOCE-r150_FDR.pdf", 
    height = 7, width = 10)
paired_tumor_stroma_plot
paired_tumor_tumor_plot
dev.off()




