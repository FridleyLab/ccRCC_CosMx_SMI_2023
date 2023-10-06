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
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received Treatment")),
         Tumor.Site = ifelse(Site == "Tumor", "Primary", "Metastatic"),
         Tumor.Site = factor(Tumor.Site, levels = c("Primary", "Metastatic")),
         Sarcomatoid = factor(Sarcomatoid, levels = c("No", "Yes")),
         Sarcomatoid.Group = case_when(IT.Treatment.before.collection == "None" & Sarcomatoid == "No" ~ "Treatment Naive Non-Sarcoma",
                                       IT.Treatment.before.collection != "None" & Sarcomatoid == "No" ~ "Post IO Non-Sarcoma",
                                       T ~ "Sarcomatoid"), 
         Sarcomatoid.Group = factor(Sarcomatoid.Group, levels = c("Treatment Naive Non-Sarcoma", "Post IO Non-Sarcoma", "Sarcomatoid"))) %>%
  filter(Histology == "clear cell") %>%
  filter(fov <= 20,
         !(Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive")) %>% #remove pretreatment sarcomatoid samples
  filter(!(unique_fov %in% c("RCC5_19", "RCC5_20")))

# association with abundances ---------------------------------------------

markers = grep("\ Cells$", colnames(uniK_dat), value = T) %>%
  grep("Total", ., value = T, invert = T)

pre_post_io_stroma = lapply(markers, function(marker){
  print(marker)
  tmp = uniK_dat %>%
    filter(Source == "Stroma",
           Site == "Tumor") %>%
    select(!!marker, `Total Cells`, Pretreatment.IO, unique_fov, Tumor.Site, Source) %>%
    distinct() %>%
    rename('Positive' := marker) %>%
    mutate(Negative = `Total Cells` - Positive,
           Pretreatment.IO = as.factor(Pretreatment.IO))
  fo = as.formula(paste0("cbind(Positive, Negative) ~ Pretreatment.IO"))
  model_fit_bb = tryCatch({
    list("No Issues", vglm(fo, betabinomial(zero = 2), data = tmp))
  }, warning = function(w){
    return(list(w, vglm(fo, betabinomial(zero = 2), data = tmp)))
  })
  
  aov_bb = coefficients(summary(model_fit_bb[[2]]))
  vec = c(aov_bb[grepl("Pretreatment.IO", row.names(aov_bb)),],
          Marker = marker,
          Notes = ifelse(is(model_fit_bb[[1]], "warning"), model_fit_bb[[1]][["message"]], "No Warning/Errors"))
  return(vec)
}) %>%
  do.call(bind_rows, .) %>%
  mutate(p.adj = p.adjust(`Pr(>|z|)`, method = "fdr"), .before = "Marker")

pre_post_io_tumor = lapply(markers, function(marker){
  print(marker)
  tmp = uniK_dat %>%
    filter(Source == "Tumor",
           Site == "Tumor") %>%
    select(!!marker, `Total Cells`, Pretreatment.IO, unique_fov, Tumor.Site, Source) %>%
    distinct() %>%
    rename('Positive' := marker) %>%
    mutate(Negative = `Total Cells` - Positive,
           Pretreatment.IO = as.factor(Pretreatment.IO))
  fo = as.formula(paste0("cbind(Positive, Negative) ~ Pretreatment.IO"))
  model_fit_bb = tryCatch({
    list("No Issues", vglm(fo, betabinomial(zero = 2), data = tmp))
  }, warning = function(w){
    return(list(w, vglm(fo, betabinomial(zero = 2), data = tmp)))
  })
  
  aov_bb = coefficients(summary(model_fit_bb[[2]]))
  vec = c(aov_bb[grepl("Pretreatment.IO", row.names(aov_bb)),],
          Marker = marker,
          Notes = ifelse(is(model_fit_bb[[1]], "warning"), model_fit_bb[[1]][["message"]], "No Warning/Errors"))
  return(vec)
}) %>%
  do.call(bind_rows, .) %>%
  mutate(p.adj = p.adjust(`Pr(>|z|)`, method = "fdr"), .before = "Marker")

#showing the paired T-test results for degree of clustering between stroma and tumor
paired_tumor_stroma_plot = pre_post_io_stroma %>%
  mutate(`Higher In` = ifelse(Estimate < 0, "Treatment Naive", "Received IO\nTreatment Before"),
         `p.adj` = as.numeric(`p.adj`)) %>%
  ggplot() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_bar(aes(x = reorder(Marker, `p.adj`), y = -log10(`p.adj`), fill = `Higher In`), stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(color = "black"),
        plot.margin = margin(0, 0, 0, 2, "cm")) +
  labs(title = "Difference in Abundance of \nPre vs Post IO of Primary Tumor (Stroma Compartment)",
       x = "Marker", y = bquote('-' ~log[10] (FDR)))

paired_tumor_tumor_plot = pre_post_io_tumor %>%
  mutate(`Higher In` = ifelse(Estimate < 0, "Treatment Naive", "Received IO\nTreatment Before"),
         `p.adj` = as.numeric(`p.adj`)) %>%
  ggplot() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_bar(aes(x = reorder(Marker, `p.adj`), y = -log10(`p.adj`), fill = `Higher In`), stat = "identity", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(color = "black"),
        plot.margin = margin(0, 0, 0, 2, "cm")) +
  labs(title = "Difference in Abundance of \nPre vs Post IO of Primary Tumor (Tumor Compartment)",
       x = "Marker", y = bquote('-' ~log[10] (FDR)))

pdf("Manley_SMI/results/figures/UnivariateRipK/io_analysis/treatment-naive_vs_post-io-treatment_tumor-and-stroma_primary-abundance.pdf", 
    height = 7, width = 10)
paired_tumor_stroma_plot
paired_tumor_tumor_plot
dev.off()

wb = createWorkbook()
addWorksheet(wb, "Pre Post IO - Stroma")
addWorksheet(wb, "Pre Post IO - Tumor")
writeData(wb, sheet = 1, x = pre_post_io_stroma)
writeData(wb, sheet = 2, x = pre_post_io_tumor)
saveWorkbook(wb, "Manley_SMI/results/univar_associations/pre-post-io_primary_abundance.xlsx",
             overwrite = T)