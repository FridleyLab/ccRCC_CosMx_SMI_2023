
# clean slate -------------------------------------------------------------

rm(list=ls())

# libraries ---------------------------------------------------------------

library(openxlsx)
library(tidyverse)

# read clinical -----------------------------------------------------------

clinical = read.csv("Manley_SMI/data/manley_files/Clinical.data_acs_02.14.23_final.csv") %>%
  mutate(Site = gsub(" $", "", Site), 
         Sarcomatoid = gsub(" $", "", Sarcomatoid),
         Pretreatment.IO = ifelse(IT.Treatment.before.collection == "None", "Treatment Naive", "Received IO"),
         Pretreatment.IO = factor(Pretreatment.IO, levels = c("Treatment Naive", "Received IO")))
clin2 = clinical %>% 
  filter(FOV <= 20,
         Site == "Tumor",
         Histology == "clear cell",
         !(FOV %in% c(3, 17) & tissue == "RCC5"),
         !(Slide == 5 & FOV %in% c(19, 20)),
         !(Slide == 4 & FOV %in% c(17)))

#patient counts
length(unique(clin2$MRN))

#sample count by tissue source
clin2 %>%
  group_by(Source) %>%
  summarise(n())

#patients by treatment
clin2 %>%
  select(MRN, IT.Treatment.before.collection) %>%
  distinct() %>%
  group_by(IT.Treatment.before.collection) %>%
  summarise(length(unique(MRN)))
clin2 %>%
  select(MRN, IT.Treatment.before.collection,OS.days) %>%
  distinct() %>%
  group_by(IT.Treatment.before.collection)%>%
  summarise(n(), median(OS.days)/(365.25/12))

#samples by treatment
clin2 %>%
  select(slide, FOV, IT.Treatment.before.collection) %>%
  distinct() %>%
  group_by(IT.Treatment.before.collection) %>%
  summarise(n())

#patients by sarcomatoid
clin2 %>%
  filter(IT.Treatment.before.collection == "None") %>%
  select(MRN, Sarcomatoid) %>%
  distinct() %>%
  group_by(Sarcomatoid) %>%
  summarise(n())
#treatments of the non-sarcomatoid samples
clin2 %>%
  select(MRN, Sarcomatoid, Pretreatment.IO) %>%
  distinct() %>%
  filter(Sarcomatoid == "No") %>%
  group_by(Pretreatment.IO) %>%
  summarise(n())

#samples by sarcomatoid
clin2 %>%
  filter(IT.Treatment.before.collection == "None") %>%
  select(slide, FOV, Sarcomatoid) %>%
  distinct() %>%
  group_by(Sarcomatoid) %>%
  summarise(n())

#overall survival
clin2 %>%
  mutate(sample_group = case_when(Sarcomatoid == "No" & Pretreatment.IO == "Treatment Naive" ~ "Treatment Naive non-Sarcomatoid",
                                  Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive" ~ "Treatment Naive Sarcomatoid", 
                                  T ~ "Received IO")) %>%
  select(sample_group, MRN, OS.days) %>%
  distinct() %>%
  group_by(sample_group) %>%
  summarise(median(OS.days)/(365.25/12))


clin2 %>%
  mutate(sample_group = case_when(Sarcomatoid == "No" & Pretreatment.IO == "Treatment Naive" ~ "Treatment Naive non-Sarcomatoid",
                                  Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive" ~ "Treatment Naive Sarcomatoid", 
                                  T ~ "Received IO")) %>%
  select(sample_group, MRN, OS.days) %>%
  distinct() %>%
  group_by(sample_group) %>%
  summarise(n())

#ns for figure
clin2 %>%
  mutate(sample_group = case_when(Sarcomatoid == "No" & Pretreatment.IO == "Treatment Naive" ~ "Treatment Naive non-Sarcomatoid",
                                  Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive" ~ "Treatment Naive Sarcomatoid", 
                                  T ~ "Received IO")) %>%
  group_by(sample_group, Source) %>%
  summarise(n())

#patients?
clin2 %>%
  mutate(sample_group = case_when(Sarcomatoid == "No" & Pretreatment.IO == "Treatment Naive" ~ "Treatment Naive non-Sarcomatoid",
                                  Sarcomatoid == "Yes" & Pretreatment.IO == "Treatment Naive" ~ "Treatment Naive Sarcomatoid", 
                                  T ~ "Received IO")) %>%
  select(sample_group, MRN) %>%
  distinct() %>%
  group_by(sample_group) %>%
  summarise(n())
