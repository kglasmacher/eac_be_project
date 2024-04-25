library(cancereffectsizeR)
library(tidyverse)

cesa <- load_cesa("analysis/cesa_before_generates.rds")

samples_used <- cesa@samples

table(samples_used$maf_source)
table(samples_used$Group)

samples_table <- samples_used %>%
  mutate(Study = case_when(maf_source == "dulak" ~ "Dulak et al.",
                           maf_source == "icgc" ~ "ICGC",
                           maf_source == "janjigian" ~ "Janjigian et al.",
                           maf_source == "naeini" ~ "Naeini et al.",
                           maf_source == "nones" ~ "Nones et al.",
                           maf_source == "paulson" ~ "Paulson et al.",
                           maf_source == "tcga" ~ "TCGA")) %>%
  mutate(Sample = Unique_Patient_Identifier) %>%
  mutate(Diagnosis = case_when(Group == "BE" ~ "Barrett's Esophagus",
                               Group == "EAC" ~ "Esophageal Adenocarcinoma")) %>%
  mutate(Coverage = coverage) %>%
  select(Study, Sample, Diagnosis, Coverage)

write_tsv(samples_table, "output_data/samples.tsv")
