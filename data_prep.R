library(tidyverse)
library(data.table)
library(readxl)
library(cancereffectsizeR)
library(TCGAretriever)


## TCGA data (EAC) ####

if (!file.exists("TCGA-ESCA.maf.gz")) {
  get_TCGA_project_MAF(project = "ESCA", filename = "TCGA-ESCA.maf.gz")
}
tcga_data <- fread("TCGA-ESCA.maf.gz")

tcga_clinical_data <- TCGAretriever::get_clinical_data(case_id = "esca_tcga_all") #find ESCA clinical data from TCGA

tcga_eac <- tcga_clinical_data %>%
  filter(CANCER_TYPE_DETAILED=="Esophageal Adenocarcinoma") %>% #filter for ESCC cases
  mutate(CASE_ID=substr(CASE_ID,1,12))
tcga_data <- tcga_data %>%
  filter(Unique_Patient_Identifier %in% tcga_eac$CASE_ID) #filter tcga data for just ESCC cases

tcga_data <- tcga_data %>%
  select(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(Progression = "EAC") %>%
  mutate(Chromosome = gsub("chr","",Chromosome)) #clean up tcga data frame 

fwrite(tcga_data, "tcga_eac.maf", sep = "\t")


## Dulak et al. data (EAC) ####

dulak_link <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3678719/bin/NIHMS474888-supplement-6.xlsx"
if(!file.exists("dulak_data.xlsx")){
  download.file(dulak_link, destfile = "dulak_data.xlsx") 
}
dulak_data <- readxl::read_xlsx("dulak_data.xlsx", sheet = 1, skip = 4)
dulak_data <- dulak_data %>%
  select(Tumor_Sample_Barcode = "Sample Name", Chromosome, Start_Position = Start, Reference_Allele = "Reference Allele", Tumor_Seq_Allele2 = "Tumor Allele 2") %>%
  mutate(Group = "EAC")

# dulak_cbioportal_link <- "https://cbioportal-datahub.s3.amazonaws.com/esca_broad.tar.gz"
# download.file(dulak_cbioportal_link, destfile = "dulak_cbioportal.tar.gz")


## Nones et al. data (EAC) ####

nones_link <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4596003/bin/NIHMS64791-supplement-Supplementary_data_3.xlsx"
if(!file.exists("nones_data.xlsx")){
  download.file(nones_link, destfile = "nones_data.xlsx") 
}
nones_data <- readxl::read_xlsx("nones_data.xlsx", sheet = 1, skip = 1)
nones_data <- nones_data %>%
  select(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(Group = "EAC")


## Sihag et al. data (EAC) ####

# cbioportal link: https://www.cbioportal.org/study/summary?id=egc_msk_tp53_ccr_2022
sihag_data <- read_tsv("sihag_cbioportal/data_mutations.txt")
sihag_clinical_data <- read_tsv("sihag_cbioportal/data_clinical_sample.txt", skip = 4)
sihag_clinical_data <- sihag_clinical_data %>%
  filter(CANCER_TYPE_DETAILED == "Esophageal Adenocarcinoma")

sihag_data <- sihag_data %>%
  filter(Tumor_Sample_Barcode %in% sihag_clinical_data$SAMPLE_ID) %>%
  filter(Mutation_Status == "SOMATIC") %>%
  select(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(Group = "EAC")


## Janjigian et al. data (EAC, WES) ####

# cbioportal link: https://www.cbioportal.org/study/clinicalData?id=egc_trap_msk_2020
janjigian_data <- read_tsv("janjigian_cbioportal/data_mutations.txt")
janjigian_sample_data <- read_tsv("janjigian_cbioportal/data_clinical_sample.txt", skip = 4)
janjigian_patient_data <- read_tsv("janjigian_cbioportal/data_clinical_patient.txt", skip = 4)
janjigian_sample_data <- janjigian_sample_data %>%
  filter(TISSUE_SITE == "Esophagus") %>%
  filter(SAMPLE_TYPE == "Primary") %>%
  filter(SOMATIC_STATUS == "Matched")

janjigian_data <- janjigian_data %>%
  filter(Tumor_Sample_Barcode %in% janjigian_sample_data$SAMPLE_ID) %>%
  select(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(Group = "EAC")
  
  
## Paulson et al. data (BE, WGS) ####
# https://www.nature.com/articles/s41467-022-29767-7

paulson_data <- read.csv("Source_Data_File_4 Revised20211107/snv_plus_indels.csv", sep = ";")
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-29767-7/MediaObjects/41467_2022_29767_MOESM3_ESM.xlsx", destfile = "paulson_extra_info.xlsx") 
paulson_extra_info <- readxl::read_xlsx("paulson_extra_info.xlsx", sheet = 4, skip = 1)
paulson_samples <- paulson_extra_info %>%
  mutate(Tumor_Sample_Barcode = paste0("P", Patient_ID, "-S", Sample_ID, "-", Cancer_Outcome_Status)) %>%
  select(Tumor_Sample_Barcode, DNANum = Sample_ID)
paulson_data <- left_join(paulson_data, paulson_samples, by = "DNANum")
paulson_data <- paulson_data %>%
  mutate(Group = "BE") %>%
  select(Tumor_Sample_Barcode, Chromosome = chrom, Start_position = pos, Reference_Allele = ref, Tumor_Seq_Allele2 = alt, Group)

## Stachler et al. (BE/EAC, WES) ####

# files from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4552571/bin/NIHMS696113-supplement-5.zip
# extra info from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4552571/bin/NIHMS696113-supplement-2.xlsx
stachler_files <- list.files("stachler_data_files/", full.names=TRUE)
stachler_data_list <- lapply(stachler_files, read_tsv)
names(stachler_data_list) <- substr(stachler_files, 22, 42)

stachler_data <- bind_rows(stachler_data_list, .id = "sample")
stachler_data <- stachler_data %>%
  mutate(diagnosis = case_when(str_detect(sample, "Primary") ~ "EAC",
                               str_detect(sample, "Barrett") ~ "BE")) %>%
  mutate(patient = sub("\\-.*", "", sample)) %>%
  mutate(Tumor_Sample_Barcode = paste0("Stachler-P", patient, "-", diagnosis)) %>%
  select(Tumor_Sample_Barcode, Chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele2, Group = diagnosis)


## Newell et al. (BE, WGS) !only clinical data so far! ####

newell_link <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs12920-019-0476-9/MediaObjects/12920_2019_476_MOESM2_ESM.xlsx"
if(!file.exists("newell_data.xlsx")){
  download.file(newell_link, destfile = "newell_data.xlsx") 
}
newell_clinical_data <- readxl::read_xlsx("newell_data.xlsx", sheet = 2, skip = 1)
# probably need to filter for only one sample per patient for non-progressor samples


## Naeini et al. (EAC, WGS) ####

naeini_link <- "https://europepmc.org/articles/PMC10232490/bin/41467_2023_38891_MOESM6_ESM.xlsx"
if(!file.exists("naeini_data.xlsx")){
  download.file(naeini_link, destfile = "naeini_data.xlsx") 
}
naeini_data <- readxl::read_xlsx("naeini_data.xlsx", sheet = 1, skip = 1)
naeini_data <- naeini_data %>%
  filter(Mutation_Status == "SOMATIC") %>%
  mutate(Group = "EAC") %>%
  select(Tumor_Sample_Barcode = Donor.Publication.ID, Chromosome, Start_position = Start_Position, Reference_Allele, Tumor_Seq_Allele2, Group)


## new data ####

