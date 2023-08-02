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


## Sihag et al. data (EAC) ###

# cbioportal link: https://www.cbioportal.org/study/summary?id=egc_msk_tp53_ccr_2022
sihag_data <- read_tsv("sihag_cbioportal/data_mutations.txt")
sihag_clinical_data <- read_tsv("sihag_cbioportal/data_clinical_sample.txt", skip = 4)
sihag_clinical_data <- sihag_clinical_data %>%
  filter(CANCER_TYPE_DETAILED == "Esophageal Adenocarcinoma")

sihag_data %>%
  filter(Mutation_Status == "SOMATIC") %>%
  select(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) 
