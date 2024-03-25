library(tidyverse)
library(data.table)
library(readxl)
library(cancereffectsizeR)
library(TCGAretriever)


## TCGA data (EAC) ####
# GRCh38
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


## ICGC (PCAWG) data (EAC, WGS) ####
# GRCh37
# https://dcc.icgc.org/releases/current/Projects/ESAD-UK
icgc_data <- read_tsv("icgc_ESAD-UK/simple_somatic_mutation.open.ESAD-UK.tsv")
icgc_donor_data <- read_tsv("icgc_ESAD-UK/donor.ESAD-UK.tsv")
icgc_sample_data <- read_tsv("icgc_ESAD-UK/sample.ESAD-UK.tsv")
icgc_specimen_data <- read_tsv("icgc_ESAD-UK/specimen.ESAD-UK.tsv")

# Looks like all of the samples in icgc_data are primary tumor samples anyways

# icgc_specimens <- icgc_specimen_data %>%
#   select(icgc_specimen_id, specimen_type)
# icgc_data_test <- left_join(icgc_data, icgc_specimens, by = "icgc_specimen_id")
# icgc_data_test <- icgc_data_test %>%
#   mutate(Tumor_Sample_Barcode = paste0(icgc_donor_id, "-", icgc_sample_id)) %>%
#   select(Tumor_Sample_Barcode, Chromosome = chromosome, Start_Position = chromosome_start, Reference_Allele = reference_genome_allele, Tumor_Seq_Allele2 = mutated_to_allele, specimen_type)

icgc_data <- icgc_data %>%
  mutate(Tumor_Sample_Barcode = paste0(icgc_donor_id, "-", icgc_sample_id)) %>%
  mutate(Group = "EAC") %>%
  select(Tumor_Sample_Barcode, Chromosome = chromosome, Start_Position = chromosome_start, Reference_Allele = reference_genome_allele, Tumor_Seq_Allele2 = mutated_to_allele, Group)

# deduplicate samples in icgc from same donor with over 20% overlap (kept the ones with most variants)
icgc_duplicates <- c("DO234158-SA130954", "DO234328-SA130951", "DO234165-SA596893", 
                     "DO234155-SA130958", "DO234162-SA594906", "DO234327-SA130945", 
                     "DO234150-SA596951", "DO234163-SA130940", "DO234165-SA596918", 
                     "DO234150-SA596942", "DO234160-SA130949", "DO234326-SA130942", 
                     "DO234153-SA130961", "DO234444-SA594875")

icgc_data <- icgc_data %>%
  filter(! Tumor_Sample_Barcode %in% icgc_duplicates)

## Dulak et al. data (EAC) ####
# GRCh37
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
# GRCh37
nones_link <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4596003/bin/NIHMS64791-supplement-Supplementary_data_3.xlsx"
if(!file.exists("nones_data.xlsx")){
  download.file(nones_link, destfile = "nones_data.xlsx") 
}
nones_data <- readxl::read_xlsx("nones_data.xlsx", sheet = 1, skip = 1)
nones_data <- nones_data %>%
  select(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(Group = "EAC")


## Sihag et al. data (EAC) ####
# GRCh37
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
# GRCh37
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
  select(Tumor_Sample_Barcode, Chromosome = chrom, Start_Position = pos, Reference_Allele = ref, Tumor_Seq_Allele2 = alt, Group)

## Stachler et al. (BE/EAC, WES) maybe get rid of this, look at notes for details####
# GRCh37
# downloaded files from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4552571/bin/NIHMS696113-supplement-5.zip to stachler_data_files/
# (extra info from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4552571/bin/NIHMS696113-supplement-2.xlsx)
stachler_files <- list.files("stachler_data_files/", full.names=TRUE)
stachler_data_list <- lapply(stachler_files, read_tsv)
names(stachler_data_list) <- substr(stachler_files, 22, 42)

stachler_data <- bind_rows(stachler_data_list, .id = "sample")
stachler_data <- stachler_data %>%
  mutate(diagnosis = case_when(str_detect(sample, "Primary") ~ "EAC",
                               str_detect(sample, "Barrett") ~ "BE")) %>%
  mutate(patient = sub("\\-.*", "", sample)) %>%
  mutate(Tumor_Sample_Barcode = paste0("Stachler-P", patient, "-", diagnosis)) %>%
  select(Tumor_Sample_Barcode, Chromosome, Start_Position = Start_position, Reference_Allele, Tumor_Seq_Allele2, Group = diagnosis)


## Newell et al. (BE, WGS) !only clinical data so far! ####

newell_link <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs12920-019-0476-9/MediaObjects/12920_2019_476_MOESM2_ESM.xlsx"
if(!file.exists("newell_data.xlsx")){
  download.file(newell_link, destfile = "newell_data.xlsx") 
}
newell_clinical_data <- readxl::read_xlsx("newell_data.xlsx", sheet = 2, skip = 1)
# probably need to filter for only one sample per patient for non-progressor samples


## Naeini et al. (EAC, WGS) ####
# GRCh37
naeini_link <- "https://europepmc.org/articles/PMC10232490/bin/41467_2023_38891_MOESM6_ESM.xlsx"
if(!file.exists("naeini_data.xlsx")){
  download.file(naeini_link, destfile = "naeini_data.xlsx") 
}
naeini_data <- readxl::read_xlsx("naeini_data.xlsx", sheet = 1, skip = 1)
naeini_data <- naeini_data %>%
  filter(Mutation_Status == "SOMATIC") %>%
  mutate(Group = "EAC") %>%
  select(Tumor_Sample_Barcode = Donor.Publication.ID, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, Group)





## Write MAFs for all ####
data.table::fwrite(tcga_data, "input_data/tcga.maf", sep = "\t")
data.table::fwrite(icgc_data, "input_data/icgc.maf", sep = "\t")
data.table::fwrite(dulak_data, "input_data/dulak.maf", sep = "\t")
data.table::fwrite(nones_data, "input_data/nones.maf", sep = "\t")
data.table::fwrite(sihag_data, "input_data/sihag.maf", sep = "\t")
data.table::fwrite(janjigian_data, "input_data/janjigian.maf", sep = "\t")
data.table::fwrite(paulson_data, "input_data/paulson.maf", sep = "\t")
data.table::fwrite(stachler_data, "input_data/stachler.maf", sep = "\t")
data.table::fwrite(naeini_data, "input_data/naeini.maf", sep = "\t")



## Preload data ####
tcga_data_maf <- preload_maf(maf = tcga_data, refset = "ces.refset.hg19", keep_extra_columns = "Group", chain_file = "~/../data/genome_data/hg38ToHg19.over.chain")
tcga_data_maf <- tcga_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3)

icgc_data_maf <- preload_maf(maf = icgc_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
icgc_data_maf <- icgc_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3)

dulak_data_maf <- preload_maf(maf = dulak_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
dulak_data_maf <- dulak_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3)

nones_data_maf <- preload_maf(maf = nones_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
nones_data_maf <- nones_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3)

sihag_data_maf <- preload_maf(maf = sihag_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
sihag_data_maf <- sihag_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3)

janjigian_data_maf <- preload_maf(maf = janjigian_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
janjigian_data_maf <- janjigian_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3)

paulson_data_maf <- preload_maf(maf = paulson_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
paulson_data_maf <- paulson_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3)

stachler_data_maf <- preload_maf(maf = stachler_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
stachler_data_maf <- stachler_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3)

naeini_data_maf <- preload_maf(maf = naeini_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
naeini_data_maf <- naeini_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3)


mafs <- list(tcga = tcga_data_maf, 
             icgc = icgc_data_maf,
             dulak = dulak_data_maf,
             nones = nones_data_maf,
             sihag = sihag_data_maf,
             janjigian = janjigian_data_maf,
             paulson = paulson_data_maf,
             stachler = stachler_data_maf,
             naeini = naeini_data_maf) 
possible_dups <- check_sample_overlap(maf_list = mafs) 




