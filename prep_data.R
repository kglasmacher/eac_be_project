library(tidyverse)
library(data.table)
library(readxl)
library(cancereffectsizeR)
library(TCGAretriever)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ces.refset.hg19)

## TCGA data (EAC) ####
# reference genome: GRCh38
if (!file.exists("TCGA-ESCA.maf.gz")) {
  get_TCGA_project_MAF(project = "ESCA", filename = "TCGA-ESCA.maf.gz")
}
tcga_data <- fread("TCGA-ESCA.maf.gz")

tcga_clinical_data <- TCGAretriever::get_clinical_data(csid = "stes_tcga_pub", case_list_id = "stes_tcga_pub_sequenced") #find ESCA clinical data from TCGA

tcga_eac <- tcga_clinical_data %>%
  filter(CANCER_TYPE_DETAILED=="Esophageal Adenocarcinoma") #filter for ESCC cases
tcga_data <- tcga_data %>%
  filter(Unique_Patient_Identifier %in% tcga_eac$patientId) #filter tcga data for just ESCC cases

tcga_data <- tcga_data %>%
  select(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(Progression = "EAC") %>%
  mutate(Chromosome = gsub("chr","",Chromosome)) #clean up tcga data frame 
#tcga only has a single variant shared with others or if more, only 0.01 greater overlap
#tcga only has one pair with 1.3% greater overlap, so no further exclusion of samples


## ICGC (PCAWG) data (EAC, WGS) ####
# reference genome: GRCh37
# https://dcc.icgc.org/releases/current/Projects/ESAD-UK
icgc_data <- read_tsv("icgc_ESAD-UK/simple_somatic_mutation.open.ESAD-UK.tsv")
icgc_donor_data <- read_tsv("icgc_ESAD-UK/donor.ESAD-UK.tsv")
icgc_sample_data <- read_tsv("icgc_ESAD-UK/sample.ESAD-UK.tsv")
icgc_specimen_data <- read_tsv("icgc_ESAD-UK/specimen.ESAD-UK.tsv")

# all of the samples in icgc_data are primary tumor samples, so the following lines of code are not necessary
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

#icgc no more than 3 variants shared with any other samples from other data sources

# deduplicate samples in icgc from same donor with over 20% overlap (kept the ones with most variants)
icgc_duplicates <- c("DO234158-SA130954", "DO234328-SA130951", "DO234165-SA596893", 
                     "DO234155-SA130958", "DO234162-SA594906", "DO234327-SA130945", 
                     "DO234150-SA596951", "DO234163-SA130940", "DO234165-SA596918", 
                     "DO234150-SA596942", "DO234160-SA130949", "DO234326-SA130942", 
                     "DO234153-SA130961", "DO234444-SA594875", "DO234284-SA599230")

icgc_data <- icgc_data %>%
  filter(! Tumor_Sample_Barcode %in% icgc_duplicates)

## Dulak et al. data (EAC) ####
# GRCh37
dulak_link <- "https://pmc.ncbi.nlm.nih.gov/articles/instance/3678719/bin/NIHMS474888-supplement-6.xlsx"
if(!file.exists("dulak_data.xlsx")){
  download.file(dulak_link, destfile = "dulak_data.xlsx") 
}
dulak_data <- readxl::read_xlsx("dulak_data.xlsx", sheet = 1, skip = 4)
dulak_data <- dulak_data %>%
  select(Tumor_Sample_Barcode = "Sample Name", Chromosome, Start_Position = Start, Reference_Allele = "Reference Allele", Tumor_Seq_Allele2 = "Tumor Allele 2") %>%
  mutate(Group = "EAC")

# dulak_cbioportal_link <- "https://cbioportal-datahub.s3.amazonaws.com/esca_broad.tar.gz"
# download.file(dulak_cbioportal_link, destfile = "dulak_cbioportal.tar.gz")

#dulak no more than 6 variants shared with any other samples in other sources
#within dulak data: a cluster of overlaps -> take only the sample with the most variants 
dulak_duplicates <- c("ESO-141-Tumor", "ESO-160-Tumor", "ESO-0125-Tumor", "ESO-0009-Tumor", "ESO-118-Tumor")

#overlap within the dataset and same patients, so filtering out samples is necessary
dulak_data <- dulak_data %>%
  filter(! Tumor_Sample_Barcode %in% dulak_duplicates)

## Nones et al. data (EAC) ####
# GRCh37
nones_link <- "https://pmc.ncbi.nlm.nih.gov/articles/instance/4596003/bin/NIHMS64791-supplement-Supplementary_data_3.xlsx"
if(!file.exists("nones_data.xlsx")){
  download.file(nones_link, destfile = "nones_data.xlsx") 
}
nones_data <- readxl::read_xlsx("nones_data.xlsx", sheet = 1, skip = 1)
nones_data <- nones_data %>%
  select(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  mutate(Group = "EAC")
#nones no more than 1 variant shared with any other samples in other data sources and no overlap within


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
# no more than 1 variant shared with any other samples from other data sources and also no overlaps within


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

# duplicates within paulson dataset
paulson_duplicates <- c("P396-S24461-CO", "P88-S23701-NCO", "P163-S24735-CO", "P521-S23579-CO", "P184-S24325-NCO", "P184-S24328-NCO", 
  "P381-S23842-NCO",  "P954-S24220-CO",   "P652-S23878-NCO",  "P184-S24331-NCO",  "P635-S24783-NCO",  "P483-S23824-NCO", 
  "P126-S24349-NCO", "P852-S23908-CO",   "P635-S24040-NCO",  "P286-S24364-CO",   "P568-S24885-CO",   "P611-S24843-NCO", 
  "P954-S24205-CO","P43-S23593-CO",    "P42-S24283-CO",    "P997-S24595-CO",   "P572-S23785-NCO",  "P891-S24873-CO",  
  "P954-S24214-CO",   "P909-S24831-NCO",  "P381-S23836-NCO",  "P740-S24250-CO",   "P295-S24184-CO",   "P126-S24355-NCO", 
  "P512-S23971-CO",   "P126-S24631-NCO",  "P779-S24085-CO",   "P222-S24373-NCO",  "P55-S23614-NCO",   "P433-S23896-NCO", 
  "P740-S24253-CO",   "P286-S24361-CO",   "P126-S24343-NCO",  "P88-S23692-NCO",   "P626-S24028-CO",   "P126-S24352-NCO", 
  "P686-S24663-CO",   "P911-S24628-NCO",  "P184-S24340-NCO",  "P772-S24571-CO",   "P660-S24058-CO",   "P387-S23656-CO",  
  "P381-S23851-NCO",  "P572-S23779-NCO",  "P169-S23764-CO",   "P977-S24864-CO",   "P422-S25123-CO",   "P74-S23677-CO",   
  "P977-S23860-CO",   "P997-S24592-CO",   "P512-S23977-CO",   "P512-S23968-CO",   "P551-S24436-CO",   "P909-S24514-NCO", 
  "P222-S24376-NCO",  "P422-S25117-CO",   "P660-S24052-CO",   "P609-S24915-NCO",  "P483-S23830-NCO",  "P852-S23905-CO",  
  "P360-S24810-NCO",  "P951-S23806-CO",   "P729-S23752-NCO",  "P862-S24589-NCO",  "P779-S24091-CO",   "P295-S24175-CO",  
  "P997-S24601-CO",   "P900-S23926-NCO",  "P422-S25129-CO",   "P160-S23719-CO",   "P772-S24577-CO",   "P909-S24523-NCO", 
  "P303-S23605-NCO",  "P609-S24007-NCO",  "P478-S23953-NCO",  "P635-S24046-NCO",  "P381-S23845-NCO",  "P635-S24049-NCO", 
  "P450-S23920-CO",   "P478-S23944-NCO",  "P1047-S24277-NCO", "P672-S25156-CO",   "P184-S24334-NCO",  "P911-S24619-NCO", 
  "P222-S24379-NCO",  "P660-S24897-CO",   "P639-S24235-NCO",  "P360-S24409-NCO",  "P619-S24463-NCO",  "P74-S23668-CO",   
  "P833-S24774-CO",   "P295-S24172-CO",   "P686-S24669-CO",   "P169-S23773-CO",   "P609-S23992-NCO",  "P1005-S24106-NCO",
  "P631-S23728-NCO",  "P619-S24466-NCO",  "P303-S23602-NCO",  "P450-S24753-CO",   "P623-S23713-CO",   "P956-S24765-NCO", 
  "P160-S23725-CO",   "P995-S24945-CO",   "P1005-S25105-NCO", "P279-S24924-CO",   "P672-S24553-CO",   "P977-S23857-CO",  
  "P483-S23833-NCO",  "P433-S24678-NCO",  "P609-S24004-NCO",  "P391-S23869-CO",   "P403-S24388-CO",   "P163-S23743-CO",  
  "P740-S24256-CO",   "P915-S25171-NCO",  "P74-S24717-CO",    "P1047-S24274-NCO", "P672-S24550-CO",   "P163-S23749-CO",  
  "P900-S23932-NCO",  "P995-S24130-CO",   "P541-S24232-NCO",  "P387-S23665-CO",   "P360-S24412-NCO",  "P631-S23731-NCO", 
  "P862-S24580-NCO",  "P652-S23884-NCO",  "P833-S24777-CO",   "P728-S25162-CO",   "P286-S24367-CO",   "P909-S24520-NCO", 
  "P909-S24529-NCO",  "P59-S23638-NCO",   "P942-S24906-NCO",  "P59-S23653-NCO",   "P686-S25096-CO",   "P956-S23980-NCO", 
  "P900-S23929-NCO",  "P639-S24244-NCO",  "P575-S25150-NCO",  "P779-S24082-CO",   "P968-S24193-NCO",  "P387-S23662-CO",  
  "P609-S24768-NCO",  "P55-S23611-NCO",   "P170-S23788-NCO",  "P833-S24016-CO",   "P391-S23866-CO",   "P322-S24807-CO",  
  "P55-S23623-NCO",   "P1005-S24097-NCO", "P575-S25147-NCO",  "P391-S23521-CO",   "P55-S23617-NCO",   "P862-S24586-NCO", 
  "P631-S24729-NCO",  "P575-S25141-NCO",  "P1047-S24280-NCO", "P626-S24031-CO",   "P170-S23797-NCO",  "P160-S23722-CO",  
  "P17-S23561-NCO",   "P521-S24654-CO",   "P652-S25111-NCO",  "P17-S24933-NCO",   "P551-S24433-CO",   "P951-S23800-CO",  
  "P170-S23791-NCO",  "P130-S24304-NCO",  "P915-S24490-NCO",  "P806-S23341-NCO",  "P951-S23803-CO",   "P450-S23923-CO",  
  "P17-S23564-NCO",   "P995-S24136-CO",   "P266-S24139-CO",   "P626-S24780-CO",   "P728-S24265-CO",   "P130-S24301-NCO", 
  "P865-S24070-NCO",  "P381-S23839-NCO",  "P911-S24625-NCO",  "P806-S24160-NCO",  "P59-S23644-NCO",   "P385-S24816-NCO", 
  "P772-S24574-CO",   "P968-S24196-NCO",  "P728-S25165-CO",   "P865-S24789-NCO",  "P856-S24756-CO",   "P385-S24819-NCO", 
  "P266-S24145-CO",   "P88-S23698-NCO",   "P541-S24229-NCO",  "P396-S24448-CO",   "P865-S24079-NCO",  "P59-S23650-NCO",  
  "P942-S24124-NCO",  "P623-S24672-CO",   "P891-S24616-CO",   "P942-S24795-NCO",  "P611-S24837-NCO",  "P303-S24660-NCO", 
  "P551-S24439-CO",   "P856-S24891-CO",   "P541-S24223-NCO",  "P478-S23950-NCO",  "P729-S23758-NCO",  "P619-S25183-NCO", 
  "P322-S24403-CO",   "P568-S24690-CO",   "P42-S24286-CO",    "P729-S23761-NCO",  "P403-S24385-CO",   "P42-S24292-CO",   
  "P865-S24786-NCO",  "P611-S24559-NCO",  "P478-S23938-NCO",  "P403-S24383-CO",   "P43-S24696-CO",    "P521-S24657-CO",  
  "P385-S24424-NCO",  "P59-S23641-NCO",   "P623-S24879-CO",   "P130-S24295-NCO",  "P515-S23689-NCO",  "P702-S24855-CO",  
  "P322-S24394-CO",   "P515-S23686-NCO",  "P515-S24723-NCO",  "P891-S24613-CO",   "P806-S24154-NCO",  "P956-S23986-NCO", 
  "P55-S23620-NCO",   "P635-S24043-NCO",  "P266-S24148-CO",   "P915-S25168-NCO",  "P865-S24076-NCO",  "P702-S24858-CO")
paulson_data <- paulson_data %>%
  filter(! Tumor_Sample_Barcode %in% paulson_duplicates)
# no more than 1 variant shared with any other samples from other data sets or less than 4% overlap


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
# no data overlap within dataset and no more than 1 variant shared with any other samples or less than 1.2% overlap




## Write MAFs for all ####
data.table::fwrite(tcga_data, "input_data/tcga.maf", sep = "\t")
data.table::fwrite(icgc_data, "input_data/icgc.maf", sep = "\t")
data.table::fwrite(dulak_data, "input_data/dulak.maf", sep = "\t")
data.table::fwrite(nones_data, "input_data/nones.maf", sep = "\t")
data.table::fwrite(janjigian_data, "input_data/janjigian.maf", sep = "\t")
data.table::fwrite(paulson_data, "input_data/paulson.maf", sep = "\t")
data.table::fwrite(naeini_data, "input_data/naeini.maf", sep = "\t")
