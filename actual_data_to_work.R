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
                     "DO234153-SA130961", "DO234444-SA594875", "DO234284-SA599230")

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

dulak_data <- dulak_data %>%
filter(! Tumor_Sample_Barcode %in% dulak_duplicates)

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

paulson_data <- paulson_data %>%
  filter(! Tumor_Sample_Barcode %in% paulson_duplicates)

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
data.table::fwrite(janjigian_data, "input_data/janjigian.maf", sep = "\t")
data.table::fwrite(paulson_data, "input_data/paulson.maf", sep = "\t")
data.table::fwrite(naeini_data, "input_data/naeini.maf", sep = "\t")

## Get normal data from ESCC step analysis and preload ####
yokoyama_data <- read_tsv("~/../glasmacherk/Desktop/ESCC_stage_epistasis/input_data/yokoyama.maf")
yuan_data <- read_tsv("~/../glasmacherk/Desktop/ESCC_stage_epistasis/input_data/yuan.maf")
martincorena_data <- read_tsv("~/../glasmacherk/Desktop/ESCC_stage_epistasis/input_data/martincorena.maf")

covered_regions_yoko <- "~/../glasmacherk/Desktop/ESCC_stage_epistasis/targeted_regions/yokoyama_covered_regions.bed"
covered_regions_yuan <- "~/../glasmacherk/Desktop/ESCC_stage_epistasis/targeted_regions/SureSelect_All_Exon_V5_S04380110_Covered_hg19.bed"
covered_regions_mart <- "~/../glasmacherk/Desktop/ESCC_stage_epistasis/targeted_regions/martincorena_covered_regions.bed"

yokoyama_data_maf <- preload_maf(maf = yokoyama_data, refset = "ces.refset.hg19", keep_extra_columns = "Pre_or_Pri")
yokoyama_data_maf <- yokoyama_data_maf %>%
  filter(Pre_or_Pri == "Pre") %>%
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Group = "PN")

yuan_data_maf <- preload_maf(maf = yuan_data, refset = "ces.refset.hg19", keep_extra_columns = "Pre_or_Pri")
yuan_data_maf <- yuan_data_maf %>%
  filter(Pre_or_Pri == "Pre") %>%
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Group = "PN")

martincorena_data_maf <- preload_maf(maf = martincorena_data, refset = "ces.refset.hg19", keep_extra_columns = "Pre_or_Pri")
martincorena_data_maf <- martincorena_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Group = "PN")

## Preload data ####
tcga_data_maf <- preload_maf(maf = tcga_data, refset = "ces.refset.hg19", keep_extra_columns = "Group", chain_file = "~/../data/genome_data/hg38ToHg19.over.chain")
tcga_data_maf <- tcga_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Group = "EAC")

icgc_data_maf <- preload_maf(maf = icgc_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
icgc_data_maf <- icgc_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Group = "EAC")

dulak_data_maf <- preload_maf(maf = dulak_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
dulak_data_maf <- dulak_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Group = "EAC")

nones_data_maf <- preload_maf(maf = nones_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
nones_data_maf <- nones_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Group = "EAC")

janjigian_data_maf <- preload_maf(maf = janjigian_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
janjigian_data_maf <- janjigian_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Group = "EAC")

paulson_data_maf <- preload_maf(maf = paulson_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
paulson_data_maf <- paulson_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Group = "BE")

naeini_data_maf <- preload_maf(maf = naeini_data, refset = "ces.refset.hg19", keep_extra_columns = "Group")
naeini_data_maf <- naeini_data_maf %>% 
  filter(germline_variant_site == F) %>%
  filter(repetitive_region == F | cosmic_site_tier %in% 1:3) %>%
  mutate(Group = "EAC")


mafs <- list(tcga = tcga_data_maf, 
             icgc = icgc_data_maf,
             dulak = dulak_data_maf,
             nones = nones_data_maf,
             janjigian = janjigian_data_maf,
             paulson = paulson_data_maf,
             naeini = naeini_data_maf,
             yokoyama = yokoyama_data_maf,
             yuan = yuan_data_maf,
             martincorena = martincorena_data_maf) 
possible_dups <- check_sample_overlap(maf_list = mafs) # less than 3% overlap between any samples


# Set up cancereffectsizeR analysis ####
cesa <- CESAnalysis(refset = "ces.refset.hg19")

# Load in all maf files 
cesa <- load_maf(cesa, maf = tcga_data_maf, coverage = "exome", sample_data_cols = "Group", maf_name = "tcga")
cesa <- load_maf(cesa, maf = icgc_data_maf, coverage = "genome", sample_data_cols = "Group", maf_name = "icgc")
cesa <- load_maf(cesa, maf = dulak_data_maf, coverage = "exome", sample_data_cols = "Group", maf_name = "dulak")
cesa <- load_maf(cesa, maf = nones_data_maf, coverage = "genome", sample_data_cols = "Group", maf_name = "nones")
cesa <- load_maf(cesa, maf = janjigian_data_maf, coverage = "exome", sample_data_cols = "Group", maf_name = "janjigian")
cesa <- load_maf(cesa, maf = paulson_data_maf, coverage = "genome", sample_data_cols = "Group", maf_name = "paulson")
cesa <- load_maf(cesa, maf = naeini_data_maf, coverage = "genome", sample_data_cols = "Group", maf_name = "naeini")
cesa <- load_maf(cesa, maf = yuan_data_maf, coverage = "exome", sample_data_cols = "Group", maf_name = "yuan",
                 covered_regions_name = "yuan", covered_regions = covered_regions_yuan, covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = yokoyama_data_maf, coverage = "targeted", sample_data_cols = "Group", maf_name = "yokoyama",
                 covered_regions_name = "yoko", covered_regions = covered_regions_yoko, covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = martincorena_data_maf, coverage = "targeted", sample_data_cols = "Group", maf_name = "martincorena",
                 covered_regions_name = "mart", covered_regions = covered_regions_mart, covered_regions_padding = 100)

save_cesa(cesa = cesa,file = "cesa_with_normal_samples_before_generates.rds")





cesa <- gene_mutation_rates(cesa, covariates = "ESCA", samples = cesa$samples[Group=="PN"], save_all_dndscv_output = T)
cesa <- gene_mutation_rates(cesa, covariates = "ESCA", samples = cesa$samples[Group=="BE"], save_all_dndscv_output = T)
cesa <- gene_mutation_rates(cesa, covariates = "ESCA", samples = cesa$samples[Group=="EAC"], save_all_dndscv_output = T)

# Use hg19 reference set
RefCDS <- ces.refset.hg19$RefCDS

dndscv_gene_names <- cesa$gene_rates$gene

# Find number of synonymous sites
nsyn_sites <- sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

# Find number of samples in normal and tumor tissue
samples_in_PN <- length(unique(cesa$dNdScv_results$rate_grp_1$annotmuts$sampleID ))
samples_in_BE <- length(unique(cesa$dNdScv_results$rate_grp_2$annotmuts$sampleID ))
samples_in_EAC <- length(unique(cesa$dNdScv_results$rate_grp_3$annotmuts$sampleID ))

mut_rate_df <- tibble(gene = cesa$dNdScv_results$rate_grp_1$genemuts$gene_name,
                      exp_PN_mu = cesa$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
                      exp_BE_mu = cesa$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv,
                      exp_EAC_mu = cesa$dNdScv_results$rate_grp_3$genemuts$exp_syn_cv)

mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

# Mu = ((expected mu)/(number of synonymous sites))/(samples)
mut_rate_df <- mut_rate_df %>% 
  mutate(PN_mu = (exp_PN_mu / n_syn_sites) / samples_in_PN) %>%
  mutate(BE_mu = (exp_BE_mu / n_syn_sites) / samples_in_BE) %>%
  mutate(EAC_mu = (exp_EAC_mu / n_syn_sites) / samples_in_EAC) %>%
  mutate(EAC_greater_than_BE = EAC_mu > BE_mu) %>%
  mutate(BE_greater_than_PN = BE_mu > PN_mu)

# Check that mutation rate in tumor tissue (stage 0->2) is greater than mutation rate in normal tissue (stage 0->1)
mut_rate_df$EAC_greater_than_BE %>% table()
mut_rate_df$BE_greater_than_PN %>% table()

# Create mutation rates
mut_rate_df <- mut_rate_df %>%
  mutate(mut_rate_PN = PN_mu) %>% # mutation rate from "stage 0->1"
  mutate(mut_rate_BE = BE_mu - PN_mu) %>% # mutation rate from "stage 1->2"
  mutate(mut_rate_EAC = EAC_mu - BE_mu) # mutation rate from "stage 2->3"

mutation_rates <- mut_rate_df


# Get proportions of mutation rates ----
mut_rates_for_p <- mut_rate_df %>%
  select(gene, PN_mu, BE_mu, EAC_mu, mut_rate_PN, mut_rate_BE, mut_rate_EAC) %>% 
  mutate(p_1 = mut_rate_PN / EAC_mu) %>% 
  mutate(p_2 = mut_rate_BE / EAC_mu) %>%
  mutate(p_3 = mut_rate_EAC / EAC_mu)

set_EAC_rates <- mut_rate_df %>%
  select(gene,EAC_mu) %>%
  data.table::setDT()

cesa <- clear_gene_rates(cesa = cesa)

# Now all samples have highest rates
cesa <- set_gene_rates(cesa = cesa, rates = set_EAC_rates, missing_genes_take_nearest = T) 


# Infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis ----
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "Eso-AdenoCA")
cesa <- trinuc_mutation_rates(cesa, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)





# Create compound variant table ----
# Get consensus coverage across whichever samples you want to include.
# Here, we use all WES/TGS, but you could choose to exclude some if they don't cover the genes of interest well.
all_cov = c(cesa$coverage_ranges$exome, cesa$coverage_ranges$targeted)

# Exclude "exome", since typically "exome+" is what's applicable
all_cov = all_cov[! names(all_cov) == 'exome'] 
all_cov = Reduce(GenomicRanges::intersect, all_cov)

# gr argument will make only universally covered variants get returned
variants <- select_variants(cesa, gr = all_cov)

# Further filter variants table based on COSMIC oncogene/TSG classification (exclude nonrecurrent except nonsense for TSGs).
recurrent_variants <- variants[maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP") & intergenic == F]
compound <- recurrent_variants

cesa <- ces_variant(cesa = cesa, variants = recurrent_variants, run_name = "recurrent_general")
cesa <- ces_variant(cesa = cesa, variants = recurrent_variants, samples = cesa$samples[Group == "PN"], run_name = "recurrent_PN")
cesa <- ces_variant(cesa = cesa, variants = recurrent_variants, samples = cesa$samples[Group == "BE"], run_name = "recurrent_BE")
cesa <- ces_variant(cesa = cesa, variants = recurrent_variants, samples = cesa$samples[Group == "EAC"], run_name = "recurrent_EAC")







# New sequential selection method ----
source("~/../glasmacherk/Desktop/ESCC_stage_epistasis/analysis/new_sequential_lik.R")

for(comp_ind in 1:length(compound)){
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  these_props <- mut_rates_for_p[mut_rates_for_p$gene %in% this_gene, c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  if(length(this_gene) != 1){
    this_gene <- unlist(str_split(this_gene[1], "\\."))
    this_gene <- this_gene[1]
  }
  
  cesa <- ces_variant(cesa = cesa, variants = recurrent_variants, model = sequential_lik_dev, 
                      ordering_col = 'Group', ordering = c('PN', 'BE', 'EAC'), 
                      lik_args = list(sequential_mut_prop = these_props), run_name = this_gene)
  
}
