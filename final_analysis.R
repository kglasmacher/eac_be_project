library(tidyverse)
library(data.table)
library(readxl)
library(cancereffectsizeR)
library(TCGAretriever)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ces.refset.hg19)

# Read in clean data
tcga_data <- read_tsv("input_data/tcga.maf")
icgc_data <- read_tsv("input_data/icgc.maf")
dulak_data <- read_tsv("input_data/dulak.maf")
nones_data <- read_tsv("input_data/nones.maf")
janjigian_data <- read_tsv("input_data/janjigian.maf")
paulson_data <- read_tsv("input_data/paulson.maf")
naeini_data <- read_tsv("input_data/naeini.maf")



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
             naeini = naeini_data_maf) 
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

save_cesa(cesa = cesa,file = "analysis/cesa_before_generates.rds")





# calculate gene mutation rates separately for both groups
cesa <- gene_mutation_rates(cesa, covariates = "ESCA", samples = cesa$samples[Group=="BE"], save_all_dndscv_output = T)
cesa <- gene_mutation_rates(cesa, covariates = "ESCA", samples = cesa$samples[Group=="EAC"], save_all_dndscv_output = T)

# Use hg19 reference set
RefCDS <- ces.refset.hg19$RefCDS

dndscv_gene_names <- cesa$gene_rates$gene

# Find number of synonymous sites
nsyn_sites <- sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

# Find number of samples in normal and tumor tissue
samples_in_BE <- length(unique(cesa$dNdScv_results$rate_grp_1$annotmuts$sampleID ))
samples_in_EAC <- length(unique(cesa$dNdScv_results$rate_grp_2$annotmuts$sampleID ))

mut_rate_df <- tibble(gene = cesa$dNdScv_results$rate_grp_1$genemuts$gene_name,
                      exp_BE_mu = cesa$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
                      exp_EAC_mu = cesa$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv)

mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

# Mu = ((expected mu)/(number of synonymous sites))/(samples)
mut_rate_df <- mut_rate_df %>% 
  mutate(BE_mu = (exp_BE_mu / n_syn_sites) / samples_in_BE) %>%
  mutate(EAC_mu = (exp_EAC_mu / n_syn_sites) / samples_in_EAC) %>%
  mutate(EAC_greater_than_BE = EAC_mu > BE_mu) #%>%

# Check that mutation rate in tumor tissue (stage 0->2) is greater than mutation rate in normal tissue (stage 0->1)
mut_rate_df$EAC_greater_than_BE %>% table()

# Create mutation rates
mut_rate_df <- mut_rate_df %>%
  mutate(mut_rate_BE = BE_mu) %>% # mutation rate from "stage 0->1"
  mutate(mut_rate_EAC = EAC_mu - BE_mu) # mutation rate from "stage 1->2"

mutation_rates <- mut_rate_df


# Get proportions of mutation rates ----
mut_rates_for_p <- mut_rate_df %>%
  select(gene, BE_mu, EAC_mu, mut_rate_BE, mut_rate_EAC) %>% 
  mutate(p_1 = mut_rate_BE / EAC_mu) %>%
  mutate(p_2 = mut_rate_EAC / EAC_mu)

set_EAC_rates <- mut_rate_df %>%
  select(gene,EAC_mu) %>%
  data.table::setDT()

cesa <- clear_gene_rates(cesa = cesa)

# Now all samples have highest rates
cesa <- set_gene_rates(cesa = cesa, rates = set_EAC_rates, missing_genes_take_nearest = T) 


# Infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis ----
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "Eso-AdenoCA")
cesa <- trinuc_mutation_rates(cesa, signature_set = "COSMIC_v3.2", signature_exclusions = signature_exclusions)





# Create recurrent variant table ----
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

# calculate effect sizes for each variant in general
cesa <- ces_variant(cesa = cesa, variants = recurrent_variants, run_name = "recurrent_across_groups")



# calculate effect sizes for groups separately
cesa <- ces_variant(cesa = cesa, variants = recurrent_variants, samples = cesa$samples[Group == "BE"], run_name = "recurrent_BE")
cesa <- ces_variant(cesa = cesa, variants = recurrent_variants, samples = cesa$samples[Group == "EAC"], run_name = "recurrent_EAC")



# Create compound variant table ----
# Get consensus coverage across whichever samples you want to include.
# Here, we use all WES
all_cov = c(cesa$coverage_ranges$exome)

# Exclude "exome", since typically "exome+" is what's applicable
all_cov = all_cov[! names(all_cov) == 'exome'] 
all_cov = Reduce(GenomicRanges::intersect, all_cov)

# Define however many genes of interest you want
be_eac_genes = c("TP53", 
                 "NOTCH1", 
                 "NOTCH2", 
                 "ERBB4", 
                 "NFE2L2", 
                 "PIK3CA", 
                 "CDKN2A.p16INK4a",
                 "CDKN2A.p14arf",
                 "ARID1A", 
                 "FAT1", 
                 "EGFR",
                 "ERBB2",
                 "FBXW7",
                 "FGFR3",
                 "RB1",
                 "SMAD4",
                 "SOX2",
                 "KRAS",
                 "CUL3",
                 "ADAMTS18",
                 "ERBB3",
                 "FAT2",
                 "FAT3",
                 "PPP1R3A",
                 "PREX2",
                 "SALL1",
                 "SETD2",
                 "SMO")

# gr argument will make only universally covered variants get returned
variants <- select_variants(cesa, genes = be_eac_genes, gr = all_cov)

# Further filter variants table based on COSMIC oncogene/TSG classification (exclude nonrecurrent except nonsense for TSGs).
top_TP53 <- variants[gene == "TP53" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_NOTCH1 <- variants[gene == "NOTCH1" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_NOTCH2 <- variants[gene == "NOTCH2" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_ERBB4 <- variants[gene == "ERBB4" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_NFE2L2 <- variants[gene == "NFE2L2" & (maf_prevalence > 1)]
top_PIK3CA <- variants[gene == "PIK3CA" & (maf_prevalence >1)]
top_CDKN2A.p16INK4a <- variants[gene == "CDKN2A.p16INK4A" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_CDKN2A.p14arf <- variants[gene == "CDKN2A.p14arf" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_FAT1 <- variants[gene == "FAT1" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_EGFR <- variants[gene == "EGFR" & (maf_prevalence >1)]
top_ERBB2 <- variants[gene == "ERBB2" & (maf_prevalence >1)]
top_FBXW7 <- variants[gene == "FBXW7" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_FGFR3 <- variants[gene == "FGFR3" & (maf_prevalence >1)]
top_RB1 <- variants[gene == "RB1" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_SMAD4 <- variants[gene == "SMAD4" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_SOX2 <- variants[gene == "SOX2" & (maf_prevalence >1)]
top_KRAS <- variants[gene == "KRAS" & (maf_prevalence >1)]
top_CUL3 <- variants[gene == "CUL3" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_ADAMTS18 <- variants[gene == "ADAMTS18" & (maf_prevalence >1)] # unsure about function
top_ERBB3 <- variants[gene == "ERBB3" & (maf_prevalence >1)]
top_FAT4 <- variants[gene == "FAT4" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_PPP1R3A <- variants[gene == "PPP1R3A" & (maf_prevalence >1)] # unsure about function
top_PREX2 <- variants[gene == "PREX2" & (maf_prevalence >1)]
top_SALL1 <- variants[gene == "SALL1" & (maf_prevalence >1)] # unsure about function
top_SETD2 <- variants[gene == "SETD2" & (maf_prevalence > 1 | (aa_ref != "STOP" & aa_alt == "STOP") | (aa_ref == "STOP" & aa_alt != "STOP")) & intergenic == F]
top_SMO <- variants[gene == "SMO" & (maf_prevalence >1)]



for_comp <- rbind(top_TP53, 
                  top_NOTCH1, 
                  top_NOTCH2, 
                  top_ERBB4, 
                  top_NFE2L2, 
                  top_PIK3CA, 
                  top_CDKN2A.p16INK4a, 
                  top_CDKN2A.p14arf, 
                  top_FAT1, 
                  top_EGFR,
                  top_ERBB2,
                  top_FBXW7,
                  top_FGFR3,
                  top_RB1,
                  top_SMAD4,
                  top_SOX2,
                  top_KRAS,
                  top_CUL3,
                  top_ADAMTS18,
                  top_ERBB3,
                  top_FAT4,
                  top_PPP1R3A,
                  top_PREX2,
                  top_SALL1,
                  top_SETD2,
                  top_SMO)

# Filter out genes with maf prevalence less than 25 (looking at the data, 25 seemed like a sensible threshold to exclude non-significant results)
# comp_genes_high_prevalence <- for_comp %>% 
#   group_by(gene) %>% 
#   summarize(occurrence = sum(maf_prevalence)) 
#   # filter(occurrence > 25) 
# for_comp <- for_comp %>%
#   filter(gene %in% comp_genes_high_prevalence$gene)

# Define compound variants to find cancer effect sizes at the gene level and not for individual variants
compound <- define_compound_variants(cesa = cesa, variant_table = for_comp, by = "gene", merge_distance = Inf)



# New sequential selection method ----
source("analysis/new_sequential_lik.R")

for(comp_ind in 1:length(compound)){
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  these_props <- mut_rates_for_p[mut_rates_for_p$gene %in% this_gene, c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  if(length(this_gene) != 1){
    this_gene <- unlist(str_split(this_gene[1], "\\."))
    this_gene <- this_gene[1]
  }
  
  cesa <- ces_variant(cesa = cesa, variants = compound, model = sequential_lik_dev, 
                      ordering_col = 'Group', ordering = c('BE', 'EAC'), 
                      lik_args = list(sequential_mut_prop = these_props), run_name = this_gene)
  
}




# Clear gene rates and calculate gene rates for all samples (not separated by progression) for epistasis ----
cesa <- clear_gene_rates(cesa)
cesa <- gene_mutation_rates(cesa, covariates = "ESCA", samples = cesa$samples, save_all_dndscv_output = T)

dndscv_gene_names <- cesa$gene_rates$gene
nsyn_sites <- sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

samples_in_all <- length(unique(cesa$dNdScv_results$rate_grp_1$annotmuts$sampleID ))

mut_rate_df <- tibble(gene = cesa$dNdScv_results$rate_grp_1$genemuts$gene_name,
                      exp_mu = cesa$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv)

mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df <- mut_rate_df %>% 
  mutate(total_mu = (exp_mu / n_syn_sites) / samples_in_all) %>%
  select(gene, total_mu) %>%
  data.table::setDT()

cesa <- clear_gene_rates(cesa = cesa)
cesa <- set_gene_rates(cesa = cesa, rates = mut_rate_df, missing_genes_take_nearest = T) 


cesa <- ces_epistasis(cesa, variants = compound, run_name = "epistasis_compound_variants_all_samples")

save_cesa(cesa = cesa,file = "analysis/cesa_after_generates.rds")
