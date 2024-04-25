library(cancereffectsizeR)
library(data.table)

escc_cesa <- load_cesa("~/../glasmacherk/Desktop/ESCC_stage_epistasis/analysis/eso_cesa_before_generates.rds")

silent_variants = cesa$variants[aa_ref == aa_alt & essential_splice == F]

sample_groups = cesa$samples[, .(num_samples = .N, coverage = coverage[1]), by = c('Group', 'covered_regions')]

# These numbers are somewhat off since they included some noncoding regions
# target_cr_sizes = sum(width(reduce(GRangesList(cesa$coverage_ranges$targeted))))
# sample_groups[, size_Mb := target_cr_sizes[covered_regions] / 1e6]

# For WXS/WGS data, assume ~35 Mb coding genome coverage
# sample_groups[coverage != 'targeted', size_Mb := 35]

# # Some gymnastics
# counts = variant_counts(cesa, variant_ids = silent_variants$variant_id, by = c('covered_regions', 'Group'))
# counts = melt(counts[, -'variant_type'], id.vars = c('variant_id'), variable.factor = F, value.name = 'N')[variable %like% 'prevalence' & N > 0]
# counts[, covered_regions := gsub('_.*', '', variable)]
# counts = counts[covered_regions != 'total']
# counts[, Group := 'BE']
# counts[variable %like% '_EAC_', Group := 'EAC']
# final_counts = counts[, .(N = sum(N)), by = c('Group', 'covered_regions')]
# 
# sample_groups[final_counts, n_silent := N, on = c("Group", 'covered_regions')]
# sample_groups[, mean_sTMB := n_silent / num_samples / size_Mb]
# 
# final = sample_groups[, .(source = covered_regions, grp = Group, num_samples, coverage, size_Mb, mean_sTMB)][order(source)]
# final[, size_Mb := round(size_Mb, 2)]
# final[, mean_sTMB := round(mean_sTMB, 3)]
# 
















all_variants = cesa$variants

sample_groups = cesa$samples[, .(num_samples = .N, coverage = coverage[1]), by = c('Group', 'maf_source')]

sample_groups[coverage != 'targeted', size_Mb := 35]

# Some gymnastics
counts = variant_counts(cesa, variant_ids = all_variants$variant_id, by = c('maf_source', 'Group'))
counts = melt(counts[, -'variant_type'], id.vars = c('variant_id'), variable.factor = F, value.name = 'N')[variable %like% 'prevalence' & N > 0]
counts[, maf_source := gsub('_.*', '', variable)]
counts = counts[maf_source != 'total']
counts[, Group := 'BE']
counts[variable %like% '_EAC_', Group := 'EAC']
final_counts = counts[, .(N = sum(N)), by = c('Group', 'maf_source')]





all_escc_variants = escc_cesa$variants

# Some gymnastics
counts = variant_counts(escc_cesa, variant_ids = all_escc_variants$variant_id, by = c('maf_source', 'Pre_or_Pri'))
counts = melt(counts[, -'variant_type'], id.vars = c('variant_id'), variable.factor = F, value.name = 'N')[variable %like% 'prevalence' & N > 0]
counts[, maf_source := gsub('_.*', '', variable)]
counts = counts[maf_source != 'total']
counts[, Pre_or_Pri := 'Pre']
counts[variable %like% '_Pri_', Pre_or_Pri := 'Pri']
final_counts = counts[, .(N = sum(N)), by = c('Pre_or_Pri', 'maf_source')]




