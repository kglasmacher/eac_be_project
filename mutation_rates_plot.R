library(cancereffectsizeR)
library(data.table)
library(scales)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(tidyverse)



# Mutation rates plot ####
scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}

genes_of_interest <- c("TP53", 
                       "NOTCH1", 
                       "NOTCH2", 
                       "NFE2L2", 
                       "PIK3CA", 
                       "FAT1", 
                       "FBXW7",
                       "RB1") #define genes of interest
mut_rates_specific_genes <- mutation_rates %>%
  mutate(highlight = ifelse(gene %in% genes_of_interest, TRUE, FALSE)) #highlight genes of interest
selected_mut_rates <- mutation_rates %>% 
  filter(gene %in% genes_of_interest)

selected_mut_rates_longer <- selected_mut_rates %>%
  pivot_longer(cols = c("PN_mu", "BE_mu", "EAC_mu"), names_to = "progression", values_to = "mutation_rate")
selected_mut_rates_longer$progression <- factor(selected_mut_rates_longer$progression, levels = unique(selected_mut_rates_longer$progression))
only_cancer_rates <- selected_mut_rates_longer %>%
  filter(progression == "EAC_mu")

change_mut_rates_plot <- ggplot() + 
  geom_point(data = selected_mut_rates_longer, mapping = aes(x = progression, y = mutation_rate, group = gene, color = gene)) +
  geom_line(data = selected_mut_rates_longer, mapping = aes(x = progression, y = mutation_rate, group = gene, color = gene)) +
  geom_text_repel(data = only_cancer_rates, mapping = aes(x = 3.05, y = mutation_rate, label = gene, color = gene), hjust = -1, direction = "y", size = 5) +
  scale_y_continuous(labels=scientific) + 
  # scale_x_discrete(labels=c("Conception to \nadult normal tissue", "Conception \nto tumor tissue")) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Evolutionary trajectory \n", y = "Mutation rate", color = "Gene") +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        axis.title.x = element_text(vjust = -1.5), 
        legend.position = "none")




# Cancer effect size plot ####

all_ces_plot <- ggplot(data = cesa$selection$recurrent_general[selection_intensity > 5000]) +
  geom_point(aes(x = variant_name, y = selection_intensity, color = included_with_variant)) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_continuous(type = "viridis") +
  labs(x = "Variant name", y = "Selection intensity", color = "Samples with variant", title = "All samples")
  

PN_ces_plot <- ggplot(data = cesa$selection$recurrent_PN[selection_intensity > 10000]) +
  geom_point(aes(x = variant_name, y = selection_intensity, color = included_with_variant)) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_continuous(type = "viridis") +
  labs(x = "Variant name", y = "Selection intensity", color = "Samples with variant", title = "Normal samples")

BE_ces_plot <- ggplot(data = cesa$selection$recurrent_BE[selection_intensity > 2000]) +
  geom_point(aes(x = variant_name, y = selection_intensity, color = included_with_variant)) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_continuous(type = "viridis") +
  labs(x = "Variant name", y = "Selection intensity", color = "Samples with variant", title = "Barretts samples")

EAC_ces_plot <- ggplot(data = cesa$selection$recurrent_EAC[selection_intensity > 10000]) +
  geom_point(aes(x = variant_name, y = selection_intensity, color = included_with_variant)) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_continuous(type = "viridis") +
  labs(x = "Variant name", y = "Selection intensity", color = "Samples with variant", title = "EAC samples")

all_plots <- ((all_ces_plot+PN_ces_plot)/(BE_ces_plot+EAC_ces_plot))


# Compound variants ####
selection_results <- rbind(cesa$selection$TP53,
                           cesa$selection$NOTCH1,
                           cesa$selection$NOTCH2,
                           cesa$selection$ERBB4,
                           cesa$selection$NFE2L2,
                           cesa$selection$PIK3CA,
                           cesa$selection$CDKN2A,
                           cesa$selection$ARID1A,
                           cesa$selection$FAT1,
                           cesa$selection$EGFR,
                           cesa$selection$ERBB2,
                           cesa$selection$FBXW7,
                           cesa$selection$FGFR3,
                           cesa$selection$RB1,
                           cesa$selection$SMAD4,
                           cesa$selection$SOX2,
                           cesa$selection$KRAS,
                           cesa$selection$CUL3,
                           cesa$selection$ADAMTS18,
                           cesa$selection$ERBB3,
                           cesa$selection$FAT4,
                           cesa$selection$PPP1R3A,
                           cesa$selection$PREX2,
                           cesa$selection$SALL1,
                           cesa$selection$SETD2,
                           cesa$selection$SMO)

selection_results <- selection_results %>%
  filter(si_PN !=1000)
# selection_results <- selection_results %>%
#   mutate(progression = cesa$samples$Group) #take selection results and add progression column

selection_results <- selection_results %>%
  mutate(gene = gsub("\\.1.*","",variant_name)) #extract gene name from variant_name

normal <- selection_results[, .(variant = variant_name, ci_low_95_si_PN, ci_high_95_si_PN, variant_type, gene, si_PN, progression = "PN")] #filter normal samples
be <- selection_results[, .(variant = variant_name, ci_low_95_si_BE, ci_high_95_si_BE, variant_type, gene, si_BE, progression = "BE")] #filter tumor samples
eac <- selection_results[, .(variant = variant_name, ci_low_95_si_EAC, ci_high_95_si_EAC, variant_type, gene, si_EAC, progression = "EAC")] #filter tumor samples

# Rename columns and combine data frames
setnames(normal, c("si_PN", "ci_low_95_si_PN", "ci_high_95_si_PN"),
         c("selection_intensity", "ci_low_95", "ci_high_95"))
setnames(be, c("si_BE", "ci_low_95_si_BE", "ci_high_95_si_BE"),
         c("selection_intensity", "ci_low_95", "ci_high_95"))
setnames(eac, c("si_EAC", "ci_low_95_si_EAC", "ci_high_95_si_EAC"),
         c("selection_intensity", "ci_low_95", "ci_high_95"))


three_stage_results <- rbind(normal, be, eac) 


# Restrict variants recurrent within each progression
three_stage_results_normal <- three_stage_results[progression == "PN"]
three_stage_results_be <- three_stage_results[progression == "BE"]
three_stage_results_eac <- three_stage_results[progression == "EAC"]


stage_data <- data.frame(variant = three_stage_results$variant,
                         gene = three_stage_results$gene,
                         selection_intensity = three_stage_results$selection_intensity,
                         progression = three_stage_results$progression,
                         ci_low = three_stage_results$ci_low_95,
                         ci_high = three_stage_results$ci_high_95)

stage_data <- stage_data %>%
  mutate(ci_low_95 = if_else(is.na(ci_low), 0, ci_low)) #set lower bound for ci_low at 0 because it can't be negative

stage_data <- stage_data %>%
  filter(gene %in% c("TP53", 
                     "NOTCH1", 
                     "NOTCH2", 
                     "ERBB4", 
                     "NFE2L2", 
                     "PIK3CA", 
                     "CDKN2A.p16INK4a",
                     "CDKN2A.p14arf",
                     "CDKN2A",
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
                     "SMO")) #select genes of interest (genes with clear trends in this case)

stage_data$gene <- factor(stage_data$gene, levels = unique(stage_data$gene))
stage_data$gene <- reorder(stage_data$gene, -stage_data$selection_intensity)

selection_plots <- ggplot(stage_data, aes(x = progression, y = selection_intensity, color = progression)) +
  geom_point(size=2.5) +
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high, width = 0.25)) +
  labs(x="Evolutionary trajectory", y="Scaled selection coefficient", color = "Tissue type") +
  scale_color_manual(labels = c("PN", "BE", "EAC"), values = c("#F8766D", "#00BFC4", "green")) +
  theme_bw() +
  facet_wrap(~gene, ncol=4, scales = "free_y") +
  expand_limits(y = 0) +
  geom_vline(xintercept = 1.5, lwd = 0.5, color = "lightgrey") +
  geom_vline(xintercept = 2.5, lwd = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 0, lwd = 0.5, color = "lightgrey", linetype = "dotted") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size = 22),
        legend.position = "none",
        text = element_text(size = 22)) +
  scale_x_discrete(labels= expression(E %->% N, N %->% B, B %->% C))
selection_plots
