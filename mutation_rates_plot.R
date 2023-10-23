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
