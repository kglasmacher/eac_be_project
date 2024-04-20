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
                       "ERBB3", 
                       "PIK3CA", 
                       "FAT1", 
                       "ERBB2",
                       "RB1",
                       "ERBB2",
                       "SETD2",
                       "CDKN2A.p14arf",
                       "CUL3",
                       "PPP1R3A",
                       "FBXW7",
                       "SALL1",
                       "SMAD4",
                       "SOX2",
                       "ERBB4",
                       "PREX2",
                       "ADAMTS18") #define genes of interest

mut_rates_specific_genes <- mutation_rates %>%
  mutate(highlight = ifelse(gene %in% genes_of_interest, TRUE, FALSE)) #highlight genes of interest
selected_mut_rates <- mutation_rates %>% 
  filter(gene %in% genes_of_interest)

selected_mut_rates_longer <- selected_mut_rates %>%
  pivot_longer(cols = c("BE_mu", "EAC_mu"), names_to = "progression", values_to = "mutation_rate")
selected_mut_rates_longer$progression <- factor(selected_mut_rates_longer$progression, levels = unique(selected_mut_rates_longer$progression))
only_cancer_rates <- selected_mut_rates_longer %>%
  filter(progression == "EAC_mu")

ggplot() + 
  geom_point(data = selected_mut_rates_longer, mapping = aes(x = progression, y = mutation_rate, group = gene, color = gene)) +
  geom_line(data = selected_mut_rates_longer, mapping = aes(x = progression, y = mutation_rate, group = gene, color = gene), alpha =0.3) +
  geom_text_repel(data = only_cancer_rates, mapping = aes(x = 2.02, y = mutation_rate, label = gene, color = gene), 
                  hjust = -0.3, direction = "y", size = 4, segment.size = 0.2) +
  scale_y_continuous(labels=scientific, expand = c(0,0), limits=c(0, 0.0000032)) + 
  scale_x_discrete(labels=c("Esophageal organogenesis \nto BE", "Esophageal organogenesis \nto EAC")) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Evolutionary trajectory \n", y = "Mutation rate", color = "Gene") +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        axis.title.x = element_text(vjust = -1.5), 
        legend.position = "none")



selected_mut_rates_distinct <- selected_mut_rates %>%
  pivot_longer(cols = c("mut_rate_BE", "mut_rate_EAC"), names_to = "progression", values_to = "mutation_rate_distinct")
selected_mut_rates_distinct$progression <- factor(selected_mut_rates_distinct$progression, levels = unique(selected_mut_rates_distinct$progression))
only_cancer_rates_distinct <- selected_mut_rates_distinct %>%
  filter(progression == "mut_rate_EAC")

ggplot() + 
  geom_point(data = selected_mut_rates_distinct, mapping = aes(x = progression, y = mutation_rate_distinct, group = gene, color = gene)) +
  geom_line(data = selected_mut_rates_distinct, mapping = aes(x = progression, y = mutation_rate_distinct, group = gene, color = gene), alpha =0.3) +
  geom_text_repel(data = only_cancer_rates_distinct, mapping = aes(x = 2.02, y = mutation_rate_distinct, label = gene, color = gene), 
                  hjust = -0.25, direction = "y", size = 4, segment.size = 0.2, max.overlaps = 10) +
  scale_y_continuous(labels=scientific, expand = c(0,0), limits=c(0, 0.0000018)) + 
  scale_x_discrete(labels=c("Esophageal organogenesis \nto BE", "BE \nto EAC")) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Evolutionary trajectory \n", y = "Mutation rate", color = "Gene") +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        axis.title.x = element_text(vjust = -1.5), 
        legend.position = "none")

