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

change_mut_rates_plot <- ggplot() + 
  geom_point(data = selected_mut_rates_longer, mapping = aes(x = progression, y = mutation_rate, group = gene, color = gene)) +
  geom_line(data = selected_mut_rates_longer, mapping = aes(x = progression, y = mutation_rate, group = gene, color = gene)) +
  geom_text_repel(data = only_cancer_rates, mapping = aes(x = 2.05, y = mutation_rate, label = gene, color = gene), hjust = -1, direction = "y", size = 5) +
  scale_y_continuous(labels=scientific) + 
  scale_x_discrete(labels=c("Embryogenesis \nto BE", "BE \nto EAC")) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Evolutionary trajectory \n", y = "Mutation rate", color = "Gene") +
  theme_bw() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        axis.title.x = element_text(vjust = -1.5), 
        legend.position = "none")




