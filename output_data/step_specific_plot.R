library(cancereffectsizeR)
library(data.table)
library(scales)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(patchwork)


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

# selection_results <- selection_results %>%
#   filter(si_PN !=1000)
# selection_results <- selection_results %>%
#   mutate(progression = cesa$samples$Group) #take selection results and add progression column

selection_results <- selection_results %>%
  mutate(gene = gsub("\\.1.*","",variant_name)) #extract gene name from variant_name

be <- selection_results[, .(variant = variant_name, ci_low_95_si_BE, ci_high_95_si_BE, variant_type, gene, si_BE, progression = "BE")] #filter tumor samples
eac <- selection_results[, .(variant = variant_name, ci_low_95_si_EAC, ci_high_95_si_EAC, variant_type, gene, si_EAC, progression = "EAC")] #filter tumor samples

# Rename columns and combine data frames
setnames(be, c("si_BE", "ci_low_95_si_BE", "ci_high_95_si_BE"),
         c("selection_intensity", "ci_low_95", "ci_high_95"))
setnames(eac, c("si_EAC", "ci_low_95_si_EAC", "ci_high_95_si_EAC"),
         c("selection_intensity", "ci_low_95", "ci_high_95"))


two_stage_results <- rbind(be, eac) 


# Restrict variants recurrent within each progression
two_stage_results_be <- two_stage_results[progression == "BE"]
two_stage_results_eac <- two_stage_results[progression == "EAC"]


stage_data <- data.frame(variant = two_stage_results$variant,
                         gene = two_stage_results$gene,
                         selection_intensity = two_stage_results$selection_intensity,
                         progression = two_stage_results$progression,
                         ci_low = two_stage_results$ci_low_95,
                         ci_high = two_stage_results$ci_high_95)

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
  scale_color_manual(labels = c("BE", "EAC"), values = c("#00BFC4", "green")) +
  theme_bw() +
  facet_wrap(~gene, ncol=4, scales = "free_y") +
  expand_limits(y = 0) +
  geom_vline(xintercept = 1.5, lwd = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 0, lwd = 0.5, color = "lightgrey", linetype = "dotted") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.ticks.x = element_blank(),
        axis.title=element_text(size = 22),
        legend.position = "none",
        text = element_text(size = 22)) +
  scale_x_discrete(labels= expression(N %->% B, B %->% C))
selection_plots
