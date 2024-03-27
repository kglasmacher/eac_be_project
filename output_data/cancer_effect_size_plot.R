library(cancereffectsizeR)
library(data.table)
library(scales)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(patchwork)

# Cancer effect size plot ####

all_ces_plot <- ggplot(data = cesa$selection$recurrent_across_groups[selection_intensity > 7000 | included_with_variant > 20]) +
  geom_point(aes(x = variant_name, y = selection_intensity, color = included_with_variant)) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_continuous(type = "viridis") +
  labs(x = "Variant name", y = "Selection intensity", color = "Samples with variant", title = "Variants with highest selection intensity among all samples")


BE_ces_plot <- ggplot(data = cesa$selection$recurrent_BE[selection_intensity > 13000]) +
  geom_point(aes(x = variant_name, y = selection_intensity, color = included_with_variant)) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_continuous(type = "viridis") +
  labs(x = "Variant name", y = "Selection intensity", color = "Samples with variant", title = "Variants with highest selection intensity among Barrett's samples")


EAC_ces_plot <- ggplot(data = cesa$selection$recurrent_EAC[selection_intensity > 8000]) +
  geom_point(aes(x = variant_name, y = selection_intensity, color = included_with_variant)) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_continuous(type = "viridis") +
  labs(x = "Variant name", y = "Selection intensity", color = "Samples with variant", title = "Variants with highest selection intensity among EAC samples")

# all_plots <- (all_ces_plot/BE_ces_plot/EAC_ces_plot)
