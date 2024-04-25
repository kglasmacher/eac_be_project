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

## More ces plots ####

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}

# BE density plot
# ggplot(data = cesa$selection$recurrent_BE) +
#   geom_point(aes(x = selection_intensity, y = 0), alpha =0.1) +
#   geom_density(aes(x = selection_intensity)) +
#   geom_label_repel(data = subset(cesa$selection$recurrent_BE, selection_intensity > 30000),
#             aes(x = selection_intensity, y = 0, label = variant_name))

be_density <- ggplot(data = cesa$selection$recurrent_BE[included_with_variant > 1]) +
  # geom_point(data = subset(cesa$selection$recurrent_BE, selection_intensity > 30000), 
  #            aes(x = selection_intensity, y = 0), alpha = 1) +
  # geom_point(data = subset(cesa$selection$recurrent_BE, selection_intensity <= 30000), 
  #            aes(x = selection_intensity, y = 0), alpha = 0.01) +
  geom_density(aes(x = selection_intensity), alpha = 0.5, fill = "lightgreen") +
  geom_label_repel(data = subset(cesa$selection$recurrent_BE[included_with_variant > 1], selection_intensity > 20000),
                   aes(x = selection_intensity, y = 0.00002, label = variant_name),
                   nudge_y = 0.0001, 
                   size = 3,
                   box.padding = unit(0.5, "lines")) +
  geom_segment(data = subset(cesa$selection$recurrent_BE[included_with_variant > 1], selection_intensity > 20000),
               aes(x = selection_intensity, y = 0, xend = selection_intensity, yend = 0.00002),
               color = "darkgreen", linewidth = 1.5) +
  labs(x = "Selection intensity", y = "Density of recurrent variants") +
  # scale_x_log10(labels=scientific) +
  theme_bw() +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))

# ggplot(data = cesa$selection$recurrent_BE[included_with_variant > 1]) +
#   geom_violin(aes(y = 1, x = selection_intensity), fill = "lightgreen", alpha = 0.5) +
#   # geom_label_repel(data = subset(cesa$selection$recurrent_BE, selection_intensity > 30000),
#   #                  aes(x = selection_intensity, y = 1, label = variant_name),
#   #                  nudge_y = 0.05) +
#   labs(x = "Selection Intensity", y = "Density") +
#   scale_fill_manual(values = c("lightgreen", "lightblue"), labels = c("<= 30000", "> 30000")) +
#   scale_x_log10() +
#   theme_bw() +
#   theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())


# ggplot(cesa$selection$recurrent_BE, aes(x = variant_name, y = selection_intensity, fill = variant_name)) +
#   geom_violin(alpha = 0.5) +
#   geom_point(data = subset(cesa$selection$recurrent_B, selection_intensity > 30000),
#              aes(x = variant_name, y = selection_intensity), color = "red") +
#   scale_y_log10() +  # Use log scale for better visualization of extreme values
#   labs(x = "Variant", y = "Selection Intensity", title = "Distribution of Selection Intensity by Variant") +
#   theme_minimal()

eac_density <- ggplot(data = cesa$selection$recurrent_EAC[included_with_variant > 1]) +
  # geom_point(data = subset(cesa$selection$recurrent_BE, selection_intensity > 30000), 
  #            aes(x = selection_intensity, y = 0), alpha = 1) +
  # geom_point(data = subset(cesa$selection$recurrent_BE, selection_intensity <= 30000), 
  #            aes(x = selection_intensity, y = 0), alpha = 0.01) +
  geom_density(aes(x = selection_intensity), alpha = 0.5, fill = "lightgreen") +
  geom_label_repel(data = subset(cesa$selection$recurrent_EAC[included_with_variant > 1], selection_intensity > 20000),
                   aes(x = selection_intensity, y = 0.0002, label = variant_name),
                   nudge_y = 0.0015) +
  geom_segment(data = subset(cesa$selection$recurrent_EAC[included_with_variant > 1], selection_intensity > 20000),
               aes(x = selection_intensity, y = 0, xend = selection_intensity, yend = 0.0002),
               color = "darkgreen", linewidth = 1.5) +
  labs(x = "Selection intensity", y = "Density of recurrent variants") +
  # scale_x_log10(labels=scientific) +
  theme_bw() +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))




# Create histogram of selection intensity
# hist_plot <- ggplot(cesa$selection$recurrent_BE, aes(x = selection_intensity)) +
#   geom_histogram(fill = "lightblue", color = "black", bins = 100) +
#   labs(title = "Distribution of Selection Intensity",
#        x = "Selection Intensity",
#        y = "Frequency") +
#   scale_x_log10()

# Create scatter plot of selection intensity with highlighted extreme values
# ggplot(cesa$selection$recurrent_BE, aes(x = variant_name, y = selection_intensity)) +
#   geom_point(color = "blue", alpha = 0.5) +
#   geom_point(data = subset(cesa$selection$recurrent_BE, selection_intensity > 25000), color = "red", size = 2) +
#   geom_label_repel(data = subset(cesa$selection$recurrent_BE, selection_intensity > 25000),
#                    aes(x = variant_name, y = selection_intensity, label = variant_name)) +
#   labs(title = "Selection Intensity for Variants in BE",
#        x = "Variants",
#        y = "Selection Intensity") +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  # Remove x-axis labels and ticks for better visualization


location_info <- cesa$variants %>%
  select(variant_id, chr, start)

BE_ces_results <- cesa$selection$recurrent_BE[included_with_variant > 1] %>%
  left_join(location_info, by = "variant_id") %>%
  arrange(chr, start) %>%
  mutate(location = row_number())


be_ces_loc <- ggplot(BE_ces_results, aes(x = location, y = selection_intensity)) +
  geom_point(color = "blue", alpha = 0.5) +
  geom_point(data = subset(BE_ces_results, selection_intensity > 20000), color = "red", size = 2) +
  geom_point(data = subset(BE_ces_results, is.na(start)), color = "black", alpha = 0.3) +
  geom_label_repel(data = subset(BE_ces_results, selection_intensity > 20000),
                   aes(x = location, y = selection_intensity, label = variant_name)) +
  labs(x = "Relative index",
       y = "Selection intensity") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())  # Remove x-axis labels and ticks for better visualization










# EAC density plot
# ggplot(data = cesa$selection$recurrent_EAC) +
#   geom_point(aes(x = selection_intensity, y = 0), alpha =0.1) +
#   geom_density(aes(x = selection_intensity)) +
#   geom_label_repel(data = subset(cesa$selection$recurrent_EAC, selection_intensity > 30000),
#                    aes(x = selection_intensity, y = 0, label = variant_name))
# 
# 
# ggplot(cesa$selection$recurrent_EAC, aes(x = variant_name, y = selection_intensity)) +
#   geom_point(color = "blue", alpha = 0.5) +
#   geom_point(data = subset(cesa$selection$recurrent_EAC, selection_intensity > 25000), color = "red", size = 2) +
#   geom_label_repel(data = subset(cesa$selection$recurrent_EAC, selection_intensity > 25000),
#                    aes(x = variant_name, y = selection_intensity, label = variant_name)) +
#   labs(title = "Selection Intensity for Variants in EAC",
#        x = "Variants",
#        y = "Selection Intensity") +
#   theme(axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank()) + # Remove x-axis labels and ticks for better visualization
#   theme_bw()


EAC_ces_results <- cesa$selection$recurrent_EAC[included_with_variant > 1] %>%
  left_join(location_info, by = "variant_id") %>%
  arrange(chr, start) %>%
  mutate(location = row_number())

eac_ces_loc <- ggplot(EAC_ces_results, aes(x = location, y = selection_intensity)) +
  geom_point(color = "blue", alpha = 0.5) +
  geom_point(data = subset(EAC_ces_results, selection_intensity > 20000), color = "red", size = 2) +
  geom_point(data = subset(EAC_ces_results, is.na(start)), color = "black", alpha = 0.3) +
  geom_label_repel(data = subset(EAC_ces_results, selection_intensity > 20000),
                   aes(x = location, y = selection_intensity, label = variant_name)) +
  labs(x = "Relative index",
       y = "Selection intensity") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) # Remove x-axis labels and ticks for better visualization


# ggarrange(be_density, be_ces_loc, eac_density, eac_ces_loc,
#           labels = c("A", "B", "C", "D"),
#           ncol = 2, nrow = 2,
#           widths = c(0.6, 1))

ggarrange(be_density, be_ces_loc, eac_density, eac_ces_loc,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2,
          widths = c(0.8, 1))

# ggarrange(be_density, be_ces_loc, eac_density, eac_ces_loc,
#           labels = c("BE", "BE", "EAC", "EAC"),
#           ncol = 2, nrow = 2,
#           widths = c(0.6, 1))
# 
# ggarrange(be_density, eac_density,
#           labels = c("BE", "EAC"),
#           ncol = 1, nrow = 2)
