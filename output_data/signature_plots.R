library(cancereffectsizeR)
library(ggplot2)
library(cowplot)
library(data.table)
devtools::install_github("psyteachr/introdataviz")



signature_table_twostage <- cesa$mutational_signatures$biological_weights
signature_table_twostage <- signature_table_twostage %>%
  select(!c("total_snvs", "sig_extraction_snvs", "group_avg_blended"))

samples_BE <- cesa$samples %>%
  filter(Group == "BE") %>%
  select(Unique_Patient_Identifier)
samples_EAC <- cesa$samples %>%
  filter(Group == "EAC") %>%
  select(Unique_Patient_Identifier)

# Check which signatures don't have median weights of zero
signature_table_twostage[samples_BE, on = "Unique_Patient_Identifier"] %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(), median)) %>% 
  pivot_longer(cols = 1:18) %>% 
  arrange(desc(value)) %>%
  filter(value != 0)

signature_table_twostage[samples_EAC, on = "Unique_Patient_Identifier"] %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(), median)) %>% 
  pivot_longer(cols = 1:18) %>% 
  arrange(desc(value))%>%
  filter(value != 0)


# Remove signatures with all-zero weight
is_all_zero <- lapply(signature_table_twostage, function(x) all(x == 0))
cols_to_remove <- names(is_all_zero)[is_all_zero == TRUE]
signature_table_twostage <- signature_table_twostage[, -c(cols_to_remove), with = F]


# all tumors will be present in table, so no NAs
be_signature <- signature_table_twostage[samples_BE, on = "Unique_Patient_Identifier"]
# be_signature <- be_signature %>%
#   select(Unique_Patient_Identifier, SBS1, SBS2, SBS5, SBS13, SBS18)

longer_be <- be_signature %>%
  pivot_longer(!Unique_Patient_Identifier, names_to = "signatures", values_to = "signature_weight")
longer_be$signatures <- as.factor(longer_be$signatures) %>%
  fct_relevel(c("SBS1","SBS2","SBS3","SBS4","SBS5","SBS6","SBS7b","SBS7c","SBS7d",
                "SBS10c","SBS10d","SBS13","SBS16","SBS17a","SBS17b","SBS18","SBS20",
                "SBS26","SBS28","SBS34","SBS38","SBS40","SBS44","SBS86","SBS87","SBS90"))
longer_be <- longer_be %>%
  filter(signatures %in% c("SBS1","SBS5","SBS10c","SBS17a","SBS17b","SBS18","SBS40"))


be_signature_plot <- ggplot(data=longer_be, aes(x=signatures, y=signature_weight, fill=signatures)) +
  geom_violin(scale="width") +
  theme(plot.title = element_text(hjust = 0), 
        axis.text.x = element_text (hjust = 0.5, angle = 90),
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 18),
        title = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("COSMIC signature") +
  ylab("Signature weight") + ylim(0, 1) + scale_x_discrete() +
  ggtitle("BE") +
  theme(axis.line = element_line(color="black", linewidth = 0.5)) 


eac_signature <- signature_table_twostage[samples_EAC, on = "Unique_Patient_Identifier"]
# eac_signature <- eac_signature %>%
#   select(Unique_Patient_Identifier, SBS1, SBS2, SBS5, SBS13, SBS18)

longer_eac <- eac_signature %>%
  pivot_longer(!Unique_Patient_Identifier, names_to = "signatures", values_to = "signature_weight")
longer_eac$signatures <- as.factor(longer_eac$signatures) %>%
  fct_relevel(c("SBS1","SBS2","SBS3","SBS4","SBS5","SBS6","SBS7b","SBS7c","SBS7d",
                "SBS10c","SBS10d","SBS13","SBS16","SBS17a","SBS17b","SBS18","SBS20",
                "SBS26","SBS28","SBS34","SBS38","SBS40","SBS44","SBS86","SBS87","SBS90"))
  # fct_relevel(c("SBS1","SBS2","SBS5","SBS13","SBS18")) # only signatures without 0 median weight
#fct_relevel(c("SBS1","SBS2","SBS5","SBS9","SBS10c","SBS10d","SBS13","SBS15","SBS17a","SBS17b","SBS18","SBS22","SBS36","SBS86","SBS87","SBS90","SBS93"))
longer_eac <- longer_eac %>%
  filter(signatures %in% c("SBS1","SBS5","SBS10c","SBS17a","SBS17b","SBS18","SBS40"))

eac_signature_plot <- ggplot(data=longer_eac, aes(x=signatures, y=signature_weight, fill=signatures)) +
  geom_violin(scale="width") +
  theme(plot.title = element_text(hjust = 0), 
        axis.text.x = element_text (hjust = 0.5, angle = 90),
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 18),
        title = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  xlab("COSMIC signature") +
  ylab("Signature weight") + ylim(0, 1) + scale_x_discrete() +
  ggtitle("EAC") +
  theme(axis.line = element_line(color="black", linewidth = 0.5))


signature_plots <- plot_grid(be_signature_plot,eac_signature_plot, labels = c("", ""), ncol = 2)


longer_be <- longer_be %>% 
  mutate(group = "BE")
longer_eac <- longer_eac %>%
  mutate(group = "EAC")

combined_long_sigs <- rbind(longer_be, longer_eac)

ggplot(combined_long_sigs, aes(x = signatures, y = signature_weight, fill = group)) +
  # introdataviz::geom_split_violin(alpha = .4, trim = FALSE) +
  geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = T) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F,
               position = position_dodge(.2), size =0.1) +
  scale_x_discrete(name = "Signatures") +
  scale_y_continuous(name = "Signature weights") +
  scale_fill_brewer(palette = "Set1", name = "Progression") +
  theme_bw()

ggplot(combined_long_sigs, aes(x = signatures, y = signature_weight, fill = group)) +
  introdataviz::geom_split_violin(alpha = .4, trim = FALSE) +
  # geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  # stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F,
  #              position = position_dodge(.2), size =0.1) +
  scale_x_discrete(name = "Signatures") +
  scale_y_continuous(name = "Signature weights") +
  scale_fill_brewer(palette = "Set1", name = "Progression") +
  theme_bw()





EAC_ces_results %>%
  filter(variant_name == "PIK3CA_H1047L")
EAC_ces_results %>%
  filter(variant_name == "KRAS_G12D")
EAC_ces_results %>%
  filter(variant_name == "PIK3CA_E545K")
  
EAC_ces_results %>%
  filter(variant_name %in% c("PIK3CA_H1047L","KRAS_G12D","PIK3CA_E545K"))

