library(plotly)

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


ggplot(selection_results, aes(x=si_BE, y=si_EAC, text = paste(gene))) +
  annotate('rect', xmin=0, xmax=500, ymin=0, ymax=500, alpha=0.05, fill='green') +
  annotate('rect', xmin=0, xmax=500, ymin=500, ymax=20000, alpha=0.05, fill='red') +
  annotate('rect', xmin=500, xmax=4000, ymin=0, ymax=500, alpha=0.05, fill='blue') +
  annotate('rect', xmin=500, xmax=4000, ymin=500, ymax=20000, alpha=0.1, fill='grey') +
  geom_hline(yintercept = 500, color = "darkgrey") +
  geom_vline(xintercept = 500, color = "darkgrey") +
  geom_point(size = 2) +
  geom_label_repel(data = subset(selection_results, si_EAC > 500 & si_BE > 500),
                   aes(x=si_BE, y=si_EAC, label=gene), color = "black",
                   size = 5) +
  geom_label_repel(data = subset(selection_results, si_EAC > 500 & si_BE < 500),
                   aes(x=si_BE, y=si_EAC, label=gene), color = "darkred",
                   nudge_x = 8, nudge_y = 10,
                   size = 5) +
  geom_label_repel(data = subset(selection_results, si_EAC < 500 & si_BE > 500),
                   aes(x=si_BE, y=si_EAC, label=gene), color = "darkblue",
                   nudge_x = 5, nudge_y = 15,
                   size = 5) +
  geom_label_repel(data = subset(selection_results, si_EAC < 500 & si_BE < 500),
                   aes(x=si_BE, y=si_EAC, label=gene), color = "darkgreen",
                   nudge_x = -5, nudge_y = 10,
                   size = 5) +
  labs(x="Selection intensity in BE", 
       y="Selection intensity in EAC", 
       # color = "Frequency of variants"
       ) +
  scale_x_sqrt() +
  scale_y_sqrt() +
  scale_color_viridis_c() +
  theme_bw() +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12))
  


ggplotly(compare_si)
