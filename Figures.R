#Figures:

lancet_colors <- c("grey70","#76b5c5","#f28e79","#a1c181","#f7bb5f","#4ebdb6","#ff9da7","#f4d4a4","#b699d7","#88d8b0")





##------------------------------------------------------------------------------------------------------##
##------------------------------------------------------------------------------------------------------##
##------------------------------------------------------------------------------------------------------##

#MRSA FIRST DATA MANAGEMENT
results_matrix_dea1 <- as.data.frame(cbind(id, results_epi_MRSAdead))
colnames(results_matrix_dea1) <- c("id", "S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_ep_dea1 <- melt(results_matrix_dea1, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))
results_ep_dea1 <- results_ep_dea1[!grepl("S7|S8|S9|S10|S11|S12", results_ep_dea1$variable), ]
results_ep_dea1$variable <- gsub("S13", "S7", results_ep_dea1$variable)
results_ep_dea1$variable <- gsub("S14", "S8", results_ep_dea1$variable)
results_ep_dea1$variable <- gsub("S15", "S9", results_ep_dea1$variable)
deaths_p365_mrsa <- (subset(results_ep_dea1, id == 366)$value/(1000*365))*100000
deaths_p365_mrsax <- data.frame(
  Strategy = paste("S", 0:9, sep=""),
  start_val = 19, # Replace these with your actual start values
  end_val = deaths_p365_mrsa   # Replace these with your actual end values
)
deaths_p365_mrsax$StrategyNum <- 1:nrow(deaths_p365_mrsax)
dif_CI_mrsaD2<- ((dif_CI_mrsaD)/(1000*365))*100000
deaths_p365_mrsax$ci_lower = deaths_p365_mrsax$end_val - (dif_CI_mrsaD2)
deaths_p365_mrsax$ci_upper = deaths_p365_mrsax$end_val + (dif_CI_mrsaD2)
#CRE FIRST DATA MANAGEMENT
results_epi_CREdead<- as.data.frame(results_epi_CREdead)
results_matrix_dea2 <- as.data.frame(cbind(id, results_epi_CREdead))
colnames(results_matrix_dea2) <- c("id", "S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_ep_dea2 <- melt(results_matrix_dea2, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))
results_ep_dea2 <- results_ep_dea2[!grepl("S7|S8|S9|S10|S11|S12", results_ep_dea2$variable), ]
results_ep_dea2$variable <- gsub("S13", "S7", results_ep_dea2$variable)
results_ep_dea2$variable <- gsub("S14", "S8", results_ep_dea2$variable)
results_ep_dea2$variable <- gsub("S15", "S9", results_ep_dea2$variable)
deaths_p365_cre <- ((subset(results_ep_dea2, id == 366)$value)/(1000*365))*100000
deaths_p365_crex <- data.frame(
  Strategy = paste("S", 0:9, sep=""),
  start_val = 20, # Replace these with your actual start values
  end_val = deaths_p365_cre   # Replace these with your actual end values
)
deaths_p365_crex$StrategyNum <- 1:nrow(deaths_p365_crex)
dif_CI_creD2<- ((dif_CI_creD)/(1000*365))*100000
deaths_p365_crex$ci_lower = deaths_p365_crex$end_val - dif_CI_creD2
deaths_p365_crex$ci_upper = deaths_p365_crex$end_val + dif_CI_creD2
#MRSA GRAPH:
abc2<-ggplot(deaths_p365_mrsax, aes(x = StrategyNum, ymin=start_val,  ymax = end_val, fill = Strategy)) +
  geom_bar(stat="identity", aes(y=end_val), position="dodge") +
  geom_rect(aes(xmin = StrategyNum - 0.44, xmax = StrategyNum + 0.44), color = "black") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.4) +
  #geom_bar(stat = 'identity', color= 'black') + # Use 'identity' to use the actual values in the data
  theme_minimal() +
  theme_lancet()+
  labs(title = "", x = "Strategy", y = "End Value") +
  scale_fill_manual(values = lancet_colors, labels = c("S0: Do-nothing", expression("S1: T+D"["agar1"]), expression("S2: T+D"["agar2"]),
                                                       expression("S3: T+D"["pcr"]), expression("S4: T+I"["agar1"]), expression("S5: T+I"["agar2"]),
                                                       expression("S6: T+I"["pcr"]), expression("S7: Pre-emptive I"["all"]),
                                                       expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"])))+
  labs(subtitle = "(B) MRSA") +
  theme(plot.subtitle = element_text(hjust = 0, size = 16, face="bold"), legend.position = "right", legend.box.just = "left", legend.justification = "left", legend.title = element_text(size = 15), legend.text = element_text(size = 10, hjust = 0)) +
  scale_y_continuous(name = "Annual mortality rate per\n100,000 hospital bed-days", labels = scales::label_number(accuracy = 1), breaks = seq(20, 70, by = 5))+
  scale_x_discrete(name = "Strategy") +
  geom_hline(yintercept = 53.1, color = "grey", linetype = "dashed", size = 1.1)+
  coord_cartesian(ylim = c(21.3, NA))
  
#CRE GRAPH:
abc1<-ggplot(deaths_p365_crex, aes(x = StrategyNum, ymin=start_val,  ymax = end_val, fill = Strategy))+
  geom_bar(stat="identity", aes(y=end_val), position="dodge") +
  geom_rect(aes(xmin = StrategyNum - 0.44, xmax = StrategyNum + 0.44), color = "black") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.4) +
  theme_minimal() +
  theme_lancet()+
  labs(title = "", x = "Strategy", y = "End Value") +
  scale_fill_manual(values = lancet_colors, labels = c("S0: Do-nothing", expression("S1: T+D"["agar1"]), expression("S2: T+D"["agar2"]),
                                                       expression("S3: T+D"["pcr"]), expression("S4: T+I"["agar1"]), expression("S5: T+I"["agar2"]),
                                                       expression("S6: T+I"["pcr"]), expression("S7: Pre-emptive I"["all"]),
                                                       expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"])))+
  labs(subtitle = "(A) CRE") +
  theme(plot.subtitle = element_text(hjust = 0, size = 16, face="bold"), legend.position = "right", legend.box.just = "left", legend.justification = "left", legend.title = element_text(size = 15), legend.text = element_text(size = 10, hjust = 0)) +
  scale_y_continuous(name = "Annual mortality rate per\n100,000 hospital bed-days",labels = scales::label_number(accuracy = 1), breaks = seq(20, 80, by = 5))+
  scale_x_discrete(name = "Strategy") +
  geom_hline(yintercept = 57, color = "grey", linetype = "dashed", size = 1.1)+
  coord_cartesian(ylim = c(22.6, NA))

library(patchwork)
abc1 <- abc1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
# Combine the figures
combined_plot2x <- abc1 / abc2 + plot_layout(guides = 'collect') 
ggsave("pxox_deathrateHB.tiff",combined_plot2x, width = 11, height = 12, dpi = 1000)



####################################################################################################
####INFECTIONS WITH 95% CI below:################################################################################################
####################################################################################################
####################################################################################################
#Infections:
#B) Infections MRSS/
# For example:
dfxox <- data.frame(
  Strategy = paste("S", 0:15, sep=""),
  start_val = 4, # Replace these with your actual start values
  end_val = infecMRSA_a   # Replace these with your actual end values
)
dfxox <- na.omit(dfxox)
dfxox$StrategyNum <- 1:nrow(dfxox)
dfxox$Strategy <- factor(dfxox$Strategy, levels = unique(dfxox$Strategy))
###NEW FIGURES WITH 9 interventions only:
dfxoxf <- dfxox[!grepl("S7|S8|S9|S10|S11|S12", dfxox$Strategy), ]
dfxoxf$Strategy <- gsub("S13", "S7", dfxoxf$Strategy)
dfxoxf$Strategy <- gsub("S14", "S8", dfxoxf$Strategy)
dfxoxf$Strategy <- gsub("S15", "S9", dfxoxf$Strategy)
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 14] <- 8  # Replace '99' with your desired new value
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 15] <- 9
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 16] <- 10 
dfxoxf$ci_lower = dfxoxf$end_val - dif_CI_mrsaf
dfxoxf$ci_upper = dfxoxf$end_val + dif_CI_mrsaf
pxox_mrsa2<-ggplot(dfxoxf, aes(x = StrategyNum, ymin = start_val, ymax = end_val, fill = Strategy)) +
  geom_rect(aes(xmin = StrategyNum - 0.4, xmax = StrategyNum + 0.4), color = "black") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.4) +
  scale_fill_manual(values = lancet_colors, labels = c("S0: Do-nothing", expression("S1: T+D"["agar1"]), expression("S2: T+D"["agar2"]),
                                                       expression("S3: T+D"["pcr"]), expression("S4: T+I"["agar1"]), expression("S5: T+I"["agar2"]),
                                                       expression("S6: T+I"["pcr"]), expression("S7: Pre-emptive I"["all"]),
                                                       expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"])),
                    name = "Strategy") +
  scale_x_continuous(name = "Strategy", breaks = dfxoxf$StrategyNum, labels = dfxoxf$Strategy) +
  scale_y_continuous(name = "Hospital MRSA-infected population\nper 100 new admissions", limits = c(4, 17), breaks = seq(4, 17, by = 1)) +
  theme_minimal() +
  theme_lancet()+
  labs(subtitle = "(B) MRSA infected individuals") +
  theme(plot.subtitle = element_text(hjust = 0, size = 16, face="bold"), legend.position = "right", legend.box.just = "left", legend.justification = "left", legend.title = element_text(size = 12), legend.text = element_text(size = 10, hjust = 0)) +
  geom_hline(yintercept = infecMRSA_a[1]+0.02, color = "grey", linetype = "dashed", size = 1.1)+
  annotate("text", x = 1, y = Inf, label = "", fontface = "bold", size = 6, hjust = 1.1, vjust = 2) # Adjust 'x' and 'y' as needed, hjust and vjust to move outside
ggsave("pxox_mrsa2.tiff",pxox_mrsa2, width = 11, height = 7, dpi = 800)


#B) Infections CRE/
# For example:
dfxox <- data.frame(
  Strategy = paste("S", 0:15, sep=""),
  start_val = 10, # Replace these with your actual start values
  end_val = infecCRE_a   # Replace these with your actual end values
)
dfxox <- na.omit(dfxox)
dfxox$StrategyNum <- 1:nrow(dfxox)
dfxox$Strategy <- factor(dfxox$Strategy, levels = unique(dfxox$Strategy))
dfxoxf <- dfxox[!grepl("S7|S8|S9|S10|S11|S12", dfxox$Strategy), ]
dfxoxf$Strategy <- gsub("S13", "S7", dfxoxf$Strategy)
dfxoxf$Strategy <- gsub("S14", "S8", dfxoxf$Strategy)
dfxoxf$Strategy <- gsub("S15", "S9", dfxoxf$Strategy)
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 14] <- 8  # Replace '99' with your desired new value
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 15] <- 9
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 16] <- 10 

dfxoxf$ci_lower = dfxoxf$end_val - dif_CI_cref
dfxoxf$ci_upper = dfxoxf$end_val + dif_CI_cref
pxox_cre2<-ggplot(dfxoxf, aes(x = StrategyNum, ymin = start_val, ymax = end_val, fill = Strategy)) +
  geom_rect(aes(xmin = StrategyNum - 0.4, xmax = StrategyNum + 0.4), color = "black") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.4) +
  scale_fill_manual(values = lancet_colors, labels = c("S0: Do-nothing", expression("S1: T+D"["agar1"]), expression("S2: T+D"["agar2"]),
                                                       expression("S3: T+D"["pcr"]), expression("S4: T+I"["agar1"]), expression("S5: T+I"["agar2"]),
                                                       expression("S6: T+I"["pcr"]), expression("S7: Pre-emptive I"["all"]),
                                                       expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"])),
                    name = "Strategy") +
  scale_x_continuous(name = "Strategy", breaks = dfxoxf$StrategyNum, labels = dfxoxf$Strategy) +
  scale_y_continuous(name = "Hospital CRE-infected population\nper 100 new admissions", limits = c(10, 22), breaks = seq(10, 22, by = 1)) +
  theme_minimal() +
  theme_lancet()+
  labs(subtitle = "(A) CRE infected individuals") +
  theme(plot.subtitle = element_text(hjust = 0, size = 16, face="bold"), legend.position = "right", legend.box.just = "left", legend.justification = "left", legend.title = element_text(size = 12), legend.text = element_text(hjust = 0, size = 10)) +
  geom_hline(yintercept = infecCRE_a[1]+0.02, color = "grey", linetype = "dashed", size = 1.1)+
  annotate("text", x = 1, y = Inf, label = "", fontface = "bold", size = 6, hjust = 1.1, vjust = 2) # Adjust 'x' and 'y' as needed, hjust and vjust to move outside
ggsave("pxox_cre2.tiff",pxox_cre2, width = 11, height = 7, dpi = 800)

#Combinning two figures: 
library(patchwork)
pxox_cre_updated <- pxox_cre2 + 
  theme(axis.title.x = element_blank(),  # Removes the y-axis title
        axis.text.x = element_blank(),    # Removes the y-axis text/labels
        axis.ticks.y = element_blank())
pxox_cre2_updated <- pxox_cre_updated + theme(legend.position = "none")

# Combine the plots and add back a single legend




#-----------------------------------#
#TWO PANEL FIGURE here below.
#-----------------------------------#
four_panel_plotwx <- (pxox_cre2_updated) /
  (pxox_mrsa2 ) +
  plot_layout(guides = 'collect') 
# Display the combined plot
#print(four_panel_plot)
ggsave("four_panel_figure2.tiff",four_panel_plotwx, width = 11, height = 9, dpi = 1000)
#-----------------------------------#
#FOUR PANEL FIGURE here below.
#-----------------------------------#
pxox_cre2_updated <- pxox_cre2_updated+ theme(legend.position = "none")
abc1 <- abc1 + labs(subtitle = "(B) CRE-associated deaths") + theme(legend.position = "none")
pxox_mrsa2 <- pxox_mrsa2 +labs(subtitle = "(C) MRSA infected individuals") 
abc2 <- abc2 + labs(subtitle = "(D) MRSA-associated deaths") + theme(legend.position = "none")


combined_plot2x2 <- (pxox_cre2_updated |abc1 )/ (pxox_mrsa2 |abc2) + plot_layout(guides = 'collect') 
ggsave("pxox_fulldeathinfect_p4.tiff",combined_plot2x2, width = 12, height = 8, dpi = 1000)



