# Load required library
library(ggplot2)
library(deSolve)
library(reshape)
library(dplyr)
library(patchwork)
library(grid) 
library(reshape2)
library(ggpattern)
library(scales) 
library(RColorBrewer)
library(gridExtra) 
par(family = "Times New Roman")
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis/0_Figures_model")
#bar_colors <- c("grey0","#297b33","#3cb44b","#c4e8c9","#ffce19","#ffe954","#fff5b6", "#2746ba","#94a6e9","#c5cff3","#9a4624","#d18b3b","#e8c49b","#f36b0c","#ffa819","#fabb8f")
lancet_colors <- c("grey70","#76b5c5","#f28e79","#a1c181","#f7bb5f","#4ebdb6","#ff9da7","#f4d4a4","#b699d7","#88d8b0")
#lancet_colors <- c("grey70","#E41A1C","#ffa819","#FC8D62","#984EA3","#d12b90","#F781BF", "#A65628","#d18b3b","#E5C494","#4DAF4A","#66C2A5","#A6D854","#2746ba","#8DA0CB","#377EB8" )

# Define a custom theme
theme_lancet <- function(base_size = 14, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(size = base_size, color = "black"),
      plot.title = element_text(size = base_size * 1.1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0, face = "bold"),
      axis.title = element_text(size = base_size),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text = element_text(size = base_size * 0.8),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = base_size * 0.8),
      legend.title = element_text(size = base_size, face = "bold"),
      panel.grid.major = element_line(color = "#eaeaea"), #"#eaeaea"
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

######################################################################
######################################################################
#### Figure 1. Dynamics MRSA CRE over time
######################################################################
######################################################################
#bar_colors <- c("grey","#297b33","#3cb44b","#c4e8c9","#ffce19","#ffe954","#fff5b6", "#2746ba","#94a6e9","#c5cff3","#9a4624","#d18b3b","#e8c49b","#f36b0c","#ffa819","#fabb8f")
#Graphs MRSA transmissions, Deaths MRSA, and MRSA infections PER STRATEGY and in /100 bed-days.
time <- 1:366
id <- 1:366
results_matrix <- as.data.frame(cbind(id, results_epi_MRSAprev))
results_matrix[1:200,3:17]<- NA
results_matrix[199:200,3:17]<- results_matrix[199:200,2]
colnames(results_matrix) <- c("id", "S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_epi_MRSAprev_l <- melt(results_matrix, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))
#FIRST FIGURE
par(family = "Times New Roman")
library(ggplot2)

#---------------------------------------------------------#

#FILTER INFORMATION TO ONLY 9 interventions:####
results_epi_MRSAprev_lf <- results_epi_MRSAprev_l[!grepl("S7|S8|S9|S10|S11|S12", results_epi_MRSAprev_l$variable), ]
results_epi_MRSAprev_lf$variable <- gsub("S13", "S7", results_epi_MRSAprev_lf$variable)
results_epi_MRSAprev_lf$variable <- gsub("S14", "S8", results_epi_MRSAprev_lf$variable)
results_epi_MRSAprev_lf$variable <- gsub("S15", "S9", results_epi_MRSAprev_lf$variable)

#MRSA COLONISED
MRSA_colonisation_strategies <- ggplot(results_epi_MRSAprev_lf, aes(x = id, y = value, colour = variable)) +
  geom_line(linewidth = 2) +
  scale_color_manual(values = lancet_colors, 
                     labels = c("S0: Do-nothing", expression(S1: T+D["agar1"]), expression(S2: T+D["agar2"]), expression(S3: T+D["pcr"]), expression(S4: T+I["agar1"]), expression(S5: T+I["agar2"]), expression(S6: T+I["pcr"]),  expression("S7: Pre-emptive I"["all"]), expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"]))) +
  labs(x = "Time (days)", 
       y = "Number of MRSA colonised individuals", 
       title = "(C)",
       subtitle = "",
       colour = "Strategy") +
  scale_x_continuous(breaks = seq(170, 365, by = 20), limits = c(170, 365)) +
  scale_y_continuous(breaks = seq(45, 115, by = 10), limits = c(45, 115)) +
  geom_hline(yintercept = 45, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 170, linetype = "solid", color = "grey", size = 0.5) +
  theme_lancet()+
  theme(legend.position = "right",
        legend.justification = "left",
        legend.title = element_text(size = 16),  # Increase size of legend title
        legend.text = element_text(hjust = 0, size = 14),  # Increase size of legend text
        axis.title.x = element_text(size = 21),  # Increase size of x-axis title
        axis.title.y = element_text(size = 21),  # Increase size of y-axis title
        axis.text.x = element_text(size = 19),  # Increase size of x-axis text/labels
        axis.text.y = element_text(size = 19),  # Increase size of y-axis text/labels
        plot.title = element_text(size = 16),   # Increase size of main title
        plot.subtitle = element_text(size = 14))  # Increase size of subtitle
# Print the plot
#print(MRSA_colonisation_strategies)
ggsave("MRSA_colonisation_strategiess.tiff", MRSA_colonisation_strategies, width = 11, height = 7, dpi = 800)

#MRSA INFECTED
results_matrix <- as.data.frame(cbind(id, results_epi_MRSAinfe))
results_matrix[1:200,3:17]<- NA
results_matrix[199:200,3:17]<- results_matrix[199:200,2]
colnames(results_matrix) <- c("id","S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_epi_MRSAinfe_l <- melt(results_matrix, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))
#Filter information:
results_epi_MRSAinfe_lf <- results_epi_MRSAinfe_l[!grepl("S7|S8|S9|S10|S11|S12", results_epi_MRSAinfe_l$variable), ]
results_epi_MRSAinfe_lf$variable <- gsub("S13", "S7", results_epi_MRSAinfe_lf$variable)
results_epi_MRSAinfe_lf$variable <- gsub("S14", "S8", results_epi_MRSAinfe_lf$variable)
results_epi_MRSAinfe_lf$variable <- gsub("S15", "S9", results_epi_MRSAinfe_lf$variable)

MRSA_inf_strategies<-ggplot(results_epi_MRSAinfe_lf, aes(x = id, y = value, colour = variable)) + 
  geom_line(linewidth = 2) +
  scale_color_manual(values = lancet_colors, 
                     labels = c("S0: Do-nothing", expression(S1: T+D["agar1"]), expression(S2: T+D["agar2"]), expression(S3: T+D["pcr"]), expression(S4: T+I["agar1"]), expression(S5: T+I["agar2"]), expression(S6: T+I["pcr"]),  expression("S7: Pre-emptive I"["all"]), expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"]))) +
  labs(x = "Time (days)", 
       y = "Number of MRSA infected individuals", 
       title = "(A)",
       subtitle = "",
       colour = "Strategy   ") +
  scale_x_continuous(breaks = seq(170, 365, by = 20), limits = c(170, 365)) +
  scale_y_continuous(breaks = seq(8, 20, by = 1), limits = c(8, 20)) +
  geom_hline(yintercept = 8, linetype = "solid", color = "grey", size = 0.5) +  # Add horizontal line at y=75
  geom_vline(xintercept = 170, linetype = "solid", color = "grey", size = 0.5) +
  theme_lancet() +
  theme(legend.position = "right",
        legend.justification = "left",
        legend.title = element_text(size = 16),  # Increase size of legend title
        legend.text = element_text(hjust = 0, size = 14),  # Increase size of legend text
        axis.title.x = element_text(size = 21),  # Increase size of x-axis title
        axis.title.y = element_text(size = 21),  # Increase size of y-axis title
        axis.text.x = element_text(size = 19),  # Increase size of x-axis text/labels
        axis.text.y = element_text(size = 19),  # Increase size of y-axis text/labels
        plot.title = element_text(size = 16),   # Increase size of main title
        plot.subtitle = element_text(size = 14))  # Increase size of subtitle
ggsave("MRSA_inf_strategies.tiff", MRSA_inf_strategies, width = 11, height = 7, dpi = 800)
#  annotate("text", x = 150, y = 115, label = "(A)", size = 5, fontface = "bold", colour = "black",
#           hjust = 1.8, vjust = -1.2)


#CRE:
results_matrix <- as.data.frame(cbind(id, results_epi_CREprev))
results_matrix[1:200,3:17]<- NA
results_matrix[199:200,3:17]<- results_matrix[199:200,2]
colnames(results_matrix) <- c("id", "S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_epi_CREprev_l <- melt(results_matrix, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))
results_epi_CREprev_l$value[results_epi_CREprev_l$variable == "S6" & results_epi_CREprev_l$id > 200] <- results_epi_CREprev_l$value[results_epi_CREprev_l$variable == "S6" & results_epi_CREprev_l$id > 200] + 1

#CRE COLONISED
#FILTER INFORMATION TO ONLY 9 interventions:####
results_epi_CREprev_lf <- results_epi_CREprev_l[!grepl("S7|S8|S9|S10|S11|S12", results_epi_CREprev_l$variable), ]
results_epi_CREprev_lf$variable <- gsub("S13", "S7", results_epi_CREprev_lf$variable)
results_epi_CREprev_lf$variable<- gsub("S14", "S8", results_epi_CREprev_lf$variable)
results_epi_CREprev_lf$variable <- gsub("S15", "S9", results_epi_CREprev_lf$variable)

CRE_colon_strategies<- ggplot(results_epi_CREprev_lf, aes(x = id, y = value, colour = variable)) +
  geom_line(linewidth = 2) +
  scale_color_manual(values = lancet_colors, 
                     labels = c("S0: Do-nothing", expression(S1: T+D["agar1"]), expression(S2: T+D["agar2"]), expression(S3: T+D["pcr"]), expression(S4: T+I["agar1"]), expression(S5: T+I["agar2"]), expression(S6: T+I["pcr"]),  expression("S7: Pre-emptive I"["all"]), expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"]))) +
  labs(x = "Time (days)", 
       y = "Number of CRE colonised individuals", 
       title = "(B)",
       subtitle = "",
       colour = "Strategy") +
  scale_x_continuous(breaks = seq(170, 365, by = 20), limits = c(170, 365)) +
  scale_y_continuous(breaks = seq(100, 150, by = 10), limits = c(100, 150)) +
  geom_hline(yintercept = 100, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 170, linetype = "solid", color = "grey", size = 0.5) +
  theme_lancet()+
  theme(legend.position = "right",
        legend.justification = "left",
        legend.title = element_text(size = 16),  # Increase size of legend title
        legend.text = element_text(hjust = 0, size = 14),  # Increase size of legend text
        axis.title.x = element_text(size = 21),  # Increase size of x-axis title
        axis.title.y = element_text(size = 21),  # Increase size of y-axis title
        axis.text.x = element_text(size = 19),  # Increase size of x-axis text/labels
        axis.text.y = element_text(size = 19),  # Increase size of y-axis text/labels
        plot.title = element_text(size = 16),   # Increase size of main title
        plot.subtitle = element_text(size = 14))  # Increase size of subtitle
ggsave("CRE_colon_strategies.tiff",CRE_colon_strategies, width = 11, height = 7, dpi = 800)

#CRE INFECTED:
results_matrix <- as.data.frame(cbind(id, results_epi_CREinfe))
results_matrix[1:200,3:17]<- NA
results_matrix[199:200,3:17]<- results_matrix[199:200,2]
colnames(results_matrix) <- c("id", "S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_epi_CREinfe_l <- melt(results_matrix, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))
results_epi_CREinfe_l$value[results_epi_CREinfe_l$variable == "S6" & results_epi_CREinfe_l$id > 200] <- results_epi_CREinfe_l$value[results_epi_CREinfe_l$variable == "S6" & results_epi_CREinfe_l$id > 200] + 0.15

#FILTER INFORMATION TO ONLY 9 interventions:####
results_epi_CREinfe_lf <- results_epi_CREinfe_l[!grepl("S7|S8|S9|S10|S11|S12", results_epi_CREinfe_l$variable), ]
results_epi_CREinfe_lf$variable <- gsub("S13", "S7", results_epi_CREinfe_lf$variable)
results_epi_CREinfe_lf$variable<- gsub("S14", "S8", results_epi_CREinfe_lf$variable)
results_epi_CREinfe_lf$variable <- gsub("S15", "S9", results_epi_CREinfe_lf$variable)

CRE_inf_strategies<- ggplot(results_epi_CREinfe_lf, aes(x = id, y = value, colour = variable)) +
  geom_line(linewidth = 2) +
  scale_color_manual(values = lancet_colors, 
                     labels = c("S0: Do-nothing", expression(S1: T+D["agar1"]), expression(S2: T+D["agar2"]), expression(S3: T+D["pcr"]), expression(S4: T+I["agar1"]), expression(S5: T+I["agar2"]), expression(S6: T+I["pcr"]),  expression("S7: Pre-emptive I"["all"]), expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"]))) +
  labs(x = "Time (days)", 
       y = "Number of CRE infected individuals", 
       title = "(B)",
       subtitle = "",
       colour = "Strategy") +
  scale_x_continuous(breaks = seq(170, 365, by = 20), limits = c(170, 365)) +
  scale_y_continuous(breaks = seq(20, 28, by = 1), limits = c(20, 28)) +
  geom_hline(yintercept = 20, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 170, linetype = "solid", color = "grey", size = 0.5) +
  theme_lancet()+
  theme(legend.position = "right",
        legend.justification = "left",
        legend.title = element_text(size = 16),  # Increase size of legend title
        legend.text = element_text(hjust = 0, size = 14),  # Increase size of legend text
        axis.title.x = element_text(size = 21),  # Increase size of x-axis title
        axis.title.y = element_text(size = 21),  # Increase size of y-axis title
        axis.text.x = element_text(size = 19),  # Increase size of x-axis text/labels
        axis.text.y = element_text(size = 19),  # Increase size of y-axis text/labels
        plot.title = element_text(size = 16),   # Increase size of main title
        plot.subtitle = element_text(size = 14))  # Increase size of subtitle
ggsave("CRE_inf_strategies.tiff",CRE_inf_strategies, width = 11, height = 7, dpi = 800)


#---------------------------------------------------------#
#FOUR PANEL FIGURE OF THE DYNAMICS!
# Arrange the plots in a 2x2 grid
adjusted_theme <- theme(
  plot.margin = unit(c(1, 1, 1, 1), "mm"), # Reduce plot margins
  plot.background = element_blank(),
  panel.background = element_blank(),
  legend.background = element_blank(),
  panel.spacing = unit(0.5, "lines"), # Reduce spacing between panels
  strip.background = element_blank())

MRSA_colonisation_strategiesF <- MRSA_colonisation_strategies +geom_line(linewidth = 2.5)+ labs(title = "(B) MRSA colonised individuals") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 16),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x= element_blank())+adjusted_theme

MRSA_inf_strategiesF <- MRSA_inf_strategies +geom_line(linewidth = 2.5) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 15),  # Adjust x-axis title size
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15)) + labs(title = "(D) MRSA infected individuals")+adjusted_theme
CRE_colon_strategiesF <- CRE_colon_strategies +geom_line(linewidth = 2.5) + labs(title = "(A) CRE colonised individuals") +   
  theme(axis.title.x = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),
        plot.title = element_text(size = 16),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +adjusted_theme
CRE_inf_strategiesF <- CRE_inf_strategies +geom_line(linewidth = 2.5)+ labs(title = "(C) CRE infected individuals") + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 16),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+adjusted_theme

CRE_dead_strategiesF <- CRE_dead_strategies +geom_line(linewidth = 2.5)+ labs(title = "(C) CRE-associated deaths") + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 16),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+adjusted_theme

MRSA_dead_strategiesF <-MRSA_dead_strategies +geom_line(linewidth = 2.5)+ labs(title = "(D) MRSA-associated deaths") + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 16),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+adjusted_theme

#---------------------------------------------------------#
# Combine the plots and ensure equal sizes 
#---------------------------------------------------------# 
combined_MRSACRE_4pan <- (CRE_colon_strategiesF | MRSA_colonisation_strategiesF )  / 
  (CRE_inf_strategiesF | MRSA_inf_strategiesF) & 
  plot_layout(guides = 'collect', 
              heights = c(1, 1), # Equal height for both rows
              widths = c(1, 1))  +
  theme_void()# Equal width for both columns
# Save the combined plot with adjusted sizes
legend <- get_legend(MRSA_colonisation_strategies)
final_layout <- combined_MRSACRE_4pan | legend 
final_layout <- final_layout + plot_layout(widths = c(5, 1))
# Adjust as necessary for size
ggsave("combined_MRSACRE_4pan_with_legend.tiff", final_layout, width = 15, height = 9, dpi = 1000)



#New figures replacing colonisatins with deaths
MRSA_inf_strategiesF <- MRSA_inf_strategies +geom_line(linewidth = 2.5) + 
  theme(plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 15),  # Adjust x-axis title size
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15)) + labs(title = "(B) MRSA infected individuals")+adjusted_theme
CRE_inf_strategiesF <- CRE_inf_strategies +geom_line(linewidth = 2.5)+ labs(title = "(A) CRE infected individuals") + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15),
        plot.title = element_text(size = 16),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))+adjusted_theme
combined_MRSACRE_4pan2 <- (CRE_inf_strategiesF | MRSA_inf_strategiesF)  / 
  (CRE_dead_strategiesF | MRSA_dead_strategiesF )  & 
  plot_layout(guides = 'collect', 
              heights = c(1, 1), # Equal height for both rows
              widths = c(1, 1))  +
  theme_void()# Equal width for both columns
# Save the combined plot with adjusted sizes
legend <- get_legend(MRSA_colonisation_strategies)
final_layout2 <- combined_MRSACRE_4pan2 | legend 
final_layout2 <- final_layout2 + plot_layout(widths = c(5, 1))

ggsave("combined_MRSACRE_4pan_with_legend_NEWdeaths.tiff", final_layout2, width = 15, height = 9, dpi = 1000)


combined_MRSACRE_2pan2infec <- (CRE_inf_strategiesF )  /(MRSA_inf_strategiesF)  + 
  plot_layout(guides = 'collect')  
  theme_void()# Equal width for both columns

ggsave("combined_MRSACRE_infec.tiff",combined_MRSACRE_2pan2infec, width = 12, height = 15, dpi = 800)


combined_MRSACRE_death <- CRE_dead_strategies /MRSA_dead_strategies # Stack them vertically
ggsave("combined_MRSACRE_death.tiff",combined_MRSACRE_death, width = 12, height = 15, dpi = 800)
#

#CRE COLONISATION

MRSA_colonisation_strategiesF <- MRSA_colonisation_strategies +geom_line(linewidth = 2.6)+ labs(title = "(B) MRSA colonised individuals") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 16),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15))+adjusted_theme
CRE_colon_strategiesF <- CRE_colon_strategies +geom_line(linewidth = 2.6) + labs(title = "(A) CRE colonised individuals") +   
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 16),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +adjusted_theme
combined_MRSACRE_col <- CRE_colon_strategiesF /MRSA_colonisation_strategiesF + plot_layout(guides = 'collect')# Stack them vertically
ggsave("combined_MRSACRE_COL.tiff",combined_MRSACRE_col, width = 10, height = 12, dpi = 800)


######################################################################
######################################################################
#### Figure 2. Hospital-colonised population per 100 new admissions
######################################################################
######################################################################

# Assuming you have a data frame where 'start_val' and 'end_val' are your starting and ending values for each 'Strategy'
# For example:
dfxox <- data.frame(
  Strategy = paste("S", 0:15, sep=""),
  start_val = 25, # Replace these with your actual start values
  end_val = colonMRSA_a   # Replace these with your actual end values
)
dfxox <- na.omit(dfxox)
dfxox$StrategyNum <- 1:nrow(dfxox)
# Now, let's create a floating bar plot using geom_segment
library(ggplot2)
dfxox$Strategy <- factor(dfxox$Strategy, levels = unique(dfxox$Strategy))
###NEW FIGURES WITH 9 interventions only:
dfxoxf <- dfxox[!grepl("S7|S8|S9|S10|S11|S12", dfxox$Strategy), ]
dfxoxf$Strategy <- gsub("S13", "S7", dfxoxf$Strategy)
dfxoxf$Strategy <- gsub("S14", "S8", dfxoxf$Strategy)
dfxoxf$Strategy <- gsub("S15", "S9", dfxoxf$Strategy)
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 14] <- 8  # Replace '99' with your desired new value
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 15] <- 9
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 16] <- 10 

#arrows(bp1, lower_boundCRE_colf, bp1, Upper_boundCRE_colf, angle = 90, code = 3, length = 0.05)

# Replot
pxox_mrsa<-ggplot(dfxoxf, aes(x = StrategyNum, ymin = start_val, ymax = end_val, fill = Strategy)) +
  geom_rect(aes(xmin = StrategyNum - 0.4, xmax = StrategyNum + 0.4), color = "black") +
  scale_fill_manual(values = lancet_colors, labels = c("S0: Do-nothing", expression("S1: T+D"["agar1"]), expression("S2: T+D"["agar2"]),
                                                       expression("S3: T+D"["pcr"]), expression("S4: T+I"["agar1"]), expression("S5: T+I"["agar2"]),
                                                       expression("S6: T+I"["pcr"]), expression("S7: Pre-emptive I"["all"]),
                                                       expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"])),
                    name = "Strategy") +
  scale_x_continuous(name = "Strategy", breaks = dfxoxf$StrategyNum, labels = dfxoxf$Strategy) +
  scale_y_continuous(name = "Hospital MRSA-colonised population\nper 100 new admissions", limits = c(25, 65), breaks = seq(25, 65, by = 5)) +
  theme_minimal() +
  theme_lancet()+
  labs(subtitle = "(B) MRSA colonised individuals") +
  theme(plot.subtitle = element_text(hjust = 0, size = 16, face="bold"), legend.position = "right", legend.box.just = "left", legend.justification = "left", legend.title = element_text(size = 12), legend.text = element_text(size = 10, hjust = 0)) +
  geom_hline(yintercept = colonMRSA_a[1]+0.08, color = "grey", linetype = "dashed", size = 1.1)+
  annotate("text", x = 1, y = Inf, label = "", fontface = "bold", size = 6, hjust = 1.1, vjust = 2)
ggsave("pxox_mrsa.tiff",pxox_mrsa, width = 11, height = 7, dpi = 800)

#B) Infections MRSS/
# For example:
dfxox <- data.frame(
  Strategy = paste("S", 0:15, sep=""),
  start_val = 6, # Replace these with your actual start values
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

pxox_mrsa2<-ggplot(dfxoxf, aes(x = StrategyNum, ymin = start_val, ymax = end_val, fill = Strategy)) +
  geom_rect(aes(xmin = StrategyNum - 0.4, xmax = StrategyNum + 0.4), color = "black") +
  scale_fill_manual(values = lancet_colors, labels = c("S0: Do-nothing", expression("S1: T+D"["agar1"]), expression("S2: T+D"["agar2"]),
                                                       expression("S3: T+D"["pcr"]), expression("S4: T+I"["agar1"]), expression("S5: T+I"["agar2"]),
                                                       expression("S6: T+I"["pcr"]), expression("S7: Pre-emptive I"["all"]),
                                                       expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"])),
                    name = "Strategy") +
  scale_x_continuous(name = "Strategy", breaks = dfxoxf$StrategyNum, labels = dfxoxf$Strategy) +
  scale_y_continuous(name = "Hospital MRSA-infected population\nper 100 new admissions", limits = c(6, 13), breaks = seq(6, 13, by = 1)) +
  theme_minimal() +
  theme_lancet()+
  labs(subtitle = "(D) MRSA infected individuals") +
  theme(plot.subtitle = element_text(hjust = 0, size = 16, face="bold"), legend.position = "right", legend.box.just = "left", legend.justification = "left", legend.title = element_text(size = 12), legend.text = element_text(size = 10, hjust = 0)) +
  geom_hline(yintercept = infecMRSA_a[1]+0.02, color = "grey", linetype = "dashed", size = 1.1)+
  annotate("text", x = 1, y = Inf, label = "", fontface = "bold", size = 6, hjust = 1.1, vjust = 2) # Adjust 'x' and 'y' as needed, hjust and vjust to move outside
ggsave("pxox_mrsa2.tiff",pxox_mrsa2, width = 11, height = 7, dpi = 800)

#Combinning two figures: 
library(patchwork)
pxox_mrsa_updated <- pxox_mrsa + 
  theme(axis.title.x = element_blank(),  # Removes the y-axis title
        axis.text.x = element_blank(),    # Removes the y-axis text/labels
        axis.ticks.y = element_blank())
pxox_mrsa2_updated <- pxox_mrsa2 + theme(legend.position = "none")
# Combine the plots and add back a single legend
combined_plot20x2 <- pxox_mrsa_updated / pxox_mrsa2_updated + plot_layout(guides = 'collect') 
ggsave("combined_plot20x_mrsa.tiff",combined_plot20x2, width = 11, height = 8, dpi = 800)
# View the combined plot
combined_plot20x2

#------------------------------------#
#CRE:

dfxox <- data.frame(
  Strategy = paste("S", 0:15, sep=""),
  start_val = 55, # Replace these with your actual start values
  end_val = colonCRE_a   # Replace these with your actual end values
)
dfxox <- na.omit(dfxox)
dfxox$StrategyNum <- 1:nrow(dfxox)
# Now, let's create a floating bar plot using geom_segment
library(ggplot2)
dfxox$Strategy <- factor(dfxox$Strategy, levels = unique(dfxox$Strategy))
dfxoxf <- dfxox[!grepl("S7|S8|S9|S10|S11|S12", dfxox$Strategy), ]
dfxoxf$Strategy <- gsub("S13", "S7", dfxoxf$Strategy)
dfxoxf$Strategy <- gsub("S14", "S8", dfxoxf$Strategy)
dfxoxf$Strategy <- gsub("S15", "S9", dfxoxf$Strategy)
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 14] <- 8  # Replace '99' with your desired new value
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 15] <- 9
dfxoxf$StrategyNum[dfxoxf$StrategyNum == 16] <- 10 
# Replot
pxox_cre<-ggplot(dfxoxf, aes(x = StrategyNum, ymin = start_val, ymax = end_val, fill = Strategy)) +
  geom_rect(aes(xmin = StrategyNum - 0.4, xmax = StrategyNum + 0.4), color = "black") +
  scale_fill_manual(values = lancet_colors, labels = c("S0: Do-nothing", expression("S1: T+D"["agar1"]), expression("S2: T+D"["agar2"]),
                                                       expression("S3: T+D"["pcr"]), expression("S4: T+I"["agar1"]), expression("S5: T+I"["agar2"]),
                                                       expression("S6: T+I"["pcr"]), expression("S7: Pre-emptive I"["all"]),
                                                       expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"])),
                    name = "Strategy") +
  scale_x_continuous(name = "Strategy", breaks = dfxoxf$StrategyNum, labels = dfxoxf$Strategy) +
  scale_y_continuous(name = "Hospital CRE-colonised population\nper 100 new admissions", limits = c(55, 85), breaks = seq(55, 85, by = 5)) +
  theme_minimal() +
  theme_lancet()+
  labs(subtitle = "(A) CRE colonised individuals") +
  theme(plot.subtitle = element_text(hjust = 0, size = 16, face="bold"), legend.position = "right", legend.box.just = "left", legend.justification = "left", legend.title = element_text(size = 12), legend.text = element_text(size = 10, hjust = 0)) +
  geom_hline(yintercept = colonCRE_a[1]+0.08, color = "grey", linetype = "dashed", size = 1.1)+
  annotate("text", x = 1, y = Inf, label = "", fontface = "bold", size = 6, hjust = 1.1, vjust = 2) # Adjust 'x' and 'y' as needed, hjust and vjust to move outside
ggsave("pxox_cre.tiff",pxox_cre, width = 11, height = 7, dpi = 800)

#B) Infections CRE/
# For example:
dfxox <- data.frame(
  Strategy = paste("S", 0:15, sep=""),
  start_val = 13, # Replace these with your actual start values
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
pxox_cre2<-ggplot(dfxoxf, aes(x = StrategyNum, ymin = start_val, ymax = end_val, fill = Strategy)) +
  geom_rect(aes(xmin = StrategyNum - 0.4, xmax = StrategyNum + 0.4), color = "black") +
  scale_fill_manual(values = lancet_colors, labels = c("S0: Do-nothing", expression("S1: T+D"["agar1"]), expression("S2: T+D"["agar2"]),
                                                       expression("S3: T+D"["pcr"]), expression("S4: T+I"["agar1"]), expression("S5: T+I"["agar2"]),
                                                       expression("S6: T+I"["pcr"]), expression("S7: Pre-emptive I"["all"]),
                                                       expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"])),
                    name = "Strategy") +
  scale_x_continuous(name = "Strategy", breaks = dfxoxf$StrategyNum, labels = dfxoxf$Strategy) +
  scale_y_continuous(name = "Hospital CRE-infected population\nper 100 new admissions", limits = c(13, 18.2), breaks = seq(13, 18, by = 1)) +
  theme_minimal() +
  theme_lancet()+
  labs(subtitle = "(C) CRE infected individuals") +
  theme(plot.subtitle = element_text(hjust = 0, size = 16, face="bold"), legend.position = "right", legend.box.just = "left", legend.justification = "left", legend.title = element_text(size = 12), legend.text = element_text(hjust = 0, size = 10)) +
  geom_hline(yintercept = infecCRE_a[1]+0.02, color = "grey", linetype = "dashed", size = 1.1)+
  annotate("text", x = 1, y = Inf, label = "", fontface = "bold", size = 6, hjust = 1.1, vjust = 2) # Adjust 'x' and 'y' as needed, hjust and vjust to move outside
ggsave("pxox_cre2.tiff",pxox_cre2, width = 11, height = 7, dpi = 800)

#Combinning two figures: 
library(patchwork)
pxox_cre_updated <- pxox_cre + 
  theme(axis.title.x = element_blank(),  # Removes the y-axis title
        axis.text.x = element_blank(),    # Removes the y-axis text/labels
        axis.ticks.y = element_blank())
pxox_cre2_updated <- pxox_cre2 + theme(legend.position = "none")

# Combine the plots and add back a single legend
combined_plot20x <- pxox_cre_updated / pxox_cre2_updated + plot_layout(guides = 'collect') 
ggsave("combined_plot20x.tiff",combined_plot20x, width = 11, height = 8, dpi = 800)
# View the combined plot
combined_plot20x



#-----------------------------------#
#FOUR PANEL FIGURE here below.
#-----------------------------------#
four_panel_plot <- (pxox_cre_updated| pxox_mrsa_updated ) /
   (pxox_cre2_updated |pxox_mrsa2_updated) +
   plot_layout(guides = 'collect') 
 # Display the combined plot
 #print(four_panel_plot)
 ggsave("four_panel_figure2.tiff",four_panel_plot, width = 15, height = 8, dpi = 1000)

 
 
 
 ######-----------------------------------######
 ####DEAD CRE:######
 ######-----------------------------------######
 dfxox <- data.frame(
   Strategy = paste("S", 0:15, sep=""),
   start_val = 0, # Replace these with your actual start values
   end_val = DeadCRE_a   # Replace these with your actual end values
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
 pxox_cre2a<-ggplot(dfxoxf, aes(x = StrategyNum, ymin = start_val, ymax = end_val, fill = Strategy)) +
   geom_rect(aes(xmin = StrategyNum - 0.4, xmax = StrategyNum + 0.4), color = "black") +
   scale_fill_manual(values = lancet_colors, labels = c("S0: Do-nothing", expression("S1: T+D"["agar1"]), expression("S2: T+D"["agar2"]),
                                                        expression("S3: T+D"["pcr"]), expression("S4: T+I"["agar1"]), expression("S5: T+I"["agar2"]),
                                                        expression("S6: T+I"["pcr"]), expression("S7: Pre-emptive I"["all"]),
                                                        expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"])),
                     name = "Strategy") +
   scale_x_continuous(name = "Strategy", breaks = dfxoxf$StrategyNum, labels = dfxoxf$Strategy) +
   scale_y_continuous(name = "Hospital CRE-infected population\nper 100 new admissions", limits = c(0, 10), breaks = seq(0, 10, by = 1)) +
   theme_minimal() +
   theme_lancet()+
   labs(subtitle = "(C) CRE infected individuals") +
   theme(plot.subtitle = element_text(hjust = 0, size = 16, face="bold"), legend.position = "right", legend.box.just = "left", legend.justification = "left", legend.title = element_text(size = 12), legend.text = element_text(hjust = 0, size = 10)) +
   geom_hline(yintercept = DeadCRE_a[1]+0.02, color = "grey", linetype = "dashed", size = 1.1)+
   annotate("text", x = 1, y = Inf, label = "", fontface = "bold", size = 6, hjust = 1.1, vjust = 2) # Adjust 'x' and 'y' as needed, hjust and vjust to move outside
 ggsave("pxox_cre2-deaths.tiff",pxox_cre2a, width = 11, height = 7, dpi = 800)
 
 #Combinning two figures: 
 library(patchwork)
 pxox_cre_updateda <- pxox_crea + 
   theme(axis.title.x = element_blank(),  # Removes the y-axis title
         axis.text.x = element_blank(),    # Removes the y-axis text/labels
         axis.ticks.y = element_blank())
 pxox_cre2_updateda <- pxox_cre2 + theme(legend.position = "none")
 

 
 
 
################################################################################
################################################################################
#.   Figure SM: for health gains and cost per admisison
################################################################################
################################################################################
 
qaly_padmin <- matrix(nrow =1, ncol =  length(model_list))
cost_padmin <- matrix(nrow =1, ncol =  length(model_list))
qaly_padmin[1,1] <- 0
cost_padmin[1,1] <- 0
for (i in 2:length(model_list)) {
  qaly_padmin[1, i]<- ((results_icerM[3, i])- results_icerM[3, 1])/(results_epi_acM[6, i])
  cost_padmin[1, i]<-  (results_econM[1, i]- results_econM[1, 1])/(results_epi_acM[6, i])  
}
new_df_ecD<- as.data.frame(rbind(qaly_padmin, cost_padmin))
colnames(new_df_ecD) <- c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
#rownames(new_df_ecD) <- c("QALYs","Cost")
new_df_ecD$Category = c("QALYs", "Costs")
library(tidyr)
long_df_ecD <- pivot_longer(new_df_ecD, 
                            cols = -Category, 
                            names_to = "Strategy", 
                            values_to = "Value")
long_df_ecD <- long_df_ecD %>% 
  pivot_wider(names_from = Category, values_from = Value)
#PLOT
long_df_ecDf <- long_df_ecD[!grepl("S7|S8|S9|S10|S11|S12", long_df_ecD$Strategy), ]
long_df_ecDf$Strategy <- gsub("S13", "S7", long_df_ecDf$Strategy)
long_df_ecDf$Strategy <- gsub("S14", "S8", long_df_ecDf$Strategy)
long_df_ecDf$Strategy <- gsub("S15", "S9", long_df_ecDf$Strategy)

newAdmgains_mrsa<-ggplot(long_df_ecDf, aes(x = QALYs, y = Costs, fill = Strategy)) +  
  geom_point(shape = 21, size = 4, colour = "black", stroke = 1) +
  #geom_errorbar(aes(ymin = Costs_lower, ymax = Costs_upper), width = 0.002, colour = "black") +
  #geom_errorbarh(aes(xmin = QALYs_lower, xmax = QALYs_upper), height = 0.2, colour = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey") +
  scale_fill_manual(  # Use 'scale_fill_manual' instead of 'scale_color_manual'
    values = lancet_colors,
    labels = c("S0: Do-nothing", expression("S1: T+D"[agar1]), expression("S2: T+D"[agar2]), expression("S3: T+D"[pcr]), expression("S4: T+I"[agar1]), expression("S5: T+I"[agar2]), expression("S6: T+I"[pcr]), expression("S7: Pre-emptive I"[all]), expression("S8: Pre-emptive I"[men]), expression("S9: Pre-emptive I"[women])),
    name = "Strategy") +
  scale_y_continuous(breaks = seq(0, 40, by = 2))+
  scale_x_continuous(breaks = seq(0, 0.03, by = 0.002))+
  labs(title = "",
       x = "Health benefit per new admission (∆QALYs)",
       y = "Incremental cost per new admission (∆$)",
       subtitle = "(B) MRSA model") +
  theme_minimal() +
  theme_lancet() +
  theme(legend.position = "right",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10, hjust = 0),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
ggsave("newAdmgains_mrsa.tiff", newAdmgains_mrsa, width = 11, height = 7, dpi = 800)

#NEW BITS!
qaly_padmin <- matrix(nrow =1, ncol =  length(model_list))
cost_padmin <- matrix(nrow =1, ncol =  length(model_list))
qaly_padmin[1,1] <- 0
cost_padmin[1,1] <- 0
for (i in 2:length(model_list)) {
  qaly_padmin[1, i]<- ((results_icerC[3, i])- results_icerC[3, 1])/(results_epi_acC[6, i])
  cost_padmin[1, i]<-  (results_icerC[2, i]- results_icerC[2, 1])/(results_epi_acC[6, i])  
}
new_df_ecD<- as.data.frame(rbind(qaly_padmin, cost_padmin))
colnames(new_df_ecD) <- c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
#rownames(new_df_ecD) <- c("QALYs","Cost")
new_df_ecD$Category = c("QALYs", "Costs")
library(tidyr)
long_df_ecD <- pivot_longer(new_df_ecD, 
                            cols = -Category, 
                            names_to = "Strategy", 
                            values_to = "Value")
long_df_ecD <- long_df_ecD %>% 
  pivot_wider(names_from = Category, values_from = Value)
long_df_ecDf <- long_df_ecD[!grepl("S7|S8|S9|S10|S11|S12", long_df_ecD$Strategy), ]
long_df_ecDf$Strategy <- gsub("S13", "S7", long_df_ecDf$Strategy)
long_df_ecDf$Strategy <- gsub("S14", "S8", long_df_ecDf$Strategy)
long_df_ecDf$Strategy <- gsub("S15", "S9", long_df_ecDf$Strategy)
#PLOT
newAdmgains_cre<-ggplot(long_df_ecDf, aes(x = QALYs, y = Costs, fill = Strategy)) +  
  geom_point(shape = 21, size = 4, colour = "black", stroke = 1) +
  #geom_errorbar(aes(ymin = Costs_lower, ymax = Costs_upper), width = 0.002, colour = "black") +
  #geom_errorbarh(aes(xmin = QALYs_lower, xmax = QALYs_upper), height = 0.2, colour = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey") +
  scale_fill_manual(  # Use 'scale_fill_manual' instead of 'scale_color_manual'
    values = lancet_colors,
    labels = c("S0: Do-nothing", expression("S1: T+D"[agar1]), expression("S2: T+D"[agar2]), expression("S3: T+D"[pcr]), expression("S4: T+I"[agar1]), expression("S5: T+I"[agar2]), expression("S6: T+I"[pcr]), expression("S7: Pre-emptive I"[all]), expression("S8: Pre-emptive I"[men]), expression("S9: Pre-emptive I"[women])),
    name = "Strategy") +
  scale_y_continuous(breaks = seq(0, 56, by = 2))+
  #scale_x_continuous(breaks = seq(0, 0.012, by = 0.001))+
  labs(title = "",
       x = "Health benefit per new admission (∆QALYs)",
       y = "Incremental cost per new admission (∆$)",
       subtitle = "(A) CRE model") +
  theme_minimal() +
  theme_lancet() +
  theme(legend.position = "right",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10, hjust = 0),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
ggsave("newAdmgains_cre.tiff", newAdmgains_cre, width = 11, height = 7, dpi = 800)

# Adjust the individual plots before combining
newAdmgains_cre_updated <- newAdmgains_cre +
  theme(legend.position = "none", # Remove legend from the first plot
        text = element_text(size = 14), # Increase text size globally for the first plot
        axis.text = element_text(size = 11), # Increase axis text size for the first plot
        axis.title = element_text(size = 16)) # Increase axis title size for the first plot
newAdmgains_mrsa_updated <- newAdmgains_mrsa +
  theme(text = element_text(size = 14),
        axis.title.y = element_blank(), # Increase text size globally for the second plot
        axis.text = element_text(size = 11), # Increase axis text size for the second plot
        axis.title = element_text(size = 16))+
  scale_x_continuous(breaks = seq(0, 0.03, by = 0.005))
# Use patchwork to combine the plots
combined_plots <- (newAdmgains_cre_updated) | # Add letter (a) and ensure it's bold
  (newAdmgains_mrsa_updated)  
  plot_layout(guides = "collect") # Collect guides into a single legend
newAdmgains_combo<-combined_plots + theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) # Reduce margins around the combined plot
ggsave("newAdmgains_combo.tiff", newAdmgains_combo, width = 13, height = 7, dpi = 800)


##########################################################################################
##########################################################################################
# Figure SM: Dynamics of deaths and infections MRSA/CRE over time.
##########################################################################################
##########################################################################################

#MRSA model, MRSA-deaths
results_matrix <- as.data.frame(cbind(id, results_epi_MRSAdead))
results_matrix[1:200,3:17]<- NA
results_matrix[199:200,3:17]<- results_matrix[199:200,2]
colnames(results_matrix) <- c("id", "S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_epi_MRSAdead_l <- melt(results_matrix, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))
#FILTER INFORMATION TO ONLY 9 interventions:####
results_epi_MRSAdead_lf <- results_epi_MRSAdead_l[!grepl("S7|S8|S9|S10|S11|S12", results_epi_MRSAdead_l$variable), ]
results_epi_MRSAdead_lf$variable<- gsub("S13", "S7", results_epi_MRSAdead_lf$variable)
results_epi_MRSAdead_lf$variable <- gsub("S14", "S8", results_epi_MRSAdead_lf$variable)
results_epi_MRSAdead_lf$variable <- gsub("S15", "S9", results_epi_MRSAdead_lf$variable)
MRSA_dead_strategies<-ggplot(results_epi_MRSAdead_lf, aes(x = id, y = value, colour = variable)) + 
  geom_line(linewidth = 2.3) +
  scale_color_manual(values = lancet_colors, 
                     labels = c("S0: Do-nothing", expression(S1: T+D["agar1"]), expression(S2: T+D["agar2"]), expression(S3: T+D["pcr"]), expression(S4: T+I["agar1"]), expression(S5: T+I["agar2"]), expression(S6: T+I["pcr"]), expression("S7: T+D, men"["agar1"]), expression("S8: T+D, men"["agar2"]), expression("S9: T+D, men"["pcr"]), expression("S10: T+D, women"["agar1"]), expression("S11: T+D, women"["agar2"]), expression("S12: T+D, women"["pcr"]), expression("S13: Pre-emptive I"["all"]), expression("S14: Pre-emptive I"["men"]), expression("S15: Pre-emptive I"["women"]))) +
  labs(x = "Time (days)", 
       y = "Number of MRSA-associated deaths", 
       title = "(B) MRSA model",
       subtitle = "",
       colour = "Strategy   ") +
  scale_x_continuous(breaks = seq(170, 350, by = 20), limits = c(170, 350)) +
  scale_y_continuous(breaks = seq(60, 190, by = 10), limits = c(60, 190)) +
  geom_hline(yintercept = 60, linetype = "solid", color = "grey", size = 0.5) +  # Add horizontal line at y=75
  geom_vline(xintercept = 170, linetype = "solid", color = "grey", size = 0.5) +
  theme_lancet() +
  theme(legend.position = "right",
        legend.justification = "left",
        legend.title = element_text(size = 16),  # Increase size of legend title
        legend.text = element_text(hjust = 0, size = 14),  # Increase size of legend text
        axis.title.x = element_text(size = 21),  # Increase size of x-axis title
        axis.title.y = element_text(size = 21),  # Increase size of y-axis title
        axis.text.x = element_text(size = 19),  # Increase size of x-axis text/labels
        axis.text.y = element_text(size = 19),  # Increase size of y-axis text/labels
        plot.title = element_text(size = 16),   # Increase size of main title
        plot.subtitle = element_text(size = 14))  # Increase size of subtitle
ggsave("MRSA_dead_strategies.tiff", MRSA_dead_strategies, width = 11, height = 7, dpi = 800)
#  annotate("text", x = 150, y = 115, label = "(A)", size = 5, fontface = "bold", colour = "black",
#           hjust = 1.8, vjust = -1.2)

#MRSA model, infections all
results_matrix <- as.data.frame(cbind(id, results_epi_MRSAinfe_all))
results_matrix[1:200,3:17]<- NA
results_matrix[199:200,3:17]<- results_matrix[199:200,2]
colnames(results_matrix) <- c("id","S0", "S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_epi_MRSAinfe_all_l <- melt(results_matrix, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))

results_epi_MRSAinfe_all_lf <- results_epi_MRSAinfe_all_l[!grepl("S7|S8|S9|S10|S11|S12", results_epi_MRSAinfe_all_l$variable), ]
results_epi_MRSAinfe_all_lf$variable<- gsub("S13", "S7", results_epi_MRSAinfe_all_lf$variable)
results_epi_MRSAinfe_all_lf$variable <- gsub("S14", "S8", results_epi_MRSAinfe_all_lf$variable)
results_epi_MRSAinfe_all_lf$variable <- gsub("S15", "S9", results_epi_MRSAinfe_all_lf$variable)


staph_Infectionsall_strategies<-ggplot(results_epi_MRSAinfe_all_lf, aes(x = id, y = value, colour = variable)) + 
  geom_line(linewidth = 2.3) +
  scale_color_manual(values = lancet_colors, 
                     labels = c("S0: Do-nothing", expression(S1: T+D["agar1"]), expression(S2: T+D["agar2"]), expression(S3: T+D["pcr"]), expression(S4: T+I["agar1"]), expression(S5: T+I["agar2"]), expression(S6: T+I["pcr"]),  expression("S7: Pre-emptive I"["all"]), expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"]))) +
  labs(x = "Time (days)", 
       y = "Number of total infections", 
       title = "(B) MRSA model",
       subtitle = "",
       colour = "Strategy   ") +
  scale_x_continuous(breaks = seq(170, 350, by = 20), limits = c(170, 350)) +
  scale_y_continuous(breaks = seq(30, 39, by = 1), limits = c(30, 39)) +
  geom_hline(yintercept = 30, linetype = "solid", color = "grey", size = 0.5) +  # Add horizontal line at y=75
  geom_vline(xintercept = 170, linetype = "solid", color = "grey", size = 0.5) +
  theme_lancet() +
  theme(legend.position = "right",
        legend.justification = "left",
        legend.title = element_text(size = 16),  # Increase size of legend title
        legend.text = element_text(hjust = 0, size = 14),  # Increase size of legend text
        axis.title.x = element_text(size = 21),  # Increase size of x-axis title
        axis.title.y = element_text(size = 21),  # Increase size of y-axis title
        axis.text.x = element_text(size = 19),  # Increase size of x-axis text/labels
        axis.text.y = element_text(size = 19),  # Increase size of y-axis text/labels
        plot.title = element_text(size = 16),   # Increase size of main title
        plot.subtitle = element_text(size = 14))  # Increase size of subtitle
ggsave("staph_Infectionsall_strategies.tiff", staph_Infectionsall_strategies, width = 11, height = 7, dpi = 800)

#MRSA model, all deaths
results_matrix <- as.data.frame(cbind(id, results_epi_MRSAdead_all))
results_matrix[1:200,3:17]<- NA
results_matrix[199:200,3:17]<- results_matrix[199:200,2]
colnames(results_matrix) <- c("id","S0", "S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_epi_MRSAdead_all_l <- melt(results_matrix, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))

results_epi_MRSAdead_all_lf <- results_epi_MRSAdead_all_l[!grepl("S7|S8|S9|S10|S11|S12", results_epi_MRSAinfe_all_l$variable), ]
results_epi_MRSAdead_all_lf$variable<- gsub("S13", "S7", results_epi_MRSAdead_all_lf$variable)
results_epi_MRSAdead_all_lf$variable <- gsub("S14", "S8", results_epi_MRSAdead_all_lf$variable)
results_epi_MRSAdead_all_lf$variable <- gsub("S15", "S9", results_epi_MRSAdead_all_lf$variable)

staph_Deathsall_strategies<-ggplot(results_epi_MRSAdead_all_lf, aes(x = id, y = value, colour = variable)) + 
  geom_line(linewidth = 2.3) +
  scale_color_manual(values = lancet_colors, 
                     labels = c("S0: Do-nothing", expression(S1: T+D["agar1"]), expression(S2: T+D["agar2"]), expression(S3: T+D["pcr"]), expression(S4: T+I["agar1"]), expression(S5: T+I["agar2"]), expression(S6: T+I["pcr"]), expression("S7: Pre-emptive I"["all"]), expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"]))) +
  labs(x = "Time (days)", 
       y = "Number of total deaths", 
       title = "(B) MRSA model",
       subtitle = "",
       colour = "Strategy   ") +
  scale_x_continuous(breaks = seq(170, 350, by = 20), limits = c(170, 350)) +
  scale_y_continuous(breaks = seq(160, 350, by = 20), limits = c(160, 350)) +
  geom_hline(yintercept = 160, linetype = "solid", color = "grey", size = 0.5) +  # Add horizontal line at y=75
  geom_vline(xintercept = 170, linetype = "solid", color = "grey", size = 0.5) +
  theme_lancet() +
  theme(legend.position = "right",
        legend.justification = "left",
        legend.title = element_text(size = 16),  # Increase size of legend title
        legend.text = element_text(hjust = 0, size = 14),  # Increase size of legend text
        axis.title.x = element_text(size = 21),  # Increase size of x-axis title
        axis.title.y = element_text(size = 21),  # Increase size of y-axis title
        axis.text.x = element_text(size = 19),  # Increase size of x-axis text/labels
        axis.text.y = element_text(size = 19),  # Increase size of y-axis text/labels
        plot.title = element_text(size = 16),   # Increase size of main title
        plot.subtitle = element_text(size = 14))  # Increase size of subtitle
ggsave("staph_Deathsall_strategies.tiff", staph_Deathsall_strategies, width = 11, height = 7, dpi = 800)

#
#CRE model, CRE-deaths
results_matrix <- as.data.frame(cbind(id, results_epi_CREdead))
results_matrix[1:200,3:17]<- NA
results_matrix[199:200,3:17]<- results_matrix[199:200,2]
colnames(results_matrix) <- c("id", "S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_epi_CREdead_l <- melt(results_matrix, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))
results_epi_CREdead_l$value[results_epi_CREdead_l$variable == "S6" & results_epi_CREdead_l$id > 200] <- results_epi_CREdead_l$value[results_epi_CREinfe_l$variable == "S6" & results_epi_CREdead_l$id > 200] + 1

results_epi_CREdead_lf <- results_epi_CREdead_l[!grepl("S7|S8|S9|S10|S11|S12", results_epi_CREdead_l$variable), ]
results_epi_CREdead_lf$variable <- gsub("S13", "S7", results_epi_CREdead_lf$variable)
results_epi_CREdead_lf$variable <- gsub("S14", "S8", results_epi_CREdead_lf$variable)
results_epi_CREdead_lf$variable <- gsub("S15", "S9", results_epi_CREdead_lf$variable)
CRE_dead_strategies<- ggplot(results_epi_CREdead_lf, aes(x = id, y = value, colour = variable)) +
  geom_line(linewidth = 2.3) +
  scale_color_manual(values = lancet_colors, 
                     labels = c("S0: Do-nothing", expression(S1: T+D["agar1"]), expression(S2: T+D["agar2"]), expression(S3: T+D["pcr"]), expression(S4: T+I["agar1"]), expression(S5: T+I["agar2"]), expression(S6: T+I["pcr"]), expression("S7: Pre-emptive I"["all"]), expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"]))) +
  labs(x = "Time (days)", 
       y = "Number of CRE-associated deaths", 
       title = "(A) CRE model",
       subtitle = "",
       colour = "Strategy") +
  scale_x_continuous(breaks = seq(170, 365, by = 20), limits = c(170, 365)) +
  scale_y_continuous(breaks = seq(80, 210, by = 10), limits = c(80, 210)) +
  geom_hline(yintercept = 80, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 170, linetype = "solid", color = "grey", size = 0.5) +
  theme_lancet()+
  theme(legend.position = "right",
        legend.justification = "left",
        legend.title = element_text(size = 16),  # Increase size of legend title
        legend.text = element_text(hjust = 0, size = 14),  # Increase size of legend text
        axis.title.x = element_text(size = 21),  # Increase size of x-axis title
        axis.title.y = element_text(size = 21),  # Increase size of y-axis title
        axis.text.x = element_text(size = 19),  # Increase size of x-axis text/labels
        axis.text.y = element_text(size = 19),  # Increase size of y-axis text/labels
        plot.title = element_text(size = 16),   # Increase size of main title
        plot.subtitle = element_text(size = 14))  # Increase size of subtitle
# Print the plot
#print(CRE_inf_strategies)
ggsave("CRE_dead_strategies.tiff",CRE_dead_strategies, width = 11, height = 7, dpi = 800)


#CRE model, infections all
results_matrix <- as.data.frame(cbind(id, results_epi_CREinfe_all))
results_matrix[1:200,3:17]<- NA
results_matrix[199:200,3:17]<- results_matrix[199:200,2]
colnames(results_matrix) <- c("id", "S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_epi_CREinfe_all_l <- melt(results_matrix, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))
results_epi_CREinfe_all_l$value[results_epi_CREinfe_all_l$variable == "S6" & results_epi_CREinfe_all_l$id > 200] <- results_epi_CREinfe_all_l$value[results_epi_CREinfe_all_l$variable == "S6" & results_epi_CREinfe_all_l$id > 200] + 0.1

results_epi_CREinfe_all_lf <- results_epi_CREinfe_all_l[!grepl("S7|S8|S9|S10|S11|S12", results_epi_CREinfe_all_l$variable), ]
results_epi_CREinfe_all_lf$variable <- gsub("S13", "S7", results_epi_CREinfe_all_lf$variable)
results_epi_CREinfe_all_lf$variable <- gsub("S14", "S8", results_epi_CREinfe_all_lf$variable)
results_epi_CREinfe_all_lf$variable <- gsub("S15", "S9", results_epi_CREinfe_all_lf$variable)
ent_Infectionsall_strategies<- ggplot(results_epi_CREinfe_all_lf, aes(x = id, y = value, colour = variable)) +
  geom_line(linewidth = 2.3) +
  scale_color_manual(values = lancet_colors, 
                     labels = c("S0: Do-nothing", expression(S1: T+D["agar1"]), expression(S2: T+D["agar2"]), expression(S3: T+D["pcr"]), expression(S4: T+I["agar1"]), expression(S5: T+I["agar2"]), expression(S6: T+I["pcr"]), expression("S7: Pre-emptive I"["all"]), expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"]))) +
  labs(x = "Time (days)", 
       y = "Number of total infections", 
       title = "(A) CRE model",
       subtitle = "",
       colour = "Strategy") +
  scale_x_continuous(breaks = seq(170, 365, by = 20), limits = c(170, 365)) +
  scale_y_continuous(breaks = seq(32, 40, by = 1), limits = c(32, 40)) +
  geom_hline(yintercept = 32, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 170, linetype = "solid", color = "grey", size = 0.5) +
  theme_lancet()+
  theme(legend.position = "right",
        legend.justification = "left",
        legend.title = element_text(size = 16),  # Increase size of legend title
        legend.text = element_text(hjust = 0, size = 14),  # Increase size of legend text
        axis.title.x = element_text(size = 21),  # Increase size of x-axis title
        axis.title.y = element_text(size = 21),  # Increase size of y-axis title
        axis.text.x = element_text(size = 19),  # Increase size of x-axis text/labels
        axis.text.y = element_text(size = 19),  # Increase size of y-axis text/labels
        plot.title = element_text(size = 16),   # Increase size of main title
        plot.subtitle = element_text(size = 14))  # Increase size of subtitle
# Print the plot
#print(CRE_inf_strategies)
ggsave("ent_Infectionsall_strategies.tiff",ent_Infectionsall_strategies, width = 11, height = 7, dpi = 800)

#CRE model, deaths all
results_matrix <- as.data.frame(cbind(id, results_epi_CREdead_all))
results_matrix[1:200,3:17]<- NA
results_matrix[199:200,3:17]<- results_matrix[199:200,2]
colnames(results_matrix) <- c("id","S0", "S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15")
results_epi_CREdead_all_l <- melt(results_matrix, id.vars = "id", measure.vars = c("S0","S1", "S2", "S3", "S4", "S5","S6", "S7", "S8", "S9", "S10","S11", "S12", "S13", "S14", "S15"))
results_epi_CREdead_all_l$value[results_epi_CREdead_all_l$variable == "S6" & results_epi_CREdead_all_l$id > 200] <- results_epi_CREdead_all_l$value[results_epi_CREdead_all_l$variable == "S6" & results_epi_CREdead_all_l$id > 200] + 3

results_epi_CREdead_all_lf <- results_epi_CREdead_all_l[!grepl("S7|S8|S9|S10|S11|S12", results_epi_CREdead_all_l$variable), ]
results_epi_CREdead_all_lf$variable <- gsub("S13", "S7", results_epi_CREdead_all_lf$variable)
results_epi_CREdead_all_lf$variable <- gsub("S14", "S8", results_epi_CREdead_all_lf$variable)
results_epi_CREdead_all_lf$variable <- gsub("S15", "S9", results_epi_CREdead_all_lf$variable)
ent_Deathsall_strategies<- ggplot(results_epi_CREdead_all_lf, aes(x = id, y = value, colour = variable)) +
  geom_line(linewidth = 2.3) +
  scale_color_manual(values = lancet_colors, 
                     labels = c("S0: Do-nothing", expression(S1: T+D["agar1"]), expression(S2: T+D["agar2"]), expression(S3: T+D["pcr"]), expression(S4: T+I["agar1"]), expression(S5: T+I["agar2"]), expression(S6: T+I["pcr"]), expression("S7: Pre-emptive I"["all"]), expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"]))) +
  labs(x = "Time (days)", 
       y = "Number of total deaths", 
       title = "(A) CRE model",
       subtitle = "",
       colour = "Strategy") +
  scale_x_continuous(breaks = seq(170, 365, by = 20), limits = c(170, 365)) +
  scale_y_continuous(breaks = seq(150, 340, by = 20), limits = c(150, 340)) +
  geom_hline(yintercept = 150, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 170, linetype = "solid", color = "grey", size = 0.5) +
  theme_lancet()+
  theme(legend.position = "right",
        legend.justification = "left",
        legend.title = element_text(size = 16),  # Increase size of legend title
        legend.text = element_text(hjust = 0, size = 14),  # Increase size of legend text
        axis.title.x = element_text(size = 21),  # Increase size of x-axis title
        axis.title.y = element_text(size = 21),  # Increase size of y-axis title
        axis.text.x = element_text(size = 19),  # Increase size of x-axis text/labels
        axis.text.y = element_text(size = 19),  # Increase size of y-axis text/labels
        plot.title = element_text(size = 16),   # Increase size of main title
        plot.subtitle = element_text(size = 14))  # Increase size of subtitle
# Print the plot
#print(CRE_inf_strategies)
ggsave("ent_Deathsall_strategies.tiff",ent_Deathsall_strategies, width = 11, height = 7, dpi = 800)


#----------------------------------------------------#
#Combine FIGURES.
#----------------------------------------------------#
#all infections
ent_Infectionsall_strategies<- ent_Infectionsall_strategies + theme(legend.position = "none")
combined_MRSACRE_inf_all <-  ent_Infectionsall_strategies / staph_Infectionsall_strategies +  plot_layout(guides = 'collect')  
theme_void()# Eq# Stack them vertically
ggsave("combined_MRSACRE_inf_all.tiff",combined_MRSACRE_inf_all, width = 12, height = 15, dpi = 800)
#MRSA- or CRE-associated deaths
combined_MRSACRE_death <- CRE_dead_strategies /MRSA_dead_strategies # Stack them vertically
ggsave("combined_MRSACRE_death.tiff",combined_MRSACRE_death, width = 12, height = 15, dpi = 800)
#MRSA/ CRE models, all deaths
combined_MRSACRE_death_all <- ent_Deathsall_strategies /staph_Deathsall_strategies # Stack them vertically
ggsave("combined_MRSACRE_death_all.tiff",combined_MRSACRE_death_all, width = 12, height = 15, dpi = 800)


##############################################################################
##############################################################################
# Figure SM: Heatmap for sensitivity of the test & time delay
##############################################################################
##############################################################################

#-----------------------------------------------------#
#First matrix for T+Decolonisation strategies:
#-----------------------------------------------------#
# Define the range for the parameters
sens_chrom_a_values <- seq(1, 20, by = 1)
turn_chrom_a_values <- seq(1, 6, by = 1)
# Prepare a matrix to store results
results_matrixCRE <- matrix(nrow = length(sens_chrom_a_values), ncol = length(turn_chrom_a_values))
results_matrixMRSA <- matrix(nrow = length(sens_chrom_a_values), ncol = length(turn_chrom_a_values))
# Loop through the values
for (i in 1:length(sens_chrom_a_values)) {
  for (j in 1:length(turn_chrom_a_values)) {
    # Update the parameters
    parameters2["sens_chrom_a"] <- sens_chrom_a_values[i]/20
    parameters2["turn_chrom_a"] <- turn_chrom_a_values[j]
    parameters1["sens_chrom_a"] <- sens_chrom_a_values[i]/20
    parameters1["turn_chrom_a"] <- turn_chrom_a_values[j]
    # Run the ODE
    outCRE <- ode(y = state2, times = times, func = ARB_model_2ch_td_newadm, parms = parameters2, method = "rk4")
    outMRSA <- ode(y = state1, times = times, func = ARB_model_1ch_td_newadm, parms = parameters1, method = "rk4")
    # Extract the outcome of interest, e.g., the last value or the max value of 'infected' individuals
    # Assuming 'infected' is one of the state variables, change it to whatever your model uses
    # Here, I'm assuming the outcome is the last point for simplicity; adjust as needed
    final_valueCRE <-  sum(outCRE[1:366, "IMR_m2"]) + sum(outCRE[1:366, "IMR_f2"])+ sum(outCRE[1:366, "ISR_m2"])+sum(outCRE[1:366, "ISR_f2"]) # Change 'infected' to your specific state variable name
    final_valueMRSA <- sum(outMRSA[1:366, "IMR_m1"]) + sum(outMRSA[1:366, "IMR_f1"])+ sum(outMRSA[1:366, "ISR_m1"])+sum(outMRSA[1:366, "ISR_f1"])
    results_matrixCRE[i, j] <- final_valueCRE
    results_matrixMRSA[i, j] <- final_valueMRSA
  }
}
#ARB_model_2ch_tiso_newadm ARB_model_2ch_td_newadm
# Convert the matrix to a format suitable for ggplot
results_meltedCRE <- melt(results_matrixCRE, varnames = c("sens_chrom_a", "turn_chrom_a"), value.name = "Infections")
results_meltedMRSA <- melt(results_matrixMRSA, varnames = c("sens_chrom_a", "turn_chrom_a"), value.name = "Infections")
# Generate the plot
results_meltedCRE<- as.data.frame(results_meltedCRE)
results_meltedMRSA<- as.data.frame(results_meltedMRSA)
results_meltedMRSA$sens_chrom_a<- (results_meltedMRSA$sens_chrom_a/20)*100
results_meltedCRE$sens_chrom_a<- (results_meltedCRE$sens_chrom_a/20)*100

cre_matrixSTdelay<- ggplot(results_meltedCRE, aes(x = sens_chrom_a, y = turn_chrom_a, fill = value)) +
  geom_tile(color = "black") +  # Add black borders to each tile
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +  # Ensuring breaks align with rounded values
  labs(x = "Sensitivity of the test (%)", 
       y = "Time delay until results (days)", 
       title = "(A) CRE",
       fill = "Number of\nInfections") +  # Adjusted for line break
  theme_lancet() +
  scale_y_continuous(breaks = 1:6) +
  scale_x_continuous(breaks = seq(5, 100, by = 5)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 0.5), 
        axis.text.y = element_text(angle = 0))
ggsave("cre_matrixSTdelay1.tiff",cre_matrixSTdelay, width = 11, height = 8, dpi = 800)

mrsa_matrixSTdelay <- ggplot(results_meltedMRSA, aes(x = sens_chrom_a, y = turn_chrom_a, fill = value)) +
  geom_tile(color = "black") +  # Add black borders to each tile
  scale_fill_gradientn(colors = brewer.pal(9, "Greens")) +  # Ensuring breaks align with rounded values
  labs(x = "Sensitivity of the test (%)", 
       y = "Time delay until results (days)",
       title = "(B) MRSA",
       fill = "Number of\nInfections") +  # Adjusted for line break
  theme_lancet() +
  scale_y_continuous(breaks = 1:6) +
  scale_x_continuous(breaks = seq(5, 100, by = 5)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 0.5), 
        axis.text.y = element_text(angle = 0))
ggsave("mrsa_matrixSTdelay1.tiff",mrsa_matrixSTdelay, width = 11, height = 8, dpi = 800)

mrsa_matrixSTdelay <- mrsa_matrixSTdelay  +# Stack them vertically
  theme(
    text = element_text(size = 17),  # Default text size for all text in the plot
    axis.title = element_text(size = 19),  # Increase size of axis titles
    axis.text = element_text(size = 17),   # Increase size of axis tick labels
    legend.text = element_text(size = 18), # Increase size of legend text
    legend.title = element_text(size = 19),
    title= element_text(size=22)# Increa# Increase size of legend title
  )

cre_matrixSTdelay <- cre_matrixSTdelay  +# Stack them vertically
  theme(
    text = element_text(size = 17),  # Default text size for all text in the plot
    axis.title = element_text(size = 19),  # Increase size of axis titles
    axis.text = element_text(size = 17),   # Increase size of axis tick labels
    legend.text = element_text(size = 18), # Increase size of legend text
    legend.title = element_text(size = 19),
    title= element_text(size=22)# Increase size of legend title
  )
combined_MRSACRE_TEST_TD1 <- cre_matrixSTdelay / mrsa_matrixSTdelay 
ggsave("combined_MRSACRE_TEST_TD1.tiff",combined_MRSACRE_TEST_TD1, width = 12, height = 15, dpi = 800)

write.csv(results_meltedCRE, "results_meltedCRE.csv", row.names = TRUE)
write.csv(results_meltedMRSA, "results_meltedMRSA.csv", row.names = TRUE)

#-----------------------------------------------------#
#Second matrix for test+isolation, sensit & time delay
#-----------------------------------------------------#

sens_chrom_a_values <- seq(1, 20, by = 1)
turn_chrom_a_values <- seq(1, 6, by = 1)
# Prepare a matrix to store results
results_matrixCRE <- matrix(nrow = length(sens_chrom_a_values), ncol = length(turn_chrom_a_values))
results_matrixMRSA <- matrix(nrow = length(sens_chrom_a_values), ncol = length(turn_chrom_a_values))
# Loop through the values
for (i in 1:length(sens_chrom_a_values)) {
  for (j in 1:length(turn_chrom_a_values)) {
    # Update the parameters
    parameters2["sens_chrom_a"] <- sens_chrom_a_values[i]/20
    parameters2["turn_chrom_a"] <- turn_chrom_a_values[j]
    parameters1["sens_chrom_a"] <- sens_chrom_a_values[i]/20
    parameters1["turn_chrom_a"] <- turn_chrom_a_values[j]
    # Run the ODE
    outCRE <- ode(y = state2, times = times, func = ARB_model_2ch_tiso_newadm, parms = parameters2, method = "rk4")
    outMRSA <- ode(y = state1, times = times, func = ARB_model_1ch_tiso_newadm, parms = parameters1, method = "rk4")
    # Extract the outcome of interest, e.g., the last value or the max value of 'infected' individuals
    # Assuming 'infected' is one of the state variables, change it to whatever your model uses
    # Here, I'm assuming the outcome is the last point for simplicity; adjust as needed
    final_valueCRE <-  sum(outCRE[1:366, "IMR_m2"]) + sum(outCRE[1:366, "IMR_f2"])+ sum(outCRE[1:366, "ISR_m2"])+sum(outCRE[1:366, "ISR_f2"]) # Change 'infected' to your specific state variable name
    final_valueMRSA <- sum(outMRSA[1:366, "IMR_m1"]) + sum(outMRSA[1:366, "IMR_f1"])+ sum(outMRSA[1:366, "ISR_m1"])+sum(outMRSA[1:366, "ISR_f1"])
    results_matrixCRE[i, j] <- final_valueCRE
    results_matrixMRSA[i, j] <- final_valueMRSA
  }
}
#ARB_model_2ch_tiso_newadm ARB_model_2ch_td_newadm
# Convert the matrix to a format suitable for ggplot
results_meltedCRE <- melt(results_matrixCRE, varnames = c("sens_chrom_a", "turn_chrom_a"), value.name = "Infections")
results_meltedMRSA <- melt(results_matrixMRSA, varnames = c("sens_chrom_a", "turn_chrom_a"), value.name = "Infections")
# Generate the plot
results_meltedCRE<- as.data.frame(results_meltedCRE)
results_meltedMRSA<- as.data.frame(results_meltedMRSA)
results_meltedMRSA$sens_chrom_a<- (results_meltedMRSA$sens_chrom_a/20)*100
results_meltedCRE$sens_chrom_a<- (results_meltedCRE$sens_chrom_a/20)*100
library(scales)
library(RColorBrewer)

cre_matrixSTdelay<- ggplot(results_meltedCRE, aes(x = sens_chrom_a, y = turn_chrom_a, fill = value)) +
  geom_tile(color = "black") +  # Add black borders to each tile
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +  # Ensuring breaks align with rounded values
  labs(x = "Sensitivity of the test (%)", 
       y = "Time delay until results (days)", 
       title = "(A) CRE",
       fill = "Number of\nInfections") +  # Adjusted for line break
  theme_lancet() +
  scale_y_continuous(breaks = 1:6) +
  scale_x_continuous(breaks = seq(5, 100, by = 5)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 0.5), 
        axis.text.y = element_text(angle = 0))
ggsave("cre_matrixSTdelay2.tiff",cre_matrixSTdelay, width = 11, height = 8, dpi = 800)


mrsa_matrixSTdelay <- ggplot(results_meltedMRSA, aes(x = sens_chrom_a, y = turn_chrom_a, fill = value)) +
  geom_tile(color = "black") +  # Add black borders to each tile
  scale_fill_gradientn(colors = brewer.pal(9, "Greens")) +  # Ensuring breaks align with rounded values
  labs(x = "Sensitivity of the test (%)", 
       y = "Time delay until results (days)",
       title = "(B) MRSA",
       fill = "Number of\nInfections") +  # Adjusted for line break
  theme_lancet() +
  scale_y_continuous(breaks = 1:6) +
  scale_x_continuous(breaks = seq(5, 100, by = 5)) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 0.5), 
        axis.text.y = element_text(angle = 0))
ggsave("mrsa_matrixSTdelay2.tiff",mrsa_matrixSTdelay, width = 11, height = 8, dpi = 800)
library(grid)  # for textGrob
library(ggplot2)
library(gridExtra) 
mrsa_matrixSTdelay <- mrsa_matrixSTdelay  +# Stack them vertically
  theme(
    text = element_text(size = 17),  # Default text size for all text in the plot
    axis.title = element_text(size = 19),  # Increase size of axis titles
    axis.text = element_text(size = 17),   # Increase size of axis tick labels
    legend.text = element_text(size = 18), # Increase size of legend text
    legend.title = element_text(size = 19),
    title= element_text(size=22)# Increa# Increase size of legend title
  )

cre_matrixSTdelay <- cre_matrixSTdelay  +# Stack them vertically
  theme(
    text = element_text(size = 17),  # Default text size for all text in the plot
    axis.title = element_text(size = 19),  # Increase size of axis titles
    axis.text = element_text(size = 17),   # Increase size of axis tick labels
    legend.text = element_text(size = 18), # Increase size of legend text
    legend.title = element_text(size = 19),
    title= element_text(size=22)# Increase size of legend title
  )
combined_MRSACRE_TEST_TD2 <- cre_matrixSTdelay / mrsa_matrixSTdelay 
ggsave("combined_MRSACRE_TEST_TD2.tiff",combined_MRSACRE_TEST_TD2, width = 12, height = 15, dpi = 800)

write.csv(results_meltedCRE, "results_meltedCRE2.csv", row.names = TRUE)
write.csv(results_meltedMRSA, "results_meltedMRSA2.csv", row.names = TRUE)


################################################################################################
################################################################################################
# Figure SM: Number of colonised/infected CRE/MRSA per 100 new admissions with 95%CI adjusted to +-1SD prevalence
################################################################################################
################################################################################################
#Upper bound: 0.1774853 R.   0.1234747 S ///// 0.2229124 transmission parameter.
# Define the models you want to run according to the functions (per strategy)
#UPPERBOUNDestimatesMRSA#####
model_list <- c("ARB_model_1ch_do_nothing", "ARB_model_1ch_td_newadm","ARB_model_1ch2_td_newadm","ARB_model_1pcr_td_newadm","ARB_model_1ch_tiso_newadm","ARB_model_1ch2_tiso_newadm","ARB_model_1pcr_tiso_newadm",
                "ARB_model_1ch_td_newadmHR_m","ARB_model_1ch2_td_newadmHR_m","ARB_model_1_pcr_td_newadmHR_m","ARB_model_1ch_td_newadmHR_f","ARB_model_1ch2_td_newadmHR_f","ARB_model_1_pcr_td_newadmHR_f","ARB_model_1preE_newadm_all","ARB_model_1preE_newadm_m","ARB_model_1preE_newadm_f")  # Add last interventions with high-risk groups and pre-emptive isolation
results_epi_ac <- matrix(nrow = 10, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
results_epi_MRSAprev <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAinfe <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAdead <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAinfe_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAdead_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_uncolonised <- matrix(nrow = 366, ncol = length(model_list))

results_econ <- matrix(nrow = 3, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
for (i in 1:length(model_list)) {
  #Run initial conditions first:
  times <- seq(from=0, to=365, by = 1)  # Simulate over a year
  #Staphylococcus aureus params [MRSA/MSSA]####### 
  # ----------------------------------------#
  zeta3_m_p1 <- 0.231  # mortality rate from IMS  , male
  zeta3_f_p1 = zeta3_m_p1*2.07  # mortality rate from IMS  , female
  zeta1_m_p1 = zeta3_m_p1*1.01  # mortality rate from IMR, male
  zeta1_f_p1 = zeta3_m_p1*1.22  # mortality rate from IMR, female
  zeta2_m_p1 = zeta3_m_p1*2.32  # mortality rate from ISR, male
  zeta2_f_p1 = zeta3_m_p1*2.25  # mortality rate from ISR, female
  zeta4_m_p1 = zeta3_m_p1*1.10  # mortality rate from ISS , male
  zeta4_f_p1 = zeta3_m_p1*2.34  # mortality rate from ISS , female
  #(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  kappa_p1 = (35.9-6.2)/(49.5-4.85)
  parameters1 <- c(
    delta1_p1 = 0.0016,  
    delta2_p1 = 0.0016,  
    Disch_U_f_p1 = 1/6,
    Disch_U_m_p1 = 1/6,
    Disch_CR_f_p1 = 1/6,
    Disch_CR_m_p1 = 1/6,
    Disch_CS_f_p1 = 1/6,
    Disch_CS_m_p1 = 1/6,
    mu0_p1 = 0.52,  # % of women among U
    mu1_p1 = 0.52,  # % of women among CR
    mu2_p1 = 0.52,  # % of women among CS
    mu3_p1 = 0.3165,  # % of women among IMR
    mu4_p1 = 0.3889,  # % of women among ISR
    mu5_p1 = 0.3938,  # % of women among IMS
    mu6_p1 = 0.4487,  # % of women among IMS
    #Exposure to anbiotics; treatment amongst susceptible populations
    psi_m_p1 = 0.1474, # % of  individuals exposed  to vancomycin/penicillin among males
    psi_w_p1 = 0.184,  # % of  individuals exposed to vancomycin/penicillin among females
    #Percentage of people under treatment for MRSA decolonisation
    psi_mtr_p1=0.109, #IEAT is 1.35 times higher among CR; hence if treatment is psi_m_p1; psi_mtr_p1=psi_m_p1/1.35
    psi_wtr_p1=0.136,
    #Fitness cost. c reduces the transmission rate among resistant strains [1/unit time] [%].
    c_p1=0.09,
    #Progression to the development of infection from colonisation among CR and CS states. 
    beta1_m_p1 = (1/21)*0.26,  # inverse of LOS plus progression from colonisation to infection among males CR
    beta2_m_p1 = (1/11)*0.099,  # inverse of LOS plus progression from colonisation to infection among males CS
    beta1_f_p1 = (1/29)*0.26, # inverse of LOS plus progression from colonisation to infection among females CR
    beta2_f_p1 = (1/14)*0.099,  # inverse of LOS plus progression from colonisation to infection among females CS
    #Natural clearance of mild and severe infections among CR and CS states, respectively [1/unit time] [%].
    gamma1_p1 = 0.001,  # Natural clearance among mild infections R
    gamma2_p1 = 0.001,  # Natural clearance among severe infections R
    gamma3_p1 = 0.001,  # Natural clearance among mild infections S
    gamma4_p1 = 0.001, # Natural clearance among severe infections S
    #Mean time of infection considering length of hospital stays [1/length of hospital stay]. 
    omega1_d_m_p1 = (1/20), #IMR patients who died, male
    omega1_r_m_p1 = (1/23), #IMR patients who recovered, male
    omega1_d_f_p1 = (1/13), #IMR patients who died, female
    omega1_r_f_p1 = (1/26), #IMR patients who recovered, female
    omega2_d_m_p1 = (1/11),  #ISR patients who died, male
    omega2_r_m_p1 = (1/18),  #ISR patients who recovered, male
    omega2_d_f_p1 = (1/14), #ISR patients who died, female 
    omega2_r_f_p1 = (1/19), #ISR patients who recovered, female
    omega3_d_m_p1 = (1/12), #IMS patients who died, male
    omega3_r_m_p1 = (1/11), #IMS patients who recovered, male
    omega3_d_f_p1 = (1/20), #IMS patients who died, female
    omega3_r_f_p1 = (1/16), #IMS patients who recovered, female
    omega4_d_m_p1 = (1/11), #ISS patients who died, male
    omega4_r_m_p1 = (1/19), #ISS patients who recovered, male
    omega4_d_f_p1 = (1/14), #ISS patients who died, female
    omega4_r_f_p1 = (1/17), #ISS patients who recovered, female
    #Percentage of inpatients with CR or CS, respectively, progressing to severe infection in intensive care units [1/unit time] [%].
    alpha1_f_p1 = 0.4196, #% patients with CR progressing to severe infection, males
    alpha2_f_p1 = 0.3438, #% patients with CS progressing to severe infection, males
    alpha1_m_p1 = 0.3517, #% patients with CR progressing to severe infection, females
    alpha2_m_p1 = 0.2826, #% patients with CS progressing to severe infection, females
    #Progression from mild to severe infection from IMR and IMS, respectively [1/unit time] [%].
    epsilon1_p1 = 0.01, #progression from IMR to ISR
    epsilon2_p1 = 0.01, #progression from IMS to ISS
    #Mortality rates from infection. ζ1 and ζ2 are mortality rates from mild and severe resistant infections, respectively. ζ3 and ζ4 are from mild and severe susceptible infections, respectively [1/unit time] [%].
    zeta3_m_p1 = 0.231,  # mortality rate from IMS  , male
    zeta3_f_p1 = 0.231*2.07,  # mortality rate from IMS  , female
    zeta1_m_p1 = 0.231*1.01,  # mortality rate from IMR, male
    zeta1_f_p1 = 0.231*1.22,  # mortality rate from IMR, female
    zeta2_m_p1 = 0.231*2.32,  # mortality rate from ISR, male
    zeta2_f_p1 = 0.231*2.25,  # mortality rate from ISR, female
    zeta4_m_p1 = 0.231*1.10,  # mortality rate from ISS , male
    zeta4_f_p1 = 0.231*2.34,  # mortality rate from ISS , female
    #Recovery rates from infection, including IMR, ISR, IMS and ISS due to treatment received [1/unit time] [%].
    nu1_m_p1 = (1-zeta1_m_p1), #recovery rates from IMR, males
    nu1_f_p1 = (1-zeta1_f_p1), #recovery rates from IMR, females
    nu2_m_p1 = (1-zeta2_m_p1), #recovery rates from ISR, males
    nu2_f_p1 = (1-zeta2_f_p1), #recovery rates from ISR, females
    nu3_m_p1 = (1-zeta3_m_p1), #recovery rates from IMS, males
    nu3_f_p1 = (1-zeta3_f_p1), #recovery rates from IMS, females
    nu4_m_p1 =  (1-zeta4_m_p1), #recovery rates from ISS, males
    nu4_f_p1 =  (1-zeta4_f_p1),#recovery rates from ISS, females
    #Constant background rate that captures transmission from non-human sources, horizontal transmission, or de novo emergence [1/unit time] [number].
    b_p1 = 0.01,
    # Percentage of people with resistant infections receiving inappropriate empirical antibiotic treatment [1/unit time] [%].
    phi_m_p1 = 0.459,  # 
    phi_f_p1 = 0.413,  # 
    #Factor of burden associated with inappropriate empirical antibiotic treatment and increased ICU admission among resistant infections [1/unit time] [%].
    pi_p1= 1.35,
    #Transmission parameter {update this correspondingly after calibrating it with real data}
    #tau_p1= 0.03461113, #estimated 
    tau_p1=0.2829124,
    #community-acquired infection upon hospital admission rate
    caIha_p1= 0.001,
    ###
    
    #percentage of people tested
    test_p1=0.20, 
    HR_perc1=0.2,
    or_HR_scenar1_a=1.04,
    or_HR_scenarMen_a=1.37,
    ### ### ### ### ###
    #sensitivity chrom_1
    sens_chrom_a=0.826,
    #sensitivity chrom_1
    sens_chrom2_a=0.622, 
    #sensitivity chrom_1
    sens_pcr_a=0.881,
    #turnaround chrom_1
    turn_chrom_a=3,
    #turnaround chrom_1
    turn_chrom2_a=2,
    #turnaround pcr_1
    turn_pcr_a=1,
    #isolation contact precaution transmission reduction
    reduc_conpre_a=0.365,
    #efficiency decolonisation
    eff_decol_a= 0.53,
    #effect on self-infection decolonisation
    eff_decol_selfi_a=0.32,
    #Turnaround decolonisation program in days
    turnaround_decol_a=5,
    ##
    #costs wards
    c_general_ward= 50,
    c_intermediate_ward=92,
    c_icu_ward=218,
    c_decol_1pd=6.5,
    c_isolation=42.3,
    c_chrom=10.2,
    c_chrom2=13.6, 
    c_pcr=33,
    c_bc=16.9,
    #utilities
    u_healthy=0.92,
    u_icu=0.58,
    u_gw=0.64,
    u_recovICU=0.74
  )
  
  
  
  
  # ----------------------------------------#
  #Staphylococcus aureus states [MRSA/MSSA]####### 
  # ----------------------------------------#
  N <- 1000  # Total population size
  # Initial conditions (population sizes in each group)
  U_m10 <- 0.7 * N *(1-0.52)
  CR_m10 <- 0.1775 * N *(1-0.52)
  CS_m10 <-  0.1235* N *(1-0.52)
  IMR_m10 <- 0.26 * CR_m10* (1-0.4196)
  ISR_m10 <- 0.26 * CR_m10* 0.4196
  IMS_m10 <- 0.099 * CS_m10*(1-0.3438)
  ISS_m10 <- 0.099 * CS_m10*0.3438
  RR_m10 <-0
  RS_m10  <-0
  DR_m10  <-0
  DS_m10 <-0
  N_to0 <- 1050
  utility_to0<-0
  cost_to0<-0
  new_admin0<-0
  discharge<-0
  
  U_f10 <- 0.7 * N *0.52
  CR_f10 <- 0.1775 * N *0.52
  CS_f10 <-  0.1235* N *0.52
  IMR_f10 <- 0.26 * CR_f10* (1-0.3517)
  ISR_f10 <- 0.26 * CR_f10*0.3517
  IMS_f10 <- 0.099 * CS_f10*(1-0.2826)
  ISS_f10 <- 0.099 * CS_f10*0.2826
  RR_f10 <-0
  RS_f10  <-0
  DR_f10  <-0
  DS_f10 <-0
  
  N_0m10<-  U_m10 + CR_m10 + CS_m10 + IMR_m10 + ISR_m10 + IMS_m10 + ISS_m10 +  RR_m10 + RS_m10 + DR_m10 +DS_m10
  N_0f10<-  U_f10 + CR_f10 + CS_f10 + IMR_f10 + ISR_f10 + IMS_f10 + ISS_f10 +  RR_f10 + RS_f10 + DR_f10 +DS_f10
  
  state1 <- c(U_m1 = U_m10, CR_m1=CR_m10, CS_m1= CS_m10, IMR_m1= IMR_m10, ISR_m1=ISR_m10, IMS_m1= IMS_m10, ISS_m1= ISS_m10, RR_m1= RR_m10, RS_m1=RS_m10, DR_m1= DR_m10, DS_m1=DS_m10,
              U_f1 = U_f10, CR_f1=CR_f10, CS_f1= CS_f10, IMR_f1= IMR_f10, ISR_f1=ISR_f10, IMS_f1= IMS_f10, ISS_f1= ISS_f10, RR_f1= RR_f10, RS_f1=RS_f10, DR_f1= DR_f10, DS_f1=DS_f10, N_to=N_to0, utility=utility_to0, cost=cost_to0, new_admin=new_admin0, discharge=discharge0)
  
  N_orig1<-N_0m10 + N_0f10
  N_tdif <- N_orig1
  #####
  # Construct the function name from the model_list string
  model_func_name <- get(model_list[i])
  # Solve the ODE using the 'deSolve' package's ode function
  O_solution <- ode(y = state1, times = times, func = model_func_name, parms = parameters1, method = "rk4")
  O_solution2 <- as.data.frame(O_solution)
  O_solution <- as.data.frame(O_solution)
  # Process the ODE solution: sum all values reported in the first 350 rows for your variable of interest
  results_epi_ac[1, i] <- sum(O_solution[1:366, "CR_m1"]) + sum(O_solution[1:366, "CR_f1"]) + sum(O_solution[1:366, "IMR_m1"]) + sum(O_solution[1:366, "IMR_f1"])+ sum(O_solution[1:366, "ISR_m1"])+sum(O_solution[1:366, "ISR_f1"])
  results_epi_ac[2, i] <- sum(O_solution[1:366, "IMR_m1"]) + sum(O_solution[1:366, "IMR_f1"])+ sum(O_solution[1:366, "ISR_m1"])+sum(O_solution[1:366, "ISR_f1"])
  results_epi_ac[3, i] <- (O_solution[366, "DR_f1"]) + (O_solution[366, "DR_m1"]) 
  results_epi_ac[4, i] <- sum(O_solution[1:366, "IMR_m1"]) + sum(O_solution[1:366, "IMR_f1"])+ sum(O_solution[1:366, "ISR_m1"])+sum(O_solution[1:366, "ISR_f1"]) + sum(O_solution[1:366, "IMS_m1"]) + sum(O_solution[1:366, "IMS_f1"])+ sum(O_solution[1:366, "ISS_m1"])+sum(O_solution[1:366, "ISS_f1"])
  results_epi_ac[5, i] <- (O_solution[366, "DR_f1"]) + (O_solution[366, "DR_m1"]) + (O_solution[366, "DS_f1"]) + (O_solution[366, "DS_m1"]) 
  results_epi_ac[6, i] <- sum(O_solution[366, "new_admin"]) 
  results_epi_ac[7, i] <- sum(O_solution[1:366, "U_m1"]) +sum(O_solution[1:366, "U_f1"]) 
  results_epi_ac[8, i] <- sum(O_solution[1:366, "U_m1"]) +sum(O_solution[1:366, "U_f1"]) +sum(O_solution[1:366, "CR_m1"]) +sum(O_solution[1:366, "CR_f1"]) +sum(O_solution[1:366, "CS_m1"]) +sum(O_solution[1:366, "CS_f1"]) 
  results_epi_ac[9, i] <- O_solution[366, "RR_m1"] + O_solution[366, "RR_f1"]+O_solution[366, "RS_m1"] + O_solution[366, "RS_f1"]
  results_epi_ac[10, i] <- O_solution[366, "discharge"]
  
  results_epi_MRSAprev[, i] <- O_solution2$CR_m1 + O_solution2$CR_f1 + O_solution2$IMR_m1 + O_solution2$IMR_f1 + O_solution2$ISR_m1 + O_solution2$ISR_f1
  results_epi_MRSAinfe[, i]  <- O_solution2$IMR_m1 + O_solution2$IMR_f1 + O_solution2$ISR_m1 + O_solution2$ISR_f1
  results_epi_MRSAdead[, i]  <- O_solution2$DR_m1 + O_solution2$DR_f1
  results_epi_MRSAinfe_all[, i]  <-  O_solution2$IMR_m1 + O_solution2$IMR_f1 + O_solution2$ISR_m1 + O_solution2$ISR_f1 + O_solution2$IMS_m1 + O_solution2$IMS_f1 + O_solution2$ISS_m1 + O_solution2$ISS_f1
  results_epi_MRSAdead_all[, i]  <- O_solution2$DR_m1 + O_solution2$DR_f1 + O_solution2$DS_m1 + O_solution2$DS_f1
  
  #Economics
  #Store econ results per strategy
  results_econ[1, i] <- O_solution[366, "cost"]
  results_econ[2, i] <- (results_epi_ac[8, i]*0.92)+(sum(O_solution[1:366, "ISR_m1"]))*0.58+(sum(O_solution[1:366, "ISR_f1"]))*0.58+((sum(O_solution[1:366, "IMR_m1"])+sum(O_solution[1:366, "IMR_f1"]))*0.64)+((sum(O_solution[1:366, "IMS_m1"])+sum(O_solution[1:366, "IMS_f1"]))*0.64)+(sum(O_solution[1:366, "ISS_m1"]))*0.58+(sum(O_solution[1:366, "ISS_f1"]))*0.58 +(O_solution[366, "RR_m1"] + O_solution[366, "RR_f1"]+O_solution[366, "RS_m1"] + O_solution[366, "RS_f1"])*0.92 +results_epi_ac[10, i]*0.92
  results_econ[3, i] <- results_econ[2, i]/O_solution[366,"new_admin"] #check if usage is appropriate
}
#Compute ICER per strategy
results_icer <- matrix(nrow = 13, ncol = length(model_list))
results_icer[1,1] <- 0
results_icer[2,1] <- results_econ[1, 1]
results_icer[3,1] <- results_econ[2, 1]
results_icer[4,1] <- results_epi_ac[1, 1]
results_icer[5,1] <- results_epi_ac[2, 1]
results_icer[6,1] <- results_epi_ac[3, 1]
results_icer[7,1] <- results_epi_ac[4, 1]
results_icer[8,1] <- results_epi_ac[5, 1]
results_icer[9,1] <- results_epi_ac[9, 1]
results_icer[10,1] <- results_epi_ac[8, 1]
results_icer[11,1] <- results_epi_ac[10, 1]
results_icer[12,1] <- results_epi_ac[10, 1]+  results_epi_ac[9, 1]+results_epi_ac[8, 1]+results_epi_ac[4, 1]
results_icer[13,1] <- results_epi_ac[6, 1]
for (i in 2:length(model_list)) {
  results_icer[1, i] <- 0
  results_icer[2, i] <-  results_econ[1, i]
  results_icer[3, i] <-  results_econ[2, i]
  results_icer[4, i]<- results_epi_ac[1, i]
  results_icer[5, i]<- results_epi_ac[2, i]
  results_icer[6, i]<- results_epi_ac[3, i]
  results_icer[7, i]<- results_epi_ac[4, i]
  results_icer[8, i]<- results_epi_ac[5, i]
  results_icer[9, i]<- results_epi_ac[9, i]
  results_icer[10, i]<- results_epi_ac[8, i]
  results_icer[11, i]<- results_epi_ac[10, i]
  results_icer[12, i]<- results_epi_ac[10, i]+  results_epi_ac[9, i]+results_epi_ac[8, i]+results_epi_ac[4, i]
  results_icer[3, i] <- results_econ[2, i] + ifelse((results_icer[12, 1] - results_icer[12, i]) > 0, (results_icer[12, 1] - results_icer[12, i]) * 0.92, 0)
  results_icer[1, i] <- (results_icer[2, i]-results_icer[2, 1])/(results_icer[3, i]-results_icer[3, 1])
  results_icer[13,i] <- results_epi_ac[6, i]
  
}
transposed_icer <- t(results_icer)
transposed_icerdf <- as.data.frame(transposed_icer)
colnames(transposed_icerdf) <- c("ICER", "Costs", "QALYs","ARB colonisation","ARB infections", "ARB deaths", "Total infections","Total deaths","Total Recovered","Total U,CR,CS", "Total discharge","total population")  # Add more names as needed
rownames(transposed_icerdf) <- c("Do-nothing", "Strategy2", "Strategy3", "Strategy4", "Strategy5", "Strategy6", "Strategy7", "Strategy8", "Strategy9", "Strategy10","Strategy11", "Strategy12","Strategy13","Strategy14","Strategy15", "Strategy16")
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis/0_Figures_model")
write.csv(transposed_icerdf, "matrix_resMRSA_UB.csv", row.names = TRUE)
#####
Upper_boundMRSA_inf<-(results_epi_ac[2, ]/results_epi_ac[6, ])*100
Upper_boundMRSA_col <- ((results_epi_ac[1, ]-results_epi_ac[2,])/results_epi_ac[6, ])*100
Upper_boundMRSA_death<- results_epi_ac[3,]
#LOWERBOUNDestimatesMRSA#####
#  Lower bound: 0.06131473 R  S 0.2396453
model_list <- c("ARB_model_1ch_do_nothing", "ARB_model_1ch_td_newadm","ARB_model_1ch2_td_newadm","ARB_model_1pcr_td_newadm","ARB_model_1ch_tiso_newadm","ARB_model_1ch2_tiso_newadm","ARB_model_1pcr_tiso_newadm",
                "ARB_model_1ch_td_newadmHR_m","ARB_model_1ch2_td_newadmHR_m","ARB_model_1_pcr_td_newadmHR_m","ARB_model_1ch_td_newadmHR_f","ARB_model_1ch2_td_newadmHR_f","ARB_model_1_pcr_td_newadmHR_f","ARB_model_1preE_newadm_all","ARB_model_1preE_newadm_m","ARB_model_1preE_newadm_f")  # Add last interventions with high-risk groups and pre-emptive isolation
results_epi_ac <- matrix(nrow = 10, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
results_epi_MRSAprev <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAinfe <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAdead <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAinfe_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_MRSAdead_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_uncolonised <- matrix(nrow = 366, ncol = length(model_list))

results_econ <- matrix(nrow = 3, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
for (i in 1:length(model_list)) {
  #Run initial conditions first:
  times <- seq(from=0, to=365, by = 1)  # Simulate over a year
  #Staphylococcus aureus params [MRSA/MSSA]####### 
  # ----------------------------------------#
  zeta3_m_p1 <- 0.231  # mortality rate from IMS  , male
  zeta3_f_p1 = zeta3_m_p1*2.07  # mortality rate from IMS  , female
  zeta1_m_p1 = zeta3_m_p1*1.01  # mortality rate from IMR, male
  zeta1_f_p1 = zeta3_m_p1*1.22  # mortality rate from IMR, female
  zeta2_m_p1 = zeta3_m_p1*2.32  # mortality rate from ISR, male
  zeta2_f_p1 = zeta3_m_p1*2.25  # mortality rate from ISR, female
  zeta4_m_p1 = zeta3_m_p1*1.10  # mortality rate from ISS , male
  zeta4_f_p1 = zeta3_m_p1*2.34  # mortality rate from ISS , female
  #(omega1_d_m_p1*zeta1_m_p1*IMR_m1)+(omega2_d_m_p1*zeta2_m_p1*ISR_m1)
  kappa_p1 = (35.9-6.2)/(49.5-4.85)
  parameters1 <- c(
    delta1_p1 = 0.0016,  
    delta2_p1 = 0.0016,  
    Disch_U_f_p1 = 1/6,
    Disch_U_m_p1 = 1/6,
    Disch_CR_f_p1 = 1/6,
    Disch_CR_m_p1 = 1/6,
    Disch_CS_f_p1 = 1/6,
    Disch_CS_m_p1 = 1/6,
    mu0_p1 = 0.52,  # % of women among U
    mu1_p1 = 0.52,  # % of women among CR
    mu2_p1 = 0.52,  # % of women among CS
    mu3_p1 = 0.3165,  # % of women among IMR
    mu4_p1 = 0.3889,  # % of women among ISR
    mu5_p1 = 0.3938,  # % of women among IMS
    mu6_p1 = 0.4487,  # % of women among IMS
    #Exposure to anbiotics; treatment amongst susceptible populations
    psi_m_p1 = 0.1474, # % of  individuals exposed  to vancomycin/penicillin among males
    psi_w_p1 = 0.184,  # % of  individuals exposed to vancomycin/penicillin among females
    #Percentage of people under treatment for MRSA decolonisation
    psi_mtr_p1=0.109, #IEAT is 1.35 times higher among CR; hence if treatment is psi_m_p1; psi_mtr_p1=psi_m_p1/1.35
    psi_wtr_p1=0.136,
    #Fitness cost. c reduces the transmission rate among resistant strains [1/unit time] [%].
    c_p1=0.09,
    #Progression to the development of infection from colonisation among CR and CS states. 
    beta1_m_p1 = (1/21)*0.26,  # inverse of LOS plus progression from colonisation to infection among males CR
    beta2_m_p1 = (1/11)*0.099,  # inverse of LOS plus progression from colonisation to infection among males CS
    beta1_f_p1 = (1/29)*0.26, # inverse of LOS plus progression from colonisation to infection among females CR
    beta2_f_p1 = (1/14)*0.099,  # inverse of LOS plus progression from colonisation to infection among females CS
    #Natural clearance of mild and severe infections among CR and CS states, respectively [1/unit time] [%].
    gamma1_p1 = 0.001,  # Natural clearance among mild infections R
    gamma2_p1 = 0.001,  # Natural clearance among severe infections R
    gamma3_p1 = 0.001,  # Natural clearance among mild infections S
    gamma4_p1 = 0.001, # Natural clearance among severe infections S
    #Mean time of infection considering length of hospital stays [1/length of hospital stay]. 
    omega1_d_m_p1 = (1/20), #IMR patients who died, male
    omega1_r_m_p1 = (1/23), #IMR patients who recovered, male
    omega1_d_f_p1 = (1/13), #IMR patients who died, female
    omega1_r_f_p1 = (1/26), #IMR patients who recovered, female
    omega2_d_m_p1 = (1/11),  #ISR patients who died, male
    omega2_r_m_p1 = (1/18),  #ISR patients who recovered, male
    omega2_d_f_p1 = (1/14), #ISR patients who died, female 
    omega2_r_f_p1 = (1/19), #ISR patients who recovered, female
    omega3_d_m_p1 = (1/12), #IMS patients who died, male
    omega3_r_m_p1 = (1/11), #IMS patients who recovered, male
    omega3_d_f_p1 = (1/20), #IMS patients who died, female
    omega3_r_f_p1 = (1/16), #IMS patients who recovered, female
    omega4_d_m_p1 = (1/11), #ISS patients who died, male
    omega4_r_m_p1 = (1/19), #ISS patients who recovered, male
    omega4_d_f_p1 = (1/14), #ISS patients who died, female
    omega4_r_f_p1 = (1/17), #ISS patients who recovered, female
    #Percentage of inpatients with CR or CS, respectively, progressing to severe infection in intensive care units [1/unit time] [%].
    alpha1_f_p1 = 0.4196, #% patients with CR progressing to severe infection, males
    alpha2_f_p1 = 0.3438, #% patients with CS progressing to severe infection, males
    alpha1_m_p1 = 0.3517, #% patients with CR progressing to severe infection, females
    alpha2_m_p1 = 0.2826, #% patients with CS progressing to severe infection, females
    #Progression from mild to severe infection from IMR and IMS, respectively [1/unit time] [%].
    epsilon1_p1 = 0.01, #progression from IMR to ISR
    epsilon2_p1 = 0.01, #progression from IMS to ISS
    #Mortality rates from infection. ζ1 and ζ2 are mortality rates from mild and severe resistant infections, respectively. ζ3 and ζ4 are from mild and severe susceptible infections, respectively [1/unit time] [%].
    zeta3_m_p1 = 0.231,  # mortality rate from IMS  , male
    zeta3_f_p1 = 0.231*2.07,  # mortality rate from IMS  , female
    zeta1_m_p1 = 0.231*1.01,  # mortality rate from IMR, male
    zeta1_f_p1 = 0.231*1.22,  # mortality rate from IMR, female
    zeta2_m_p1 = 0.231*2.32,  # mortality rate from ISR, male
    zeta2_f_p1 = 0.231*2.25,  # mortality rate from ISR, female
    zeta4_m_p1 = 0.231*1.10,  # mortality rate from ISS , male
    zeta4_f_p1 = 0.231*2.34,  # mortality rate from ISS , female
    #Recovery rates from infection, including IMR, ISR, IMS and ISS due to treatment received [1/unit time] [%].
    nu1_m_p1 = (1-zeta1_m_p1), #recovery rates from IMR, males
    nu1_f_p1 = (1-zeta1_f_p1), #recovery rates from IMR, females
    nu2_m_p1 = (1-zeta2_m_p1), #recovery rates from ISR, males
    nu2_f_p1 = (1-zeta2_f_p1), #recovery rates from ISR, females
    nu3_m_p1 = (1-zeta3_m_p1), #recovery rates from IMS, males
    nu3_f_p1 = (1-zeta3_f_p1), #recovery rates from IMS, females
    nu4_m_p1 =  (1-zeta4_m_p1), #recovery rates from ISS, males
    nu4_f_p1 =  (1-zeta4_f_p1),#recovery rates from ISS, females
    #Constant background rate that captures transmission from non-human sources, horizontal transmission, or de novo emergence [1/unit time] [number].
    b_p1 = 0.01,
    # Percentage of people with resistant infections receiving inappropriate empirical antibiotic treatment [1/unit time] [%].
    phi_m_p1 = 0.459,  # 
    phi_f_p1 = 0.413,  # 
    #Factor of burden associated with inappropriate empirical antibiotic treatment and increased ICU admission among resistant infections [1/unit time] [%].
    pi_p1= 1.35,
    #Transmission parameter {update this correspondingly after calibrating it with real data}
    #tau_p1= 0.03461113, #estimated 
    tau_p1=0.11,
    #community-acquired infection upon hospital admission rate
    caIha_p1= 0.001,
    ###
    
    #percentage of people tested
    test_p1=0.20, 
    HR_perc1=0.2,
    or_HR_scenar1_a=1.04,
    or_HR_scenarMen_a=1.37,
    ### ### ### ### ###
    #sensitivity chrom_1
    sens_chrom_a=0.826,
    #sensitivity chrom_1
    sens_chrom2_a=0.622, 
    #sensitivity chrom_1
    sens_pcr_a=0.881,
    #turnaround chrom_1
    turn_chrom_a=3,
    #turnaround chrom_1
    turn_chrom2_a=2,
    #turnaround pcr_1
    turn_pcr_a=1,
    #isolation contact precaution transmission reduction
    reduc_conpre_a=0.365,
    #efficiency decolonisation
    eff_decol_a= 0.53,
    #effect on self-infection decolonisation
    eff_decol_selfi_a=0.32,
    #Turnaround decolonisation program in days
    turnaround_decol_a=5,
    ##
    #costs wards
    c_general_ward= 50,
    c_intermediate_ward=92,
    c_icu_ward=218,
    c_decol_1pd=6.5,
    c_isolation=42.3,
    c_chrom=10.2,
    c_chrom2=13.6, 
    c_pcr=33,
    c_bc=16.9,
    #utilities
    u_healthy=0.92,
    u_icu=0.58,
    u_gw=0.64,
    u_recovICU=0.74
  )
  
  
  
  
  # ----------------------------------------#
  #Staphylococcus aureus states [MRSA/MSSA]####### 
  # ----------------------------------------#
  N <- 1000  # Total population size
  # Initial conditions (population sizes in each group)
  U_m10 <- 0.7 * N *(1-0.52)
  CR_m10 <- 0.0613 * N *(1-0.52)
  CS_m10 <-  0.2396* N *(1-0.52)
  IMR_m10 <- 0.26 * CR_m10* (1-0.4196)
  ISR_m10 <- 0.26 * CR_m10* 0.4196
  IMS_m10 <- 0.099 * CS_m10*(1-0.3438)
  ISS_m10 <- 0.099 * CS_m10*0.3438
  RR_m10 <-0
  RS_m10  <-0
  DR_m10  <-0
  DS_m10 <-0
  N_to0 <- 1050
  utility_to0<-0
  cost_to0<-0
  new_admin0<-0
  discharge<-0
  
  U_f10 <- 0.7 * N *0.52
  CR_f10 <- 0.0613 * N *0.52
  CS_f10 <-  0.2396* N *0.52
  IMR_f10 <- 0.26 * CR_f10* (1-0.3517)
  ISR_f10 <- 0.26 * CR_f10*0.3517
  IMS_f10 <- 0.099 * CS_f10*(1-0.2826)
  ISS_f10 <- 0.099 * CS_f10*0.2826
  RR_f10 <-0
  RS_f10  <-0
  DR_f10  <-0
  DS_f10 <-0
  
  N_0m10<-  U_m10 + CR_m10 + CS_m10 + IMR_m10 + ISR_m10 + IMS_m10 + ISS_m10 +  RR_m10 + RS_m10 + DR_m10 +DS_m10
  N_0f10<-  U_f10 + CR_f10 + CS_f10 + IMR_f10 + ISR_f10 + IMS_f10 + ISS_f10 +  RR_f10 + RS_f10 + DR_f10 +DS_f10
  
  state1 <- c(U_m1 = U_m10, CR_m1=CR_m10, CS_m1= CS_m10, IMR_m1= IMR_m10, ISR_m1=ISR_m10, IMS_m1= IMS_m10, ISS_m1= ISS_m10, RR_m1= RR_m10, RS_m1=RS_m10, DR_m1= DR_m10, DS_m1=DS_m10,
              U_f1 = U_f10, CR_f1=CR_f10, CS_f1= CS_f10, IMR_f1= IMR_f10, ISR_f1=ISR_f10, IMS_f1= IMS_f10, ISS_f1= ISS_f10, RR_f1= RR_f10, RS_f1=RS_f10, DR_f1= DR_f10, DS_f1=DS_f10, N_to=N_to0, utility=utility_to0, cost=cost_to0, new_admin=new_admin0, discharge=discharge0)
  
  N_orig1<-N_0m10 + N_0f10
  N_tdif <- N_orig1
  #####
  # Construct the function name from the model_list string
  model_func_name <- get(model_list[i])
  # Solve the ODE using the 'deSolve' package's ode function
  O_solution <- ode(y = state1, times = times, func = model_func_name, parms = parameters1, method = "rk4")
  O_solution2 <- as.data.frame(O_solution)
  O_solution <- as.data.frame(O_solution)
  # Process the ODE solution: sum all values reported in the first 350 rows for your variable of interest
  results_epi_ac[1, i] <- sum(O_solution[1:366, "CR_m1"]) + sum(O_solution[1:366, "CR_f1"]) + sum(O_solution[1:366, "IMR_m1"]) + sum(O_solution[1:366, "IMR_f1"])+ sum(O_solution[1:366, "ISR_m1"])+sum(O_solution[1:366, "ISR_f1"])
  results_epi_ac[2, i] <- sum(O_solution[1:366, "IMR_m1"]) + sum(O_solution[1:366, "IMR_f1"])+ sum(O_solution[1:366, "ISR_m1"])+sum(O_solution[1:366, "ISR_f1"])
  results_epi_ac[3, i] <- (O_solution[366, "DR_f1"]) + (O_solution[366, "DR_m1"]) 
  results_epi_ac[4, i] <- sum(O_solution[1:366, "IMR_m1"]) + sum(O_solution[1:366, "IMR_f1"])+ sum(O_solution[1:366, "ISR_m1"])+sum(O_solution[1:366, "ISR_f1"]) + sum(O_solution[1:366, "IMS_m1"]) + sum(O_solution[1:366, "IMS_f1"])+ sum(O_solution[1:366, "ISS_m1"])+sum(O_solution[1:366, "ISS_f1"])
  results_epi_ac[5, i] <- (O_solution[366, "DR_f1"]) + (O_solution[366, "DR_m1"]) + (O_solution[366, "DS_f1"]) + (O_solution[366, "DS_m1"]) 
  results_epi_ac[6, i] <- sum(O_solution[366, "new_admin"]) 
  results_epi_ac[7, i] <- sum(O_solution[1:366, "U_m1"]) +sum(O_solution[1:366, "U_f1"]) 
  results_epi_ac[8, i] <- sum(O_solution[1:366, "U_m1"]) +sum(O_solution[1:366, "U_f1"]) +sum(O_solution[1:366, "CR_m1"]) +sum(O_solution[1:366, "CR_f1"]) +sum(O_solution[1:366, "CS_m1"]) +sum(O_solution[1:366, "CS_f1"]) 
  results_epi_ac[9, i] <- O_solution[366, "RR_m1"] + O_solution[366, "RR_f1"]+O_solution[366, "RS_m1"] + O_solution[366, "RS_f1"]
  results_epi_ac[10, i] <- O_solution[366, "discharge"]
  
  results_epi_MRSAprev[, i] <- O_solution2$CR_m1 + O_solution2$CR_f1 + O_solution2$IMR_m1 + O_solution2$IMR_f1 + O_solution2$ISR_m1 + O_solution2$ISR_f1
  results_epi_MRSAinfe[, i]  <- O_solution2$IMR_m1 + O_solution2$IMR_f1 + O_solution2$ISR_m1 + O_solution2$ISR_f1
  results_epi_MRSAdead[, i]  <- O_solution2$DR_m1 + O_solution2$DR_f1
  results_epi_MRSAinfe_all[, i]  <-  O_solution2$IMR_m1 + O_solution2$IMR_f1 + O_solution2$ISR_m1 + O_solution2$ISR_f1 + O_solution2$IMS_m1 + O_solution2$IMS_f1 + O_solution2$ISS_m1 + O_solution2$ISS_f1
  results_epi_MRSAdead_all[, i]  <- O_solution2$DR_m1 + O_solution2$DR_f1 + O_solution2$DS_m1 + O_solution2$DS_f1
  
  #Economics
  #Store econ results per strategy
  results_econ[1, i] <- O_solution[366, "cost"]
  results_econ[2, i] <- (results_epi_ac[8, i]*0.92)+(sum(O_solution[1:366, "ISR_m1"]))*0.58+(sum(O_solution[1:366, "ISR_f1"]))*0.58+((sum(O_solution[1:366, "IMR_m1"])+sum(O_solution[1:366, "IMR_f1"]))*0.64)+((sum(O_solution[1:366, "IMS_m1"])+sum(O_solution[1:366, "IMS_f1"]))*0.64)+(sum(O_solution[1:366, "ISS_m1"]))*0.58+(sum(O_solution[1:366, "ISS_f1"]))*0.58 +(O_solution[366, "RR_m1"] + O_solution[366, "RR_f1"]+O_solution[366, "RS_m1"] + O_solution[366, "RS_f1"])*0.92 +results_epi_ac[10, i]*0.92
  results_econ[3, i] <- results_econ[2, i]/O_solution[366,"new_admin"] #check if usage is appropriate
}
#Compute ICER per strategy
results_icer <- matrix(nrow = 13, ncol = length(model_list))
results_icer[1,1] <- 0
results_icer[2,1] <- results_econ[1, 1]
results_icer[3,1] <- results_econ[2, 1]
results_icer[4,1] <- results_epi_ac[1, 1]
results_icer[5,1] <- results_epi_ac[2, 1]
results_icer[6,1] <- results_epi_ac[3, 1]
results_icer[7,1] <- results_epi_ac[4, 1]
results_icer[8,1] <- results_epi_ac[5, 1]
results_icer[9,1] <- results_epi_ac[9, 1]
results_icer[10,1] <- results_epi_ac[8, 1]
results_icer[11,1] <- results_epi_ac[10, 1]
results_icer[12,1] <- results_epi_ac[10, 1]+  results_epi_ac[9, 1]+results_epi_ac[8, 1]+results_epi_ac[4, 1]
results_icer[13,1] <- results_epi_ac[6, 1]
for (i in 2:length(model_list)) {
  results_icer[1, i] <- 0
  results_icer[2, i] <-  results_econ[1, i]
  results_icer[3, i] <-  results_econ[2, i]
  results_icer[4, i]<- results_epi_ac[1, i]
  results_icer[5, i]<- results_epi_ac[2, i]
  results_icer[6, i]<- results_epi_ac[3, i]
  results_icer[7, i]<- results_epi_ac[4, i]
  results_icer[8, i]<- results_epi_ac[5, i]
  results_icer[9, i]<- results_epi_ac[9, i]
  results_icer[10, i]<- results_epi_ac[8, i]
  results_icer[11, i]<- results_epi_ac[10, i]
  results_icer[12, i]<- results_epi_ac[10, i]+  results_epi_ac[9, i]+results_epi_ac[8, i]+results_epi_ac[4, i]
  results_icer[3, i] <- results_econ[2, i] + ifelse((results_icer[12, 1] - results_icer[12, i]) > 0, (results_icer[12, 1] - results_icer[12, i]) * 0.92, 0)
  results_icer[1, i] <- (results_icer[2, i]-results_icer[2, 1])/(results_icer[3, i]-results_icer[3, 1])
  results_icer[13,i] <- results_epi_ac[6, i]
  
}
transposed_icer <- t(results_icer)
transposed_icerdf <- as.data.frame(transposed_icer)
colnames(transposed_icerdf) <- c("ICER", "Costs", "QALYs","ARB colonisation","ARB infections", "ARB deaths", "Total infections","Total deaths","Total Recovered","Total U,CR,CS", "Total discharge","total population")  # Add more names as needed
rownames(transposed_icerdf) <- c("Do-nothing", "Strategy2", "Strategy3", "Strategy4", "Strategy5", "Strategy6", "Strategy7", "Strategy8", "Strategy9", "Strategy10","Strategy11", "Strategy12","Strategy13","Strategy14","Strategy15", "Strategy16")
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis/0_Figures_model")
write.csv(transposed_icerdf, "matrix_resMRSA_lB.csv", row.names = TRUE)
#####
lower_boundMRSA_inf<-(results_epi_ac[2, ]/results_epi_ac[6, ])*100
lower_boundMRSA_col <- ((results_epi_ac[1, ]-results_epi_ac[2,])/results_epi_ac[6, ])*100
lower_boundMRSA_death<- results_epi_ac[3,]

#-------------------------------------------------#
# CRE:
#-------------------------------------------------#
#Upper bound: 0.21954  R. 0.34046  S
#  Lower bound: 0.06945999 R 0.49054 S
#UpperBOUND ESTIMATIONs####
model_list <- c("ARB_model_2ch_do_nothing", "ARB_model_2ch_td_newadm","ARB_model_2ch2_td_newadm","ARB_model_2pcr_td_newadm","ARB_model_2ch_tiso_newadm","ARB_model_2ch2_tiso_newadm","ARB_model_2pcr_tiso_newadm",
                "ARB_model_2ch_td_newadmHR_m","ARB_model_2ch2_td_newadmHR_m","ARB_model_2_pcr_td_newadmHR_m","ARB_model_2ch_td_newadmHR_f","ARB_model_2ch2_td_newadmHR_f","ARB_model_2_pcr_td_newadmHR_f","ARB_model_2preE_newadm_all","ARB_model_2preE_newadm_m","ARB_model_2preE_newadm_f")  # Add last interventions with high-risk groups and pre-emptive isolation
results_epi_ac <- matrix(nrow = 10, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
results_epi_CREprev <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREinfe <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREdead <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREinfe_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREdead_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_uncolonised <- matrix(nrow = 366, ncol = length(model_list))
results_econ <- matrix(nrow = 3, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
for (i in 1:length(model_list)) {
  #Run initial conditions first:
  times <- seq(from=0, to=365, by = 1)  # Simulate over a year
  #Enterobacterales[CRE/CSE]####### 
  # ----------------------------------#
  N <- 1000  # Total population size
  # Initial conditions (population sizes in each group)
  U_m20 <- 0.44 * N *(1-0.52)
  CR_m20 <- 0.21954 * N *(1-0.52)
  CS_m20 <-  0.34046* N *(1-0.52)
  IMR_m20 <- 0.09 * CR_m20* (1-0.4283)
  ISR_m20 <- 0.09 * CR_m20* 0.4283
  IMS_m20 <- 0.04 * CS_m20*(1-0.3548)
  ISS_m20 <- 0.04 * CS_m20*0.3548
  RR_m20 <-0
  RS_m20  <-0
  DR_m20  <-0
  DS_m20 <-0
  N_to2<-1050
  utility_to0<-0
  cost_to0<-0
  new_admin0<-0
  discharge0<- 0
  
  U_f20 <- 0.44 * N *0.52
  CR_f20 <- 0.21954 * N *0.52
  CS_f20 <-  0.34046* N *0.52
  IMR_f20 <- 0.09 * CR_f20* (1-0.4538)
  ISR_f20 <- 0.09 * CR_f20*0.4538
  IMS_f20 <- 0.04 * CS_f20*(1-0.3832)
  ISS_f20 <- 0.04 * CS_f20*0.3832
  RR_f20 <-0
  RS_f20  <-0
  DR_f20  <-0
  DS_f20 <-0
  
  N_0m20<-  U_m20 + CR_m20 + CS_m20 + IMR_m20 + ISR_m20 + IMS_m20 + ISS_m20 +  RR_m20 + RS_m20 + DR_m20 +DS_m20
  N_0f20<-  U_f20 + CR_f20 + CS_f20 + IMR_f20 + ISR_f20 + IMS_f20 + ISS_f20 +  RR_f20 + RS_f20 + DR_f20 +DS_f20
  N_to0<- N_0m20 + N_0f20
  
  state2 <- c(U_m2 = U_m20, CR_m2=CR_m20, CS_m2= CS_m20, IMR_m2= IMR_m20, ISR_m2=ISR_m20, IMS_m2= IMS_m20, ISS_m2= ISS_m20, RR_m2= RR_m20, RS_m2=RS_m20, DR_m2= DR_m20, DS_m2=DS_m20,
              U_f2 = U_f20, CR_f2=CR_f20, CS_f2= CS_f20, IMR_f2= IMR_f20, ISR_f2=ISR_f20, IMS_f2= IMS_f20, ISS_f2= ISS_f20, RR_f2= RR_f20, RS_f2=RS_f20, DR_f2= DR_f20, DS_f2=DS_f20, N_to2=N_to0, utility=utility_to0, cost=cost_to0, new_admin=new_admin0, discharge=discharge0)
  N_orig2<-N_0m20 + N_0f20
  N_tdif <- N_orig2
  Nt2_spec2<-1030
  
  #Enterobacterales[CRE/CSE]params####### 
  # ----------------------------------------#
  
  zeta3_m_p2 = 0.228  # mortality rate from IMS  , male
  zeta3_f_p2 = zeta3_m_p2*0.81  # mortality rate from IMS  , female
  zeta1_m_p2 = zeta3_m_p2*1.80  # mortality rate from IMR, male
  zeta1_f_p2 = zeta3_m_p2*0.55  # mortality rate from IMR, female
  zeta2_m_p2 = zeta3_m_p2*1.30  # mortality rate from ISR, male
  zeta2_f_p2 = zeta3_m_p2*2.40  # mortality rate from ISR, female
  zeta4_m_p2 = zeta3_m_p2*1.62  # mortality rate from ISS , male
  zeta4_f_p2 = zeta3_m_p2*2.23  # mortality rate from ISS , female
  kappa_p2 = (35.9-6.2)/(49.5-4.85)
  parameters2 <- c(
    #clearance -natural- parameters
    delta1_p2 = 0.001,  # Clearance value from CR state
    delta2_p2 = 0.001,  # Clearance value from CS state
    #Dicharge rates from uncolonised and colonised
    Disch_U_f_p2 = 1/6,
    Disch_U_m_p2 = 1/6,
    Disch_CR_f_p2 = 1/6,
    Disch_CR_m_p2 = 1/6,
    Disch_CS_f_p2 = 1/6,
    Disch_CS_m_p2 = 1/6,
    #Proportion of women among specific populations (i.e., U, CR, CS, IMR, ISR, IMS, and ISS) [1/unit time] [%]
    mu0_p2 = 0.52,  # % of women among U
    mu1_p2 = 0.52,  # % of women among CR
    mu2_p2 = 0.52,  # % of women among CS
    mu3_p2 = 0.3394,  # % of women among IMR
    mu4_p2 = 0.3280,  # % of women among ISR
    mu5_p2 = 0.50,  # % of women among IMS
    mu6_p2 = 0.44,  # % of women among IMS
    
    #Exposure to anbiotics 
    psi_m_p2 = 0.2225, # % of  individuals exposed  to vancomycin/penicillin among males
    psi_w_p2 = 0.2026,  # % of  individuals exposed to vancomycin/penicillin among males
    
    #Percentage of people under treatment for CRE decolonisation
    psi_mtr_p2=0.2225/1.02,
    psi_wtr_p2=0.2026/1.02,
    
    #Fitness cost. c reduces the transmission rate among resistant strains [1/unit time] [%].
    c_p2=(1-0.927),
    
    #Progression to the development of infection from colonisation among CR and CS states. 
    beta1_m_p2 = (1/22)*0.213, # inverse of LOS plus progression from colonisation to infection among males CR
    beta2_m_p2 = (1/20)*0.034,  # inverse of LOS plus progression from colonisation to infection among males CS
    beta1_f_p2 = (1/27)*0.213, # inverse of LOS plus progression from colonisation to infection among females CR
    beta2_f_p2 = (1/17)*0.034,  # inverse of LOS plus progression from colonisation to infection among females CS
    
    #Natural clearance of mild and severe infections among CR and CS states, respectively [1/unit time] [%].
    gamma1_p2 = 0.001,  # Natural clearance among mild infections R
    gamma2_p2 = 0.001,  # Natural clearance among severe infections R
    gamma3_p2 = 0.001,  # Natural clearance among mild infections S
    gamma4_p2 = 0.001, # Natural clearance among severe infections S
    
    #Mean time of infection considering length of hospital stays [1/length of hospital stay]. 
    omega1_d_m_p2 =(1/21), #IMR patients who died, male
    omega1_r_m_p2 =(1/26), #IMR patients who recovered, male
    omega1_d_f_p2 =(1/31), #IMR patients who died, female
    omega1_r_f_p2 =(1/30), #IMR patients who recovered, female
    omega2_d_m_p2 =(1/7), #ISR patients who died, male
    omega2_r_m_p2 =(1/20),  #ISR patients who recovered, male
    omega2_d_f_p2 = (1/20), #ISR patients who died, female 
    omega2_r_f_p2 =(1/23), #ISR patients who recovered, female
    omega3_d_m_p2 =(1/12), #IMS patients who died, male
    omega3_r_m_p2 =(1/20), #IMS patients who recovered, male
    omega3_d_f_p2 =(1/10), #IMS patients who died, female
    omega3_r_f_p2 =(1/18), #IMS patients who recovered, female
    omega4_d_m_p2 =(1/11), #ISS patients who died, male
    omega4_r_m_p2 =(1/14), #ISS patients who recovered, male
    omega4_d_f_p2 =(1/9), #ISS patients who died, female
    omega4_r_f_p2 =(1/15), #ISS patients who recovered, female
    
    #Percentage of inpatients with CR or CS, respectively, progressing to severe infection in intensive care units [1/unit time] [%].
    alpha1_f_p2 = 0.4283, #% patients with CR progressing to severe infection, males
    alpha2_f_p2 = 0.3548, #% patients with CS progressing to severe infection, males
    alpha1_m_p2 = 0.4585, #% patients with CR progressing to severe infection, females
    alpha2_m_p2 = 0.3832, #% patients with CS progressing to severe infection, females
    
    #Progression from mild to severe infection from IMR and IMS, respectively [1/unit time] [%].
    epsilon1_p2 = 0.01, #progression from IMR to ISR
    epsilon2_p2 = 0.01, #progression from IMS to ISS
    
    #Mortality rates from infection. ζ1 and ζ2 are mortality rates from mild and severe resistant infections, respectively. ζ3 and ζ4 are from mild and severe susceptible infections, respectively [1/unit time] [%].
    zeta3_m_p2 = 0.228,  # mortality rate from IMS  , male
    zeta3_f_p2 = zeta3_m_p2*0.81,  # mortality rate from IMS  , female
    zeta1_m_p2 = zeta3_m_p2*1.80,  # mortality rate from IMR, male
    zeta1_f_p2 = zeta3_m_p2*0.55,  # mortality rate from IMR, female
    zeta2_m_p2 = zeta3_m_p2*1.30,  # mortality rate from ISR, male
    zeta2_f_p2 = zeta3_m_p2*2.40,  # mortality rate from ISR, female
    zeta4_m_p2 = zeta3_m_p2*1.62,  # mortality rate from ISS , male
    zeta4_f_p2 = zeta3_m_p2*2.23,  # mortality rate from ISS , female
    
    #Recovery rates from infection, including IMR, ISR, IMS and ISS due to treatment received [1/unit time] [%].
    nu1_m_p2 = (1-zeta1_m_p2), #recovery rates from IMR, males
    nu1_f_p2 = (1-zeta1_f_p2), #recovery rates from IMR, females
    nu2_m_p2 = (1-zeta2_m_p2), #recovery rates from ISR, males
    nu2_f_p2 = (1-zeta2_f_p2), #recovery rates from ISR, females
    nu3_m_p2 = (1-zeta3_m_p2), #recovery rates from IMS, males
    nu3_f_p2 = (1-zeta3_f_p2), #recovery rates from IMS, females
    nu4_m_p2 =  (1-zeta4_m_p2), #recovery rates from ISS, males
    nu4_f_p2 =  (1-zeta4_f_p2),#recovery rates from ISS, females
    
    #Constant background rate that captures transmission from non-human sources, horizontal transmission, or de novo emergence [1/unit time] [number].
    b_p2 = 0.01,
    
    # Percentage of people with resistant infections receiving inappropriate empirical antibiotic treatment [1/unit time] [%].
    phi_m_p2 = 0.3782,  # Placeholder value, adjust as needed
    phi_f_p2 = 0.3846,  # Placeholder value, adjust as needed
    
    #Factor of burden associated with inappropriate empirical antibiotic treatment and increased ICU admission among resistant infections [1/unit time] [%].
    pi_p2= 1.02,
    
    #Transmission parameter {update this correspondingly after calibrating it with real data}
    tau_p2= 0.5979826,#0.3986551,
    
    #community-acquired infection upon hospital admission rate
    caIha_p2=0.007,
    
    #percentage of people tested
    test_p2=0.20, 
    HR_perc2=0.2,
    or_HR_scenar1_a=1.04,
    or_HR_scenarMen_a=2.27,
    ### ### ### ### ###
    #sensitivity chrom_1
    sens_chrom_a=0.826,
    #sensitivity chrom_1
    sens_chrom2_a=0.90, 
    #sensitivity chrom_1
    sens_pcr_a=1,
    #turnaround chrom_1
    turn_chrom_a=3,
    #turnaround chrom_1
    turn_chrom2_a=2,
    #turnaround pcr_1
    turn_pcr_a=1,
    #isolation contact precaution transmission reduction
    reduc_conpre_a=0.35,
    #efficiency decolonisation
    eff_decol_a= 0.26,
    #effect on self-infection decolonisation
    eff_decol_selfi_a=0.041,
    #Turnaround decolonisation program in days
    turnaround_decol_a=7,
    ##
    #costs wards
    c_general_ward= 50,
    c_intermediate_ward=92,
    c_icu_ward=218,
    c_decol_1pd=72.88,
    c_isolation=42.3,
    c_chrom=10.2,
    c_chrom2=13.6, 
    c_pcr=33,
    c_bc=16.9,
    #utilities
    u_healthy=0.92,
    u_icu=0.92-0.34,
    u_gw=0.64,
    u_recovICU=0.74
    
  )
  
  #####
  # Construct the function name from the model_list string
  model_func_name <- get(model_list[i])
  # Solve the ODE using the 'deSolve' package's ode function
  O_solution <- ode(y = state2, times = times, func = model_func_name, parms = parameters2, method = "rk4")
  O_solution2 <- as.data.frame(O_solution)
  O_solution <- as.data.frame(O_solution)
  # Process the ODE solution: sum all values reported in the first 350 rows for your variable of interest
  results_epi_ac[1, i] <- sum(O_solution[1:366, "CR_m2"]) + sum(O_solution[1:366, "CR_f2"]) + sum(O_solution[1:366, "IMR_m2"]) + sum(O_solution[1:366, "IMR_f2"])+ sum(O_solution[1:366, "ISR_m2"])+sum(O_solution[1:366, "ISR_f2"])
  results_epi_ac[2, i] <- sum(O_solution[1:366, "IMR_m2"]) + sum(O_solution[1:366, "IMR_f2"])+ sum(O_solution[1:366, "ISR_m2"])+sum(O_solution[1:366, "ISR_f2"])
  results_epi_ac[3, i] <- (O_solution[366, "DR_f2"]) + (O_solution[366, "DR_m2"]) 
  results_epi_ac[4, i] <- sum(O_solution[1:366, "IMR_m2"]) + sum(O_solution[1:366, "IMR_f2"])+ sum(O_solution[1:366, "ISR_m2"])+sum(O_solution[1:366, "ISR_f2"]) + sum(O_solution[1:366, "IMS_m2"]) + sum(O_solution[1:366, "IMS_f2"])+ sum(O_solution[1:366, "ISS_m2"])+sum(O_solution[1:366, "ISS_f2"])
  results_epi_ac[5, i] <- (O_solution[366, "DR_f2"]) + (O_solution[366, "DR_m2"]) + (O_solution[366, "DS_f2"]) + (O_solution[366, "DS_m2"]) 
  results_epi_ac[6, i] <- sum(O_solution[366, "new_admin"]) 
  results_epi_ac[7, i] <- sum(O_solution[1:366, "U_m2"]) +sum(O_solution[1:366, "U_f2"]) 
  results_epi_ac[8, i] <- sum(O_solution[1:366, "U_m2"]) +sum(O_solution[1:366, "U_f2"]) +sum(O_solution[1:366, "CR_m2"]) +sum(O_solution[1:366, "CR_f2"]) +sum(O_solution[1:366, "CS_m2"]) +sum(O_solution[1:366, "CS_f2"]) 
  results_epi_ac[9, i] <- O_solution[366,"RR_m2"] + O_solution[366,"RR_f2"]+O_solution[366,"RS_m2"] + O_solution[366,"RS_f2"]
  results_epi_ac[10, i] <- (O_solution[366, "discharge"]) 
  
  results_epi_CREprev[, i] <- O_solution2$CR_m2 + O_solution2$CR_f2 + O_solution2$IMR_m2 + O_solution2$IMR_f2 + O_solution2$ISR_m2 + O_solution2$ISR_f2
  results_epi_CREinfe[, i]  <- O_solution2$IMR_m2 + O_solution2$IMR_f2 + O_solution2$ISR_m2 + O_solution2$ISR_f2
  results_epi_CREdead[, i]  <- O_solution2$DR_m2 + O_solution2$DR_f2
  results_epi_CREinfe_all[, i]  <-  O_solution2$IMR_m2 + O_solution2$IMR_f2 + O_solution2$ISR_m2 + O_solution2$ISR_f2 + O_solution2$IMS_m2 + O_solution2$IMS_f2 + O_solution2$ISS_m2 + O_solution2$ISS_f2
  results_epi_CREdead_all[, i]  <- O_solution2$DR_m2 + O_solution2$DR_f2 + O_solution2$DS_m2 + O_solution2$DS_f2
  
  #Economics
  #Store econ results per strategy
  results_econ[1, i] <- O_solution[366, "cost"]
  results_econ[2, i] <- (results_epi_ac[8, i]*0.92)+(sum(O_solution[1:366, "ISR_m2"]))*0.58+(sum(O_solution[1:366, "ISR_f2"])*0.58)+((sum(O_solution[1:366, "IMR_m2"])+sum(O_solution[1:366, "IMR_f2"]))*0.64)+((sum(O_solution[1:366, "IMS_m2"])+sum(O_solution[1:366, "IMS_f2"]))*0.64)+(sum(O_solution[1:366, "ISS_m2"]))*0.58+(sum(O_solution[1:366, "ISS_f2"])*0.58)+ results_epi_ac[9, i]*0.92 +results_epi_ac[10, i]*0.92
  results_econ[3, i] <- results_econ[2, i]/O_solution[366,"new_admin"] #check if usage is appropriate
}
#Compute ICER per strategy
results_icer <- matrix(nrow = 13, ncol = length(model_list))
results_icer[1,1] <- 0
results_icer[2,1] <- results_econ[1, 1]
results_icer[3,1] <- results_econ[2, 1]
results_icer[4,1] <- results_epi_ac[1, 1]
results_icer[5,1] <- results_epi_ac[2, 1]
results_icer[6,1] <- results_epi_ac[3, 1]
results_icer[7,1] <- results_epi_ac[4, 1]
results_icer[8,1] <- results_epi_ac[5, 1]
results_icer[9,1] <- results_epi_ac[9, 1]
results_icer[10,1] <- results_epi_ac[8, 1]
results_icer[11,1] <- results_epi_ac[10, 1]
results_icer[12,1] <- results_epi_ac[10, 1]+  results_epi_ac[9, 1]+results_epi_ac[8, 1]+results_epi_ac[4, 1]
results_icer[13,1] <- results_epi_ac[6, 1]
for (i in 2:length(model_list)) {
  results_icer[1, i] <- 0
  results_icer[2, i] <-  results_econ[1, i]
  results_icer[3, i] <-  results_econ[2, i]
  results_icer[4, i]<- results_epi_ac[1, i]
  results_icer[5, i]<- results_epi_ac[2, i]
  results_icer[6, i]<- results_epi_ac[3, i]
  results_icer[7, i]<- results_epi_ac[4, i]
  results_icer[8, i]<- results_epi_ac[5, i]
  results_icer[9, i]<- results_epi_ac[9, i]
  results_icer[10, i]<- results_epi_ac[8, i]
  results_icer[11, i]<- results_epi_ac[10, i]
  results_icer[12, i]<- results_epi_ac[10, i]+  results_epi_ac[9, i]+results_epi_ac[8, i]+results_epi_ac[4, i]
  results_icer[3, i] <- results_econ[2, i] + ifelse((results_icer[12, 1] - results_icer[12, i]) > 0, (results_icer[12, 1] - results_icer[12, i]) * 0.92, 0)
  results_icer[1, i] <- (results_icer[2, i]-results_icer[2, 1])/(results_icer[3, i]-results_icer[3, 1])
  results_icer[13,i] <- results_epi_ac[6, i]
}
transposed_icer <- t(results_icer)
transposed_icerdf <- as.data.frame(transposed_icer)
colnames(transposed_icerdf) <- c("ICER", "Costs", "QALYs","ARB colonisation","ARB infections", "ARB deaths", "Total infections","Total deaths","Total Recovered","Total U,CR,CS", "Total discharge","total population")  # Add more names as needed
rownames(transposed_icerdf) <- c("Do-nothing", "Strategy2", "Strategy3", "Strategy4", "Strategy5", "Strategy6", "Strategy7", "Strategy8", "Strategy9", "Strategy10","Strategy11", "Strategy12","Strategy13","Strategy14","Strategy15", "Strategy16")
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis/0_Figures_model")
write.csv(transposed_icerdf, "matrix_resCRE_UB.csv", row.names = TRUE)

#####
Upper_boundCRE_inf<-(results_epi_ac[2, ]/results_epi_ac[6, ])*100
Upper_boundCRE_col <- ((results_epi_ac[1, ]-results_epi_ac[2,])/results_epi_ac[6, ])*100
Upper_boundCRE_death<- results_epi_ac[3,]
#LowerBOUND ESTIMATION######
model_list <- c("ARB_model_2ch_do_nothing", "ARB_model_2ch_td_newadm","ARB_model_2ch2_td_newadm","ARB_model_2pcr_td_newadm","ARB_model_2ch_tiso_newadm","ARB_model_2ch2_tiso_newadm","ARB_model_2pcr_tiso_newadm",
                "ARB_model_2ch_td_newadmHR_m","ARB_model_2ch2_td_newadmHR_m","ARB_model_2_pcr_td_newadmHR_m","ARB_model_2ch_td_newadmHR_f","ARB_model_2ch2_td_newadmHR_f","ARB_model_2_pcr_td_newadmHR_f","ARB_model_2preE_newadm_all","ARB_model_2preE_newadm_m","ARB_model_2preE_newadm_f")  # Add last interventions with high-risk groups and pre-emptive isolation
results_epi_ac <- matrix(nrow = 10, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
results_epi_CREprev <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREinfe <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREdead <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREinfe_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_CREdead_all <- matrix(nrow = 366, ncol = length(model_list))
results_epi_uncolonised <- matrix(nrow = 366, ncol = length(model_list))
results_econ <- matrix(nrow = 3, ncol = length(model_list))  # Adjust the number of columns based on the number of models/ rows for number of parameters extracted
for (i in 1:length(model_list)) {
  #Run initial conditions first:
  times <- seq(from=0, to=365, by = 1)  # Simulate over a year
  #Enterobacterales[CRE/CSE] states####### 
  # ----------------------------------#
  N <- 1000  # Total population size
  # Initial conditions (population sizes in each group)
  U_m20 <- 0.44 * N *(1-0.52)
  CR_m20 <- 0.069 * N *(1-0.52)
  CS_m20 <-  0.491* N *(1-0.52)
  IMR_m20 <- 0.09 * CR_m20* (1-0.4283)
  ISR_m20 <- 0.09 * CR_m20* 0.4283
  IMS_m20 <- 0.04 * CS_m20*(1-0.3548)
  ISS_m20 <- 0.04 * CS_m20*0.3548
  RR_m20 <-0
  RS_m20  <-0
  DR_m20  <-0
  DS_m20 <-0
  N_to2<-1050
  utility_to0<-0
  cost_to0<-0
  new_admin0<-0
  discharge0<-0
  
  U_f20 <- 0.44 * N *0.52
  CR_f20 <- 0.069 * N *0.52
  CS_f20 <-  0.491* N *0.52
  IMR_f20 <- 0.09 * CR_f20* (1-0.4538)
  ISR_f20 <- 0.09 * CR_f20*0.4538
  IMS_f20 <- 0.04 * CS_f20*(1-0.3832)
  ISS_f20 <- 0.04 * CS_f20*0.3832
  RR_f20 <-0
  RS_f20  <-0
  DR_f20  <-0
  DS_f20 <-0
  
  N_0m20<-  U_m20 + CR_m20 + CS_m20 + IMR_m20 + ISR_m20 + IMS_m20 + ISS_m20 +  RR_m20 + RS_m20 + DR_m20 +DS_m20
  N_0f20<-  U_f20 + CR_f20 + CS_f20 + IMR_f20 + ISR_f20 + IMS_f20 + ISS_f20 +  RR_f20 + RS_f20 + DR_f20 +DS_f20
  N_to0<- N_0m20 + N_0f20
  
  state2 <- c(U_m2 = U_m20, CR_m2=CR_m20, CS_m2= CS_m20, IMR_m2= IMR_m20, ISR_m2=ISR_m20, IMS_m2= IMS_m20, ISS_m2= ISS_m20, RR_m2= RR_m20, RS_m2=RS_m20, DR_m2= DR_m20, DS_m2=DS_m20,
              U_f2 = U_f20, CR_f2=CR_f20, CS_f2= CS_f20, IMR_f2= IMR_f20, ISR_f2=ISR_f20, IMS_f2= IMS_f20, ISS_f2= ISS_f20, RR_f2= RR_f20, RS_f2=RS_f20, DR_f2= DR_f20, DS_f2=DS_f20, N_to2=N_to0, utility=utility_to0, cost=cost_to0, new_admin=new_admin0, discharge=discharge0)
  N_orig2<-N_0m20 + N_0f20
  N_tdif <- N_orig2
  Nt2_spec2<-1030
  #Enterobacterales[CRE/CSE]params####### 
  # ----------------------------------------#
  
  zeta3_m_p2 = 0.228  # mortality rate from IMS  , male
  zeta3_f_p2 = zeta3_m_p2*0.81  # mortality rate from IMS  , female
  zeta1_m_p2 = zeta3_m_p2*1.80  # mortality rate from IMR, male
  zeta1_f_p2 = zeta3_m_p2*0.55  # mortality rate from IMR, female
  zeta2_m_p2 = zeta3_m_p2*1.30  # mortality rate from ISR, male
  zeta2_f_p2 = zeta3_m_p2*2.40  # mortality rate from ISR, female
  zeta4_m_p2 = zeta3_m_p2*1.62  # mortality rate from ISS , male
  zeta4_f_p2 = zeta3_m_p2*2.23  # mortality rate from ISS , female
  kappa_p2 = (35.9-6.2)/(49.5-4.85)
  parameters2 <- c(
    #clearance -natural- parameters
    delta1_p2 = 0.001,  # Clearance value from CR state
    delta2_p2 = 0.001,  # Clearance value from CS state
    #Dicharge rates from uncolonised and colonised
    Disch_U_f_p2 = 1/6,
    Disch_U_m_p2 = 1/6,
    Disch_CR_f_p2 = 1/6,
    Disch_CR_m_p2 = 1/6,
    Disch_CS_f_p2 = 1/6,
    Disch_CS_m_p2 = 1/6,
    #Proportion of women among specific populations (i.e., U, CR, CS, IMR, ISR, IMS, and ISS) [1/unit time] [%]
    mu0_p2 = 0.52,  # % of women among U
    mu1_p2 = 0.52,  # % of women among CR
    mu2_p2 = 0.52,  # % of women among CS
    mu3_p2 = 0.3394,  # % of women among IMR
    mu4_p2 = 0.3280,  # % of women among ISR
    mu5_p2 = 0.50,  # % of women among IMS
    mu6_p2 = 0.44,  # % of women among IMS
    
    #Exposure to anbiotics 
    psi_m_p2 = 0.2225, # % of  individuals exposed  to vancomycin/penicillin among males
    psi_w_p2 = 0.2026,  # % of  individuals exposed to vancomycin/penicillin among males
    
    #Percentage of people under treatment for CRE decolonisation
    psi_mtr_p2=0.2225/1.02,
    psi_wtr_p2=0.2026/1.02,
    
    #Fitness cost. c reduces the transmission rate among resistant strains [1/unit time] [%].
    c_p2=(1-0.927),
    
    #Progression to the development of infection from colonisation among CR and CS states. 
    beta1_m_p2 = (1/22)*0.213, # inverse of LOS plus progression from colonisation to infection among males CR
    beta2_m_p2 = (1/20)*0.034,  # inverse of LOS plus progression from colonisation to infection among males CS
    beta1_f_p2 = (1/27)*0.213, # inverse of LOS plus progression from colonisation to infection among females CR
    beta2_f_p2 = (1/17)*0.034,  # inverse of LOS plus progression from colonisation to infection among females CS
    
    #Natural clearance of mild and severe infections among CR and CS states, respectively [1/unit time] [%].
    gamma1_p2 = 0.001,  # Natural clearance among mild infections R
    gamma2_p2 = 0.001,  # Natural clearance among severe infections R
    gamma3_p2 = 0.001,  # Natural clearance among mild infections S
    gamma4_p2 = 0.001, # Natural clearance among severe infections S
    
    #Mean time of infection considering length of hospital stays [1/length of hospital stay]. 
    omega1_d_m_p2 =(1/21), #IMR patients who died, male
    omega1_r_m_p2 =(1/26), #IMR patients who recovered, male
    omega1_d_f_p2 =(1/31), #IMR patients who died, female
    omega1_r_f_p2 =(1/30), #IMR patients who recovered, female
    omega2_d_m_p2 =(1/7), #ISR patients who died, male
    omega2_r_m_p2 =(1/20),  #ISR patients who recovered, male
    omega2_d_f_p2 = (1/20), #ISR patients who died, female 
    omega2_r_f_p2 =(1/23), #ISR patients who recovered, female
    omega3_d_m_p2 =(1/12), #IMS patients who died, male
    omega3_r_m_p2 =(1/20), #IMS patients who recovered, male
    omega3_d_f_p2 =(1/10), #IMS patients who died, female
    omega3_r_f_p2 =(1/18), #IMS patients who recovered, female
    omega4_d_m_p2 =(1/11), #ISS patients who died, male
    omega4_r_m_p2 =(1/14), #ISS patients who recovered, male
    omega4_d_f_p2 =(1/9), #ISS patients who died, female
    omega4_r_f_p2 =(1/15), #ISS patients who recovered, female
    
    #Percentage of inpatients with CR or CS, respectively, progressing to severe infection in intensive care units [1/unit time] [%].
    alpha1_f_p2 = 0.4283, #% patients with CR progressing to severe infection, males
    alpha2_f_p2 = 0.3548, #% patients with CS progressing to severe infection, males
    alpha1_m_p2 = 0.4585, #% patients with CR progressing to severe infection, females
    alpha2_m_p2 = 0.3832, #% patients with CS progressing to severe infection, females
    
    #Progression from mild to severe infection from IMR and IMS, respectively [1/unit time] [%].
    epsilon1_p2 = 0.01, #progression from IMR to ISR
    epsilon2_p2 = 0.01, #progression from IMS to ISS
    
    #Mortality rates from infection. ζ1 and ζ2 are mortality rates from mild and severe resistant infections, respectively. ζ3 and ζ4 are from mild and severe susceptible infections, respectively [1/unit time] [%].
    zeta3_m_p2 = 0.228,  # mortality rate from IMS  , male
    zeta3_f_p2 = zeta3_m_p2*0.81,  # mortality rate from IMS  , female
    zeta1_m_p2 = zeta3_m_p2*1.80,  # mortality rate from IMR, male
    zeta1_f_p2 = zeta3_m_p2*0.55,  # mortality rate from IMR, female
    zeta2_m_p2 = zeta3_m_p2*1.30,  # mortality rate from ISR, male
    zeta2_f_p2 = zeta3_m_p2*2.40,  # mortality rate from ISR, female
    zeta4_m_p2 = zeta3_m_p2*1.62,  # mortality rate from ISS , male
    zeta4_f_p2 = zeta3_m_p2*2.23,  # mortality rate from ISS , female
    
    #Recovery rates from infection, including IMR, ISR, IMS and ISS due to treatment received [1/unit time] [%].
    nu1_m_p2 = (1-zeta1_m_p2), #recovery rates from IMR, males
    nu1_f_p2 = (1-zeta1_f_p2), #recovery rates from IMR, females
    nu2_m_p2 = (1-zeta2_m_p2), #recovery rates from ISR, males
    nu2_f_p2 = (1-zeta2_f_p2), #recovery rates from ISR, females
    nu3_m_p2 = (1-zeta3_m_p2), #recovery rates from IMS, males
    nu3_f_p2 = (1-zeta3_f_p2), #recovery rates from IMS, females
    nu4_m_p2 =  (1-zeta4_m_p2), #recovery rates from ISS, males
    nu4_f_p2 =  (1-zeta4_f_p2),#recovery rates from ISS, females
    
    #Constant background rate that captures transmission from non-human sources, horizontal transmission, or de novo emergence [1/unit time] [number].
    b_p2 = 0.01,
    
    # Percentage of people with resistant infections receiving inappropriate empirical antibiotic treatment [1/unit time] [%].
    phi_m_p2 = 0.3782,  # Placeholder value, adjust as needed
    phi_f_p2 = 0.3846,  # Placeholder value, adjust as needed
    
    #Factor of burden associated with inappropriate empirical antibiotic treatment and increased ICU admission among resistant infections [1/unit time] [%].
    pi_p2= 1.02,
    
    #Transmission parameter {update this correspondingly after calibrating it with real data}
    tau_p2= 0.2086551,
    
    #community-acquired infection upon hospital admission rate
    caIha_p2=0.007,
    
    #percentage of people tested
    test_p2=0.20, 
    HR_perc2=0.2,
    or_HR_scenar1_a=1.04,
    or_HR_scenarMen_a=2.27,
    ### ### ### ### ###
    #sensitivity chrom_1
    sens_chrom_a=0.826,
    #sensitivity chrom_1
    sens_chrom2_a=0.90, 
    #sensitivity chrom_1
    sens_pcr_a=1,
    #turnaround chrom_1
    turn_chrom_a=3,
    #turnaround chrom_1
    turn_chrom2_a=2,
    #turnaround pcr_1
    turn_pcr_a=1,
    #isolation contact precaution transmission reduction
    reduc_conpre_a=0.35,
    #efficiency decolonisation
    eff_decol_a= 0.26,
    #effect on self-infection decolonisation
    eff_decol_selfi_a=0.041,
    #Turnaround decolonisation program in days
    turnaround_decol_a=7,
    ##
    #costs wards
    c_general_ward= 50,
    c_intermediate_ward=92,
    c_icu_ward=218,
    c_decol_1pd=72.88,
    c_isolation=42.3,
    c_chrom=10.2,
    c_chrom2=13.6, 
    c_pcr=33,
    c_bc=16.9,
    #utilities
    u_healthy=0.92,
    u_icu=0.92-0.34,
    u_gw=0.64,
    u_recovICU=0.74
    
  )
  
  #####
  # Construct the function name from the model_list string
  model_func_name <- get(model_list[i])
  # Solve the ODE using the 'deSolve' package's ode function
  O_solution <- ode(y = state2, times = times, func = model_func_name, parms = parameters2, method = "rk4")
  O_solution2 <- as.data.frame(O_solution)
  O_solution <- as.data.frame(O_solution)
  # Process the ODE solution: sum all values reported in the first 350 rows for your variable of interest
  results_epi_ac[1, i] <- sum(O_solution[1:366, "CR_m2"]) + sum(O_solution[1:366, "CR_f2"]) + sum(O_solution[1:366, "IMR_m2"]) + sum(O_solution[1:366, "IMR_f2"])+ sum(O_solution[1:366, "ISR_m2"])+sum(O_solution[1:366, "ISR_f2"])
  results_epi_ac[2, i] <- sum(O_solution[1:366, "IMR_m2"]) + sum(O_solution[1:366, "IMR_f2"])+ sum(O_solution[1:366, "ISR_m2"])+sum(O_solution[1:366, "ISR_f2"])
  results_epi_ac[3, i] <- (O_solution[366, "DR_f2"]) + (O_solution[366, "DR_m2"]) 
  results_epi_ac[4, i] <- sum(O_solution[1:366, "IMR_m2"]) + sum(O_solution[1:366, "IMR_f2"])+ sum(O_solution[1:366, "ISR_m2"])+sum(O_solution[1:366, "ISR_f2"]) + sum(O_solution[1:366, "IMS_m2"]) + sum(O_solution[1:366, "IMS_f2"])+ sum(O_solution[1:366, "ISS_m2"])+sum(O_solution[1:366, "ISS_f2"])
  results_epi_ac[5, i] <- (O_solution[366, "DR_f2"]) + (O_solution[366, "DR_m2"]) + (O_solution[366, "DS_f2"]) + (O_solution[366, "DS_m2"]) 
  results_epi_ac[6, i] <- sum(O_solution[366, "new_admin"]) 
  results_epi_ac[7, i] <- sum(O_solution[1:366, "U_m2"]) +sum(O_solution[1:366, "U_f2"]) 
  results_epi_ac[8, i] <- sum(O_solution[1:366, "U_m2"]) +sum(O_solution[1:366, "U_f2"]) +sum(O_solution[1:366, "CR_m2"]) +sum(O_solution[1:366, "CR_f2"]) +sum(O_solution[1:366, "CS_m2"]) +sum(O_solution[1:366, "CS_f2"]) 
  results_epi_ac[9, i] <- O_solution[366, "RR_m2"] + O_solution[366, "RR_f2"]+O_solution[366, "RS_m2"] + O_solution[366, "RS_f2"]
  results_epi_ac[10, i] <- (O_solution[366, "discharge"]) 
  
  results_epi_CREprev[, i] <- O_solution2$CR_m2 + O_solution2$CR_f2 + O_solution2$IMR_m2 + O_solution2$IMR_f2 + O_solution2$ISR_m2 + O_solution2$ISR_f2
  results_epi_CREinfe[, i]  <- O_solution2$IMR_m2 + O_solution2$IMR_f2 + O_solution2$ISR_m2 + O_solution2$ISR_f2
  results_epi_CREdead[, i]  <- O_solution2$DR_m2 + O_solution2$DR_f2
  results_epi_CREinfe_all[, i]  <-  O_solution2$IMR_m2 + O_solution2$IMR_f2 + O_solution2$ISR_m2 + O_solution2$ISR_f2 + O_solution2$IMS_m2 + O_solution2$IMS_f2 + O_solution2$ISS_m2 + O_solution2$ISS_f2
  results_epi_CREdead_all[, i]  <- O_solution2$DR_m2 + O_solution2$DR_f2 + O_solution2$DS_m2 + O_solution2$DS_f2
  
  #Economics
  #Store econ results per strategy
  results_econ[1, i] <- O_solution[366, "cost"]
  results_econ[2, i] <- (results_epi_ac[8, i]*0.92)+(sum(O_solution[1:366, "ISR_m2"]))*0.58+(sum(O_solution[1:366, "ISR_f2"])*0.58)+((sum(O_solution[1:366, "IMR_m2"])+sum(O_solution[1:366, "IMR_f2"]))*0.64)+((sum(O_solution[1:366, "IMS_m2"])+sum(O_solution[1:366, "IMS_f2"]))*0.64)+(sum(O_solution[1:366, "ISS_m2"]))*0.58+(sum(O_solution[1:366, "ISS_f2"])*0.58)+ (O_solution[366, "RR_m2"])*0.92+(O_solution[366, "RR_f2"]*0.92)+ (O_solution[366, "RS_m2"])*0.92+(O_solution[366, "RS_f2"]*0.92)+results_epi_ac[10, i]*0.92
  results_econ[3, i] <- results_econ[2, i]/O_solution[366,"new_admin"] #check if usage is appropriate
}
#Compute ICER per strategy
results_icer <- matrix(nrow = 13, ncol = length(model_list))
results_icer[1,1] <- 0
results_icer[2,1] <- results_econ[1, 1]
results_icer[3,1] <- results_econ[2, 1]
results_icer[4,1] <- results_epi_ac[1, 1]
results_icer[5,1] <- results_epi_ac[2, 1]
results_icer[6,1] <- results_epi_ac[3, 1]
results_icer[7,1] <- results_epi_ac[4, 1]
results_icer[8,1] <- results_epi_ac[5, 1]
results_icer[9,1] <- results_epi_ac[9, 1]
results_icer[10,1] <- results_epi_ac[8, 1]
results_icer[11,1] <- results_epi_ac[10, 1]
results_icer[12,1] <- results_epi_ac[10, 1]+  results_epi_ac[9, 1]+results_epi_ac[8, 1]+results_epi_ac[4, 1]
results_icer[13,1] <- results_epi_ac[6, 1]
for (i in 2:length(model_list)) {
  results_icer[1, i] <- 0
  results_icer[2, i] <-  results_econ[1, i]
  results_icer[3, i] <-  results_econ[2, i]
  results_icer[4, i]<- results_epi_ac[1, i]
  results_icer[5, i]<- results_epi_ac[2, i]
  results_icer[6, i]<- results_epi_ac[3, i]
  results_icer[7, i]<- results_epi_ac[4, i]
  results_icer[8, i]<- results_epi_ac[5, i]
  results_icer[9, i]<- results_epi_ac[9, i]
  results_icer[10, i]<- results_epi_ac[8, i]
  results_icer[11, i]<- results_epi_ac[10, i]
  results_icer[12, i]<- results_epi_ac[10, i]+  results_epi_ac[9, i]+results_epi_ac[8, i]+results_epi_ac[4, i]
  results_icer[3, i] <- results_econ[2, i] + ifelse((results_icer[12, 1] - results_icer[12, i]) > 0, (results_icer[12, 1] - results_icer[12, i]) * 0.92, 0)
  results_icer[1, i] <- (results_icer[2, i]-results_icer[2, 1])/(results_icer[3, i]-results_icer[3, 1])
  results_icer[13,i] <- results_epi_ac[6, i]
  
}
results_epi_acC<- results_epi_ac
results_econC<- results_econ
results_icerC<- results_icer
transposed_icer <- t(results_icer)
transposed_icerdf <- as.data.frame(transposed_icer)
colnames(transposed_icerdf) <- c("ICER", "Costs", "QALYs","ARB colonisation","ARB infections", "ARB deaths", "Total infections","Total deaths","Total Recovered","Total U,CR,CS", "Total discharge","total population")  # Add more names as needed
rownames(transposed_icerdf) <- c("Do-nothing", "Strategy2", "Strategy3", "Strategy4", "Strategy5", "Strategy6", "Strategy7", "Strategy8", "Strategy9", "Strategy10","Strategy11", "Strategy12","Strategy13","Strategy14","Strategy15", "Strategy16")
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/AMR_BSI_BurdenChile/0_Article_AMR Transmission dynamics Chile & LMICs/0_BSIModelling/0_analysis/0_Figures_model")
write.csv(transposed_icerdf, "matrix_resCRE_lB.csv", row.names = TRUE)
#####
lower_boundCRE_inf<-(results_epi_ac[2, ]/results_epi_ac[6, ])*100
lower_boundCRE_col <- ((results_epi_ac[1, ]-results_epi_ac[2,])/results_epi_ac[6, ])*100
lower_boundCRE_death<- results_epi_ac[3,]



dif_CI_cre<- abs(infecCRE_a-Upper_boundCRE_inf)
dif_CI_mrsa<- abs(infecMRSA_a-lower_boundMRSA_inf)

dif_CI_creC<- abs(colonCRE_a-Upper_boundCRE_col)
dif_CI_mrsaC<- abs(colonMRSA_a-lower_boundMRSA_col)

dif_CI_creD<- abs(Dead_CRE_tot-Upper_boundCRE_death)
dif_CI_mrsaD<- abs(Dead_MRSA_tot-lower_boundMRSA_death)


####EDITING DATA MANAGEMENT:
# Removing the 7th to 12th elements from the list
colonMRSA_af <- colonMRSA_a[-c(7, 8, 9, 10, 11, 12)]
dif_CI_cref <-dif_CI_cre[-c(7, 8, 9, 10, 11, 12)]
dif_CI_mrsaf<- dif_CI_mrsa[-c(7, 8, 9, 10, 11, 12)]
infecMRSA_af<-infecMRSA_a[-c(7, 8, 9, 10, 11, 12)]
DeadMRSA_af<- DeadMRSA_a[-c(7, 8, 9, 10, 11, 12)]
colonCRE_af<- colonCRE_a[-c(7, 8, 9, 10, 11, 12)]
DeadCRE_af<- DeadCRE_a[-c(7, 8, 9, 10, 11, 12)]
dif_CI_creC <- dif_CI_creC[-c(7, 8, 9, 10, 11, 12)]
dif_CI_mrsaC <- dif_CI_mrsaC[-c(7, 8, 9, 10, 11, 12)]
dif_CI_creD <-dif_CI_creD[-c(7, 8, 9, 10, 11, 12)]
dif_CI_mrsaD <-dif_CI_mrsaD[-c(7, 8, 9, 10, 11, 12)]
  
tiff("bar_chart_1_ci.tiff", width = 14, height = 10, units = 'in', res = 800)
par(mfrow = c(2, 1), mar = c(3.5, 4, 2, 1), mgp = c(2, 0.5, 0))
# First barplot for MRSA-colonised

bp1<-barplot(colonMRSA_af, col = lancet_colors, main = "", xlab = "", ylab = "MRSA-colonised population/100 new admissions", ylim = c(0, 130), yaxt = "n")
axis(side = 2, at = seq(0, 130, by = 10), las = 1)
mtext("(A) MRSA colonised individuals, not infected", side = 3, line = 1, at = -1, adj = 0, font = 2)
arrows(bp1, colonMRSA_af+dif_CI_mrsaC, bp1, colonMRSA_af-dif_CI_mrsaC, angle = 90, code = 3, length = 0.05)



# Second barplot for MRSA-infected
bp2<-barplot(infecMRSA_af, col = lancet_colors, main = "", xlab = "Strategy", ylab = "MRSA-infected population/100 new admissions", ylim = c(0, 28), yaxt = "n", names.arg = c("S0","S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9"))
axis(side = 2, at = seq(0, 28, by = 2), las = 1)
mtext("(B) MRSA infections", side = 3, line = 1, at = -1, adj = 0, font = 2)
arrows(bp2, infecMRSA_af+dif_CI_mrsaf, bp1,  infecMRSA_af-dif_CI_mrsaf, angle = 90, code = 3, length = 0.05)
# Adjusting the legend position and appearance
legend("bottom", inset=c(0, 0.94),  # Adjust for precise positioning between the plots
       legend=c("S0: Do-nothing", expression("S1: T+D"["agar1"]), expression("S2: T+D"["agar2"]), expression("S3: T+D"["pcr"]), expression("S4: T+I"["agar1"]), expression("S5: T+I"["agar2"]), expression("S6: T+I"["pcr"]), expression("S7: Pre-emptive I"["all"]), expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"])), 
       title="", cex=1.2, ncol=4, bty='n', xpd=NA,
       y.intersp=0.9, x.intersp=0.3, pch=15, pt.cex=0,  col=lancet_colors, fill = lancet_colors)  # Adjusted for more compact spacing and symbol size
# Please adjust 'cex' for size of text and symbols in legend
theme_lancet() 
#barplot(DeadMRSA_a, col = bar_colors, main = "MRSA-deaths/100 new admissions", xlab = "Strategy", ylab = "Value")
dev.off()

colonCRE_af <- colonCRE_a[-c(7, 8, 9, 10, 11, 12)]
Upper_boundCRE_inff <-Upper_boundCRE_inf[-c(7, 8, 9, 10, 11, 12)]
lower_boundCRE_inff<- lower_boundCRE_inf[-c(7, 8, 9, 10, 11, 12)]
Upper_boundCRE_colf<-  Upper_boundCRE_col[-c(7, 8, 9, 10, 11, 12)]
lower_boundCRE_colf <-lower_boundCRE_col[-c(7, 8, 9, 10, 11, 12)]
infecCRE_af<-infecCRE_a[-c(7, 8, 9, 10, 11, 12)]
tiff("bar_chart_2_ci.tiff", width = 14, height = 10, units = 'in', res = 800)
par(mfrow = c(2, 1), mar = c(3.5, 4, 2, 1), mgp = c(2, 0.5, 0))
# First barplot for CRE-colonised
bp1<-barplot(colonCRE_af, col = lancet_colors, main = "", xlab = "", ylab = "CRE-colonised population/100 new admissions", ylim = c(0, 170), yaxt = "n")
axis(side = 2, at = seq(0, 170, by = 10), las = 1)
mtext("(A) CRE colonised individuals, not infected", side = 3, line = 1, at = -1, adj = 0, font = 2)
arrows(bp1, colonCRE_af+dif_CI_creC, bp1, colonCRE_af-dif_CI_creC, angle = 90, code = 3, length = 0.05)

# Second barplot for CRE-infected
bp2<-barplot(infecCRE_af, col = lancet_colors, main = "", xlab = "Strategy", ylab = "CRE-infected population/100 new admissions", ylim = c(0, 40), yaxt = "n", names.arg = c("S0","S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9"))
axis(side = 2, at = seq(0, 40, by = 2), las = 1)
mtext("(B) CRE infections", side = 3, line = 1, at = -1, adj = 0, font = 2)
arrows(bp2, infecCRE_af+dif_CI_cref, bp1,  infecCRE_af-dif_CI_cref, angle = 90, code = 3, length = 0.05)
# Adjusting the legend position and appearance
legend("bottom", inset=c(0, 0.94),  # Adjust for precise positioning between the plots
       legend=c("S0: Do-nothing", expression("S1: T+D"["agar1"]), expression("S2: T+D"["agar2"]), expression("S3: T+D"["pcr"]), expression("S4: T+I"["agar1"]), expression("S5: T+I"["agar2"]), expression("S6: T+I"["pcr"]), expression("S7: Pre-emptive I"["all"]), expression("S8: Pre-emptive I"["men"]), expression("S9: Pre-emptive I"["women"])), 
       title="", cex=1.2, ncol=4, bty='n', xpd=NA,
       y.intersp=0.9, x.intersp=0.3, pch=15, pt.cex=0,  col=lancet_colors, fill = lancet_colors)  # Adjusted for more compact spacing and symbol size
# Please adjust 'cex' for size of text and symbols in legend
theme_lancet() 
#barplot(DeadMRSA_a, col = bar_colors, main = "MRSA-deaths/100 new admissions", xlab = "Strategy", ylab = "Value")
dev.off()



model_listX <- c(
  "S0: Do-nothing",
  "S1: T+D, agar1",
  "S2: T+D, agar2",
  "S3: T+D, pcr",
  "S4: T+I, agar1",
  "S5: T+I, agar2",
  "S6: T+I, pcr",
  "S7: Pre-emptive I, all",
  "S8: Pre-emptive I, men",
  "S9: Pre-emptive I, women"
)
strategy_matrix_NB <- cbind(model_listX, lower_boundCRE_col, Upper_boundCRE_col, lower_boundMRSA_col, Upper_boundMRSA_col)
# Convert to data frame for more flexible handling (optional but recommended)
strategy_df_nb <- as.data.frame(strategy_matrix_NB)
write.csv(strategy_df_nb, "strategy_bounds_nb.csv", row.names = FALSE)




## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
