library(dada2)
library(phyloseq)
library(tidyverse)
library(ggfittext)
library(seqinr)
library(cowplot)
library(exactRankTests)
library(rlist)
library(foreach)
library(RColorBrewer)
library(data.table)
library(grid)
library(scales)

setwd("../../")
source("analysis/generateMixturesDataframe.R")
plotdir <- "analysis/individualOTUPlots"

# Import family coding.
familyCoding <- fread("config/familyCodingFull.txt")

# Example ASV trajectories. -----------------------------------------------

p_OTUAbundancesExampleFlat <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XDA-XDB" & Family=="Bacteroidaceae" & OTUnum==10) %>% 
  mutate(lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2"),
         plot = "Resident\nBacteroidaceae-10 (C1/C2)") %>% 
  mutate(lineGroup = fct_relevel(lineGroup, c("parent","theoretical-1","theoretical-2","actual-1","actual-2","actual-3"))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Experimental 1",
                              "actual-2"="Experimental 2",
                              "actual-3"="Experimental 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~plot) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none")
p_OTUAbundancesExampleFlat
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleFlat.png"),
          p_OTUAbundancesExampleFlat, base_width=1.6, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/OTU_abundance_ExampleFlat.pdf"),
          p_OTUAbundancesExampleFlat, base_width=1.7, base_height=1.65)

trajectoryLegend <- get_legend(p_OTUAbundancesExampleFlat)
save_plot(paste0(plotdirParentDiff, "/OTUtrajectoriesLegend.pdf"), 
          trajectoryLegend, base_width=0.8, base_height=0.7)

p_OTUAbundancesExampleNoisy3 <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XFA-XBA" & Family=="Rikenellaceae" & OTUnum==4) %>% 
  mutate(plot="Noisy\nRikenellaceae-4 (D1/A1)") %>%
  mutate(lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2")) %>% 
  filter(!(passage==5 & combo=="XFA-XBA" & ratio=="100:1" & line=="actual-1")) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~plot)
p_OTUAbundancesExampleNoisy3
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleNoisy.png"),
          p_OTUAbundancesExampleNoisy3, base_width=1.6, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/OTU_abundance_ExampleNoisy.pdf"),
          p_OTUAbundancesExampleNoisy3, base_width=1.7, base_height=1.65)

p_OTUAbundancesExampleUndercolonizing <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XCA-XDA" & Family=="Lachnospiraceae" & OTUnum==144) %>% 
  mutate(plot = "Weak colonizer\nLachnospiraceae-144 (B1/C1)") %>% 
  mutate(lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2")) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~plot)
p_OTUAbundancesExampleUndercolonizing
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleUndercolonizing.png"),
          p_OTUAbundancesExampleUndercolonizing, base_width=1.6, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/OTU_abundance_ExampleUndercolonizing.pdf"),
          p_OTUAbundancesExampleUndercolonizing, base_width=1.7, base_height=1.65)

p_OTUAbundancesExampleOvercolonizing <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XFA-XBA" & Family=="Lachnospiraceae" & OTUnum==115) %>%
  # Filter out XFA-XBA poorly sequenced replicate/ratio.
  filter(!(ratio=="100:1" & line=="actual-1")) %>% 
  #filter(line %in% c("theoretical-1","theoretical-2")) %>% 
  mutate(plot = "Strong colonizer\nLachnospiraceae-115 (D1/A1)") %>% 
  mutate(lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2")) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~plot)
p_OTUAbundancesExampleOvercolonizing
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleOvercolonizingTheoretical.png"),
          p_OTUAbundancesExampleOvercolonizingTheoretical, base_width=1.6, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/OTU_abundance_ExampleOvercolonizing.pdf"),
          p_OTUAbundancesExampleOvercolonizing, base_width=1.7, base_height=1.65)

p_OTUAbundancesExampleDoseDependent <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XCA-XDA" & Family=="Bacteroidaceae" & OTUnum==9) %>%
  mutate(plot = "Dose-dependent\nBacteroidaceae-9 (B1/C1)") %>% 
  mutate(lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2")) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~plot)
p_OTUAbundancesExampleDoseDependent
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleDoseDependent.png"),
          p_OTUAbundancesExampleDoseDependent, base_width=1.6, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/OTU_abundance_ExampleDoseDependent.pdf"),
          p_OTUAbundancesExampleDoseDependent, base_width=1.7, base_height=1.65)

# Neutral examples. -------------------------------------------------------

# Combined example with Streptococcaceae and Enterococcaceae.
p_OTUAbundancesExampleNeutralStrepEnteC <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & ((combo=="XFA-XFB" & Family=="Streptococcaceae" & OTUnum %in% c(7,17)) | 
                         (combo=="XFA-XFB" & Family=="Enterococcaceae" & OTUnum %in% c(2,3)))) %>% 
  left_join(familyCoding %>% rename(Family=family)) %>% 
  mutate(OTUID = paste0(code,"-",OTUnum),
         lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2"),
         lineGroupStrain = case_when(lineGroup == "parent" ~ "parent", 
                                     lineGroup %in% c("theoretical-1","theoretical-2") ~ lineGroup,
                                     lineGroup %in% c("actual-1","actual-2","actual-3") ~ paste0(lineGroup, "_", OTUID))) %>% 
  ungroup() %>% 
  mutate(plot = case_when(OTUID=="Strep-7" ~ "L. garvieae",
                          OTUID=="Strep-17" ~ "L. lactis",
                          OTUID=="EnteC-2" ~ "E. faecalis",
                          OTUID=="EnteC-3" ~ "E. casseliflavus"),
         plot = fct_relevel(plot, c("E. faecalis","E. casseliflavus",
                                    "L. garvieae","L. lactis"))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, group=interaction(lineGroupStrain, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroupStrain, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical"="Theoretical",
                              "actual-1_EnteC-2"="Actual 1 EnteC2",
                              "actual-2_EnteC-2"="Actual 2 EnteC2",
                              "actual-3_EnteC-2"="Actual 3 EnteC2",
                              "actual-1_EnteC-3"="Actual 1 EnteC3",
                              "actual-2_EnteC-3"="Actual 2 EnteC3",
                              "actual-3_EnteC-3"="Actual 3 EnteC3",
                              "actual-1_Strep-7"="Actual 1 Strep7",
                              "actual-2_Strep-7"="Actual 2 Strep7",
                              "actual-3_Strep-7"="Actual 3 Strep7",
                              "actual-1_Strep-17"="Actual 1 Strep17",
                              "actual-2_Strep-17"="Actual 2 Strep17",
                              "actual-3_Strep-17"="Actual 3 Strep17"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080",
                              "theoretical-2"="#707070",
                              "actual-1_EnteC-2"="#D3BC4A", 
                              "actual-2_EnteC-2"="#C1A92F",
                              "actual-3_EnteC-2"="#A48F28",
                              "actual-1_EnteC-3"="#D38D4A",
                              "actual-2_EnteC-3"="#BF732E",
                              "actual-3_EnteC-3"="#A46428",
                              "actual-1_Strep-7"="#7EAAB4",
                              "actual-2_Strep-7"="#6398A4",
                              "actual-3_Strep-7"="#53838D",
                              "actual-1_Strep-17"="#4044B5",
                              "actual-2_Strep-17"="#343795",
                              "actual-3_Strep-17"="#2A2D79"), guide="none") +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        strip.text=element_text(face="bold.italic")) +
  theme(legend.position="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~plot, nrow=1)
p_OTUAbundancesExampleNeutralStrepEnteC
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleNeutralStrepEnteC-p5.pdf"),
          p_OTUAbundancesExampleNeutralStrepEnteC, base_width=4.45, base_height=1.6)

# B. fragilis strain in XFA/XFB mixtures.
Bacte0126 <- "TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGACTGGTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGTCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGATGCTCGAAAGTGTGGGTATCAAACAGG"
p_OTUAbundancesExampleBfrag <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XFA-XFB" & OTU==Bacte0126) %>%  
  mutate(plot = "B. fragilis",
         lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2")) %>% 
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080",
                              "theoretical-2"="#707070",
                              "actual-1"="#943874", 
                              "actual-2"="#762D5D",
                              "actual-3"="#592245"), guide="none") +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        strip.text=element_text(face="bold.italic")) +
  theme(legend.position="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~plot)
p_OTUAbundancesExampleBfrag
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleBfrag.pdf"),
          p_OTUAbundancesExampleBfrag, base_width=1.7, base_height=1.65)
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleBfrag.png"),
          p_OTUAbundancesExampleBfrag, base_width=1.6, base_height=1.6)

# P. goldsteinii strain in XFA/XFB mixtures.
Tanne0007 <- "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGTGGTTAATTAAGTCAGCGGTGAAAGTTTGTGGCTCAACCATAAAATTGCCGTTGAAACTGGTTGACTTGAGTATATTTGAGGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCTTACTAAACTATAACTGACACTGAAGCACGAAAGCGTGGGGATCAAACAGG"
p_OTUAbundancesExamplePgold <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XFA-XFB" & OTU==Tanne0007) %>% 
  mutate(plot = "P. goldsteinii",
         lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2")) %>% 
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080",
                              "theoretical-2"="#707070",
                              "actual-1"="#DA2F35", 
                              "actual-2"="#B62025",
                              "actual-3"="#9C1C20"), guide="none") +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        strip.text=element_text(face="bold.italic")) +
  theme(legend.position="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~plot)
p_OTUAbundancesExamplePgold
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExamplePgold.pdf"),
          p_OTUAbundancesExamplePgold, base_width=1.7, base_height=1.65)
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExamplePgold.png"),
          p_OTUAbundancesExamplePgold, base_width=1.6, base_height=1.6)

# Example passage 3 to 5. -------------------------------------------------

p_OTUAbundancesExampleTanne <- theoretical_mixtures_df_categorized %>%
  filter(combo=="XDA-XFA" & Family=="Tannerellaceae" & OTUnum==27) %>% 
  #left_join(familyCoding %>% rename(Family=family)) %>% 
  mutate(#OTUID = paste0(code,"-",OTUnum),
         lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2")) %>% 
  mutate(plot = ifelse(passage==3,"Passage 3\nDose dependent", 
                       "Passage 5\nNoisy")) %>% 
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~plot)
p_OTUAbundancesExampleTanne
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleTanne.pdf"),
          p_OTUAbundancesExampleTanne, base_width=2.5, base_height=1.7)
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleTanne.png"),
          p_OTUAbundancesExampleTanne, base_width=3.2, base_height=2)

# Plot all ASVs from one mixture. -----------------------------------------

p_XDAXDBDD <- theoretical_mixtures_df_categorized %>% 
  filter(passage==5 & combo=="XDA-XDB" & category=="dose-dependent") %>% 
  #left_join(familyCoding %>% rename(Family=family)) %>% 
  mutate(lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2"),
         OTUID = paste0(Family, "-", OTUnum)) %>% 
  mutate(lineGroup = fct_relevel(lineGroup, c("parent","theoretical-1","theoretical-2","actual-1","actual-2","actual-3"))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  theme(strip.background=element_rect(fill="#811517"),
        strip.text=element_text(color="white")) +
  facet_wrap(~OTUID, nrow=1)
p_XDAXDBDD
save_plot(paste0(plotdirParentDiff, "/XDAXDB-DD.pdf"), p_XDAXDBDD, base_width=5.5, base_height=1.75)

p_XDAXDBOC <- theoretical_mixtures_df_categorized %>% 
  filter(passage==5 & combo=="XDA-XDB" & category=="overcolonizing") %>% 
  #left_join(familyCoding %>% rename(Family=family)) %>% 
  mutate(lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2"),
         OTUID = paste0(Family, "-", OTUnum)) %>%  
  mutate(lineGroup = fct_relevel(lineGroup, c("parent","theoretical-1","theoretical-2","actual-1","actual-2","actual-3"))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  theme(strip.background=element_rect(fill="#D4A129")) +
  facet_wrap(~OTUID, nrow=1)
p_XDAXDBOC
save_plot(paste0(plotdirParentDiff, "/XDAXDB-OC.pdf"), p_XDAXDBOC, base_width=4.5, base_height=1.75)

p_XDAXDBUC <- theoretical_mixtures_df_categorized %>% 
  filter(passage==5 & combo=="XDA-XDB" & category=="undercolonizing") %>% 
  #left_join(familyCoding %>% rename(Family=family)) %>% 
  mutate(lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2"),
         OTUID = paste0(Family, "-", OTUnum)) %>%  
  mutate(lineGroup = fct_relevel(lineGroup, c("parent","theoretical-1","theoretical-2","actual-1","actual-2","actual-3"))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  theme(strip.background=element_rect(fill="#86B7E2")) +
  facet_wrap(~OTUID, nrow=1)
p_XDAXDBUC
save_plot(paste0(plotdirParentDiff, "/XDAXDB-UC.pdf"), p_XDAXDBUC, base_width=6.5, base_height=1.75)

p_XDAXDBFlat <- theoretical_mixtures_df_categorized %>% 
  filter(passage==5 & combo=="XDA-XDB" & category=="flat") %>% 
  #left_join(familyCoding %>% rename(Family=family)) %>% 
  mutate(lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2"),
         OTUID = paste0(Family, "-", OTUnum)) %>%  
  mutate(lineGroup = fct_relevel(lineGroup, c("parent","theoretical-1","theoretical-2","actual-1","actual-2","actual-3"))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  theme(strip.background=element_rect(fill="#616667"),
        strip.text=element_text(color="white")) +
  facet_wrap(~OTUID, nrow=4, scales="free_x")
p_XDAXDBFlat
save_plot(paste0(plotdirParentDiff, "/XDAXDB-flat.pdf"), p_XDAXDBFlat, base_width=6.4, base_height=5)

p_XDAXDBNoisy <- theoretical_mixtures_df_categorized %>% 
  filter(passage==5 & combo=="XDA-XDB" & category=="noisy") %>% 
  #left_join(familyCoding %>% rename(Family=family)) %>% 
  mutate(lineGroup = case_when(group=="control" ~ "parent",
                               ratio %in% ratios & line=="theoretical-1" ~ "theoretical-1",
                               ratio %in% ratios & line=="theoretical-2" ~ "theoretical-2",
                               ratio %in% ratios & line=="actual-1" ~ "actual-1",
                               ratio %in% ratios & line=="actual-2" ~ "actual-2",
                               ratio %in% ratios & line=="actual-3" ~ "actual-3"),
         parentGroup = case_when(group=="mixture" ~ "mixture",
                                 ratio=="1:0" & inoculationReplicate==1 ~ "leftParent1",
                                 ratio=="1:0" & inoculationReplicate==2 ~ "leftParent2",
                                 ratio=="0:1" & inoculationReplicate==1 ~ "rightParent1",
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2"),
         OTUID = paste0(Family, "-", OTUnum)) %>%  
  mutate(lineGroup = fct_relevel(lineGroup, c("parent","theoretical-1","theoretical-2","actual-1","actual-2","actual-3"))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=19,"openCircle"=1), guide="none") +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  theme(strip.background=element_rect(fill="#999999"),
        strip.text=element_text(color="white")) +
  facet_wrap(~OTUID, nrow=4, scales="free_x")
p_XDAXDBNoisy
save_plot(paste0(plotdirParentDiff, "/XDAXDB-noisy.pdf"), p_XDAXDBNoisy, base_width=5.2, base_height=5)