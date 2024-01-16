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

# Plot all OTU trajectories in all categories but lowAbundance, per combo and passage. --------

#plotdirParentDiff <- paste0(plotdir, "/trajectoryPlots-parentDiff")
plotdirParentDiff <- paste0(plotdir, "/trajectoryPlots-parentDiff-noTheoretical")

categories <- c("dose-dependent","overcolonizing","undercolonizing","noisy","flat")

# Use this function to plot ASVs in each category separately.
plotComboOTUAbundancesHighAbundanceOneParent <- function(currentCombo, passageNum, currCategory) {
  df <- theoretical_mixtures_df_categorized %>%
    filter(category==currCategory) %>% 
    filter(combo==currentCombo & passage==passageNum) %>% 
    mutate(FamilyOTU=paste(Family, OTUnum, sep="-")) %>% 
    # Filter out XFA-XBA poorly sequenced replicate/ratio.
    filter(!(passage==5 & combo=="XFA-XBA" & ratio=="100:1" & line=="actual-1"))
  if (nrow(df)>0) {
    p_OTUAbundancesCombo <- df %>%
      ggplot() +
      geom_line(aes(x=ratio, y=logAbundance, color=line, group=line), alpha=0.6) +
      geom_point(aes(x=ratio, y=logAbundance, color=line, shape=pointType), alpha=0.8) +
      geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
      xlab("Ratio") +
      ylab("Relative abundance (log)") +
      ylim(c(-3, 0)) +
      ggtitle(paste0(currentCombo, "-p", passageNum)) +
      scale_color_manual(name="Mixture type",
                         labels=c("theoretical-1"="Theoretical 1",
                                  "theoretical-2"="Theoretical 2",
                                  "actual-1"="Actual 1",
                                  "actual-2"="Actual 2",
                                  "actual-3"="Actual 3"), 
                         values=c("theoretical-1"="#313695", 
                                  "theoretical-2"="#4575B4", 
                                  "actual-1"="#D73027", 
                                  "actual-2"="#F46D43",
                                  "actual-3"="#FDAE61")) +
      scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
      DEFAULTS.THEME_PRES +
      facet_wrap(~FamilyOTU, scales="free") +
      theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
            strip.text.y=element_text(angle=0))
    plotName <- paste("OTU_abundances", currentCombo, sep="_")
    save_plot(paste0(plotdirParentDiff, "/", plotName, "-P", passageNum, "_", currCategory, ".png"),
              p_OTUAbundancesCombo, nrow=4, ncol=4, limitsize=FALSE)
  }
}

foreach(x=categories) %do% {
  sapply(combos, FUN=plotComboOTUAbundancesHighAbundanceOneParent, 3, x)
}

# OTU abundance scatterplots ----------------------------------------------

plotdirScatter <- paste0(plotdir, "/scatterplots")

# Create abundance dataframe with theoretical mean in one column and actual mean in another.
ASV_sum_abundances <- full_join(
  theoretical_mixtures_df_statsFiltered %>% 
    filter(mixtureType=="theoretical" & ratio %in% ratios) %>% 
    group_by(OTU, passage, combo, ratio) %>% 
    slice(1) %>% 
    select(OTU, OTUnum, Family, passage, combo, ratio, avgAbundance, sdAbundance) %>% 
    rename(t_avgAbundance=avgAbundance, t_sdAbundance=sdAbundance),
  theoretical_mixtures_df_statsFiltered %>% 
    filter(mixtureType=="actual") %>% 
    group_by(OTU, passage, combo, ratio) %>% 
    slice(1) %>% 
    select(OTU, passage, combo, ratio, avgAbundance, sdAbundance) %>% 
    rename(a_avgAbundance=avgAbundance, a_sdAbundance=sdAbundance),
  by=c("OTU", "passage", "combo", "ratio")) %>% 
  mutate(a_avgAbundance=replace(a_avgAbundance, is.na(a_avgAbundance), -4))

# Helper function to plot OTU abundance scatterplots for a smaller selection of OTUs.
plotOTUScatterplotSubset <- function(plotNumber, df, family, passageNum) {
  p_OTUAbundancesFamilyScatterplot <- df %>%
    filter(OTUnum==plotNumber) %>%
    ggplot() +
    geom_point(aes(x=t_avgAbundance, y=a_avgAbundance, color=ratio)) +
    geom_errorbarh(aes(y=a_avgAbundance, color=ratio, height=0,
                       xmin=t_avgAbundance-t_sdAbundance/2, xmax=t_avgAbundance+t_sdAbundance/2)) +
    geom_errorbar(aes(x=t_avgAbundance, color=ratio,
                      ymin=a_avgAbundance-a_sdAbundance/2, ymax=a_avgAbundance+a_sdAbundance/2)) +
    scale_color_brewer(palette="Reds", name="Ratio") +
    geom_abline(aes(slope=1, intercept=0)) +
    geom_hline(aes(yintercept=-4), linetype="dotted") +
    geom_vline(aes(xintercept=-4), linetype="dotted") +
    ylab("Average log actual relative abundance") +
    xlab("Average log theoretical relative abundance") +
    ylim(c(-4,0)) +
    xlim(c(-4,0)) +
    ggtitle(paste0(family, "-", plotNumber, "-p", passageNum)) +
    DEFAULTS.THEME_PRES +
    facet_wrap(~combo, scales="free")
  plotName <- paste("OTU_scatterplot", family, plotNumber, sep="_")
  save_plot(paste0(plotdirScatter, "/", plotName, "-P", passageNum, ".png"),
            p_OTUAbundancesFamilyScatterplot, nrow=1.5, ncol=2)
}

# Function to plot OTU scatterplots for each family and combo, faceted by OTU.
plotFamilyOTUScatterplots <- function(family, passageNum) {
  df <- ASV_sum_abundances %>% 
    filter(Family==family & passage==passageNum)
    # Annotate OTUs for plotting in separate plots.
  plotNums <- unique(df$OTUnum)
  sapply(plotNums, FUN=plotOTUScatterplotSubset, df, family, passageNum)
}

# Plot individual abundance scatterplots for all OTUs.
sapply(familyList, FUN=plotFamilyOTUScatterplots, 3)
sapply(familyList, FUN=plotFamilyOTUScatterplots, 5)


# Plot specific OTU figures for presentation. -----------------------------

p_OTUAbundances_Hafniaceae1P3 <- theoretical_mixtures_df_stats %>%
  filter(Family=="Hafniaceae" & passage==3) %>%
  # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
  select(OTU, OTUnum, ratio, combo, line, logAbundance) %>%
  group_by(combo, line) %>%
  complete(nesting(OTU, OTUnum), 
           ratio=ratiosFull,
           fill=list(logAbundance=NA)) %>%
  ungroup() %>%
  mutate(ratio=fct_relevel(ratio, ratiosFull)) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=line, group=line)) +
  geom_point(aes(x=ratio, y=logAbundance, color=line), alpha=0.5) +
  geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
  xlab("Ratio") +
  ylab("Relative abundance (log)") +
  ylim(c(-4.1, 0)) +
  ggtitle("Hafniaceae-1-p3") +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("theoretical-1"="chartreuse2", 
                              "theoretical-2"="green1", 
                              "actual-1"="slateblue4", 
                              "actual-2"="steelblue4",
                              "actual-3"="royalblue2")) +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
        strip.text.y=element_text(angle=0))
save_plot(paste0(plotdir, "/OTU_abundances_Hafniaceae-1_P3.png"),
          p_OTUAbundances_Hafniaceae1P3, nrow=1.5, ncol=2)

p_OTUAbundances_Hafniaceae1P5 <- theoretical_mixtures_df_stats %>%
  filter(Family=="Hafniaceae" & passage==5 & combo %in% c("XFA-XFB", "XDA-XFA", "XFA-XBA")) %>%
  # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
  select(OTU, OTUnum, ratio, combo, line, logAbundance) %>%
  group_by(combo, line) %>%
  complete(nesting(OTU, OTUnum), 
           ratio=ratiosFull,
           fill=list(logAbundance=NA)) %>%
  ungroup() %>%
  mutate(ratio=fct_relevel(ratio, ratiosFull)) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=line, group=line)) +
  geom_point(aes(x=ratio, y=logAbundance, color=line), alpha=0.5) +
  geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
  xlab("Ratio") +
  ylab("Relative abundance (log)") +
  ylim(c(-4.1, 0)) +
  ggtitle("Hafniaceae-1-p5") +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("theoretical-1"="chartreuse2", 
                              "theoretical-2"="green1", 
                              "actual-1"="slateblue4", 
                              "actual-2"="steelblue4",
                              "actual-3"="royalblue2")) +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
        strip.text.y=element_text(angle=0))
save_plot(paste0(plotdir, "/OTU_abundances_Hafniaceae-1_P5.png"),
          p_OTUAbundances_Hafniaceae1P5, nrow=1.5, ncol=2)

p_OTUAbundances_Oscillospiraceae7P3 <- theoretical_mixtures_df_stats %>%
  filter(Family=="Oscillospiraceae" & passage==3 & OTUnum==7) %>%
  # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
  select(OTU, OTUnum, ratio, combo, line, logAbundance) %>%
  group_by(combo, line) %>%
  complete(nesting(OTU, OTUnum), 
           ratio=ratiosFull,
           fill=list(logAbundance=NA)) %>%
  ungroup() %>%
  mutate(ratio=fct_relevel(ratio, ratiosFull)) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=line, group=line)) +
  geom_point(aes(x=ratio, y=logAbundance, color=line), alpha=0.5) +
  geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
  xlab("Ratio") +
  ylab("Relative abundance (log)") +
  ylim(c(-4.1, 0)) +
  ggtitle("Oscillospiraceae-7-p3") +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("theoretical-1"="chartreuse2", 
                              "theoretical-2"="green1", 
                              "actual-1"="slateblue4", 
                              "actual-2"="steelblue4",
                              "actual-3"="royalblue2")) +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
        strip.text.y=element_text(angle=0))
save_plot(paste0(plotdir, "/OTU_abundances_Oscillospiraceae-7_P3.png"),
          p_OTUAbundances_Oscillospiraceae7P3, nrow=1.5, ncol=2)

p_OTUAbundances_Oscillospiraceae7P5 <- theoretical_mixtures_df_stats %>%
  filter(Family=="Oscillospiraceae" & passage==5 & OTUnum==7 & combo %in% c("XBA-XBB", "XBA-XCA", "XFA-XBA")) %>%
  # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
  select(OTU, OTUnum, ratio, combo, line, logAbundance) %>%
  group_by(combo, line) %>%
  complete(nesting(OTU, OTUnum), 
           ratio=ratiosFull,
           fill=list(logAbundance=NA)) %>%
  ungroup() %>%
  mutate(ratio=fct_relevel(ratio, ratiosFull)) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=line, group=line)) +
  geom_point(aes(x=ratio, y=logAbundance, color=line), alpha=0.5) +
  geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
  xlab("Ratio") +
  ylab("Relative abundance (log)") +
  ylim(c(-4.1, 0)) +
  ggtitle("Oscillospiraceae-7-p5") +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("theoretical-1"="chartreuse2", 
                              "theoretical-2"="green1", 
                              "actual-1"="slateblue4", 
                              "actual-2"="steelblue4",
                              "actual-3"="royalblue2")) +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
        strip.text.y=element_text(angle=0))
save_plot(paste0(plotdir, "/OTU_abundances_Oscillospiraceae-7_P5.png"),
          p_OTUAbundances_Oscillospiraceae7P5, nrow=1.5, ncol=2)

p_OTUAbundances_Bacteroidaceae7P3_theoretical <- theoretical_mixtures_df_stats %>%
  filter(Family=="Bacteroidaceae" & passage==3 & OTUnum==7 & mixtureType=="theoretical") %>%
  # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
  select(OTU, OTUnum, ratio, combo, line, logAbundance) %>%
  group_by(combo, line) %>%
  complete(nesting(OTU, OTUnum), 
           ratio=ratiosFull,
           fill=list(logAbundance=NA)) %>%
  ungroup() %>%
  mutate(ratio=fct_relevel(ratio, ratiosFull)) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=line, group=line)) +
  geom_point(aes(x=ratio, y=logAbundance, color=line), alpha=0.5) +
  geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
  xlab("Ratio") +
  ylab("Relative abundance (log)") +
  ylim(c(-4.1, 0)) +
  ggtitle("Bacteroidaceae-7-p3") +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("theoretical-1"="chartreuse2", 
                              "theoretical-2"="green1", 
                              "actual-1"="slateblue4", 
                              "actual-2"="steelblue4",
                              "actual-3"="royalblue2")) +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
        strip.text.y=element_text(angle=0), legend.position="none")
save_plot(paste0(plotdir, "/OTU_abundances_Bacteroidaceae-7_P3_theoretical.png"),
          p_OTUAbundances_Bacteroidaceae7P3_theoretical, nrow=1.5, ncol=2)


# New filtering. ----------------------------------------------------------

# Set up trajectory plot directories.
plotdir104 <- paste0(plotdir, "/trajectoryPlots-10-4")
plotdir103 <- paste0(plotdir, "/trajectoryPlots-10-3")
plotdirOneParent <- paste0(plotdir, "/trajectoryPlots-oneParent")
plotdirHighAbundance <- paste0(plotdir, "/trajectoryPlots-highAbundance")

# Set up scatter plot directories.
plotdir104Scatter <- paste0(plotdir, "/scatterPlots-10-4")
plotdir103Scatter <- paste0(plotdir, "/scatterPlots-10-3")
plotdirOneParentScatter <- paste0(plotdir, "/scatterPlots-oneParent")
plotdirHighAbundanceScatter <- paste0(plotdir, "/scatterPlots-highAbundance")

# Function to plot OTU abundances for each family and combo, filtered to 10^-4.
plotFamilyOTUAbundances104 <- function(family, passageNum) {
  df <- theoretical_mixtures_df %>% 
    filter(filterRelAb4==TRUE & Family==family & passage==passageNum) %>%
    # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
    select(OTU, OTUnum, ratio, combo, line, logAbundance) %>%
    group_by(combo, line) %>%
    complete(nesting(OTU, OTUnum), 
             ratio=ratiosFull,
             fill=list(logAbundance=NA)) %>%
    ungroup() %>%
    mutate(ratio=fct_relevel(ratio, ratiosFull))
  OTUnums <- unique(df$OTUnum)
  foreach(x=OTUnums) %do% {
    p_OTUAbundancesFamily <- df %>%
      filter(OTUnum==x) %>%
      ggplot() +
      geom_line(aes(x=ratio, y=logAbundance, color=line, group=line)) +
      geom_point(aes(x=ratio, y=logAbundance, color=line), alpha=0.5) +
      geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
      xlab("Ratio") +
      ylab("Relative abundance (log)") +
      ylim(c(-4.1, 0)) +
      ggtitle(paste0(family, "-", x, "-p", passageNum)) +
      scale_color_manual(name="Mixture type",
                         labels=c("theoretical-1"="Theoretical 1",
                                  "theoretical-2"="Theoretical 2",
                                  "actual-1"="Actual 1",
                                  "actual-2"="Actual 2",
                                  "actual-3"="Actual 3"), 
                         values=c("theoretical-1"="chartreuse2", 
                                  "theoretical-2"="green1", 
                                  "actual-1"="slateblue4", 
                                  "actual-2"="steelblue4",
                                  "actual-3"="royalblue2")) +
      DEFAULTS.THEME_PRES +
      facet_wrap(~combo, scales="free") +
      theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
            strip.text.y=element_text(angle=0))
    plotName <- paste("OTU_abundances", family, x, sep="_")
    save_plot(paste0(plotdir104, "/", plotName, "-P", passageNum, ".png"),
              p_OTUAbundancesFamily, nrow=1.5, ncol=2)
  }
}
# Plot individual abundances for all OTUs.
sapply(familyList, FUN=plotFamilyOTUAbundances104, 3)
sapply(familyList, FUN=plotFamilyOTUAbundances104, 5)

# Function to plot OTU abundances for each family and combo, filtered to 10^-3.
plotFamilyOTUAbundances103 <- function(family, passageNum) {
  df <- theoretical_mixtures_df %>% 
    filter(filterRelAb3==TRUE & Family==family & passage==passageNum) %>%
    # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
    select(OTU, OTUnum, ratio, combo, line, logAbundance) %>%
    group_by(combo, line) %>%
    complete(nesting(OTU, OTUnum), 
             ratio=ratiosFull,
             fill=list(logAbundance=NA)) %>%
    ungroup() %>%
    mutate(ratio=fct_relevel(ratio, ratiosFull))
  OTUnums <- unique(df$OTUnum)
  foreach(x=OTUnums) %do% {
    p_OTUAbundancesFamily <- df %>%
      filter(OTUnum==x) %>%
      ggplot() +
      geom_line(aes(x=ratio, y=logAbundance, color=line, group=line)) +
      geom_point(aes(x=ratio, y=logAbundance, color=line), alpha=0.5) +
      geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
      xlab("Ratio") +
      ylab("Relative abundance (log)") +
      ylim(c(-3.1, 0)) +
      ggtitle(paste0(family, "-", x, "-p", passageNum)) +
      scale_color_manual(name="Mixture type",
                         labels=c("theoretical-1"="Theoretical 1",
                                  "theoretical-2"="Theoretical 2",
                                  "actual-1"="Actual 1",
                                  "actual-2"="Actual 2",
                                  "actual-3"="Actual 3"), 
                         values=c("theoretical-1"="chartreuse2", 
                                  "theoretical-2"="green1", 
                                  "actual-1"="slateblue4", 
                                  "actual-2"="steelblue4",
                                  "actual-3"="royalblue2")) +
      DEFAULTS.THEME_PRES +
      facet_wrap(~combo, scales="free") +
      theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
            strip.text.y=element_text(angle=0))
    plotName <- paste("OTU_abundances", family, x, sep="_")
    save_plot(paste0(plotdir103, "/", plotName, "-P", passageNum, ".png"),
              p_OTUAbundancesFamily, nrow=1.5, ncol=2)
  }
}
# Plot individual abundances for all OTUs.
sapply(familyList, FUN=plotFamilyOTUAbundances103, 3)
sapply(familyList, FUN=plotFamilyOTUAbundances103, 5)

# Function to plot OTU abundances for each family and combo, filtered to presence in one parent community.
plotFamilyOTUAbundancesOneParent <- function(family, passageNum) {
  df <- theoretical_mixtures_df %>% 
    filter(filterOneParent==TRUE & Family==family & passage==passageNum) %>%
    # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
    select(OTU, OTUnum, ratio, combo, line, logAbundance) %>%
    group_by(combo, line) %>%
    complete(nesting(OTU, OTUnum), 
             ratio=ratiosFull,
             fill=list(logAbundance=NA)) %>%
    ungroup() %>%
    mutate(ratio=fct_relevel(ratio, ratiosFull))
  OTUnums <- unique(df$OTUnum)
  foreach(x=OTUnums) %do% {
    p_OTUAbundancesFamily <- df %>%
      filter(OTUnum==x) %>%
      ggplot() +
      geom_line(aes(x=ratio, y=logAbundance, color=line, group=line)) +
      geom_point(aes(x=ratio, y=logAbundance, color=line), alpha=0.5) +
      geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
      xlab("Ratio") +
      ylab("Relative abundance (log)") +
      ylim(c(-3.1, 0)) +
      ggtitle(paste0(family, "-", x, "-p", passageNum)) +
      scale_color_manual(name="Mixture type",
                         labels=c("theoretical-1"="Theoretical 1",
                                  "theoretical-2"="Theoretical 2",
                                  "actual-1"="Actual 1",
                                  "actual-2"="Actual 2",
                                  "actual-3"="Actual 3"), 
                         values=c("theoretical-1"="chartreuse2", 
                                  "theoretical-2"="green1", 
                                  "actual-1"="slateblue4", 
                                  "actual-2"="steelblue4",
                                  "actual-3"="royalblue2")) +
      DEFAULTS.THEME_PRES +
      facet_wrap(~combo, scales="free") +
      theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
            strip.text.y=element_text(angle=0))
    plotName <- paste("OTU_abundances", family, x, sep="_")
    save_plot(paste0(plotdirOneParent, "/", plotName, "-P", passageNum, ".png"),
              p_OTUAbundancesFamily, nrow=1.5, ncol=2)
  }
}
# Plot individual abundances for all OTUs.
sapply(familyList, FUN=plotFamilyOTUAbundancesOneParent, 3)
sapply(familyList, FUN=plotFamilyOTUAbundancesOneParent, 5)

# Function to plot OTU abundances for each family and combo, filtered to high abundance in one parent.
plotFamilyOTUAbundancesHighAbundance <- function(family, passageNum) {
  df <- theoretical_mixtures_df %>% 
    filter(filterHighAbundance==TRUE & Family==family & passage==passageNum) %>%
    # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
    select(OTU, OTUnum, ratio, combo, line, logAbundance) %>%
    group_by(combo, line) %>%
    complete(nesting(OTU, OTUnum), 
             ratio=ratiosFull,
             fill=list(logAbundance=NA)) %>%
    ungroup() %>%
    mutate(ratio=fct_relevel(ratio, ratiosFull))
  OTUnums <- unique(df$OTUnum)
  foreach(x=OTUnums) %do% {
    p_OTUAbundancesFamily <- df %>%
      filter(OTUnum==x) %>%
      ggplot() +
      geom_line(aes(x=ratio, y=logAbundance, color=line, group=line)) +
      geom_point(aes(x=ratio, y=logAbundance, color=line), alpha=0.5) +
      geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
      xlab("Ratio") +
      ylab("Relative abundance (log)") +
      ylim(c(-3.1, 0)) +
      ggtitle(paste0(family, "-", x, "-p", passageNum)) +
      scale_color_manual(name="Mixture type",
                         labels=c("theoretical-1"="Theoretical 1",
                                  "theoretical-2"="Theoretical 2",
                                  "actual-1"="Actual 1",
                                  "actual-2"="Actual 2",
                                  "actual-3"="Actual 3"), 
                         values=c("theoretical-1"="chartreuse2", 
                                  "theoretical-2"="green1", 
                                  "actual-1"="slateblue4", 
                                  "actual-2"="steelblue4",
                                  "actual-3"="royalblue2")) +
      DEFAULTS.THEME_PRES +
      facet_wrap(~combo, scales="free") +
      theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
            strip.text.y=element_text(angle=0))
    plotName <- paste("OTU_abundances", family, x, sep="_")
    save_plot(paste0(plotdirHighAbundance, "/", plotName, "-P", passageNum, ".png"),
              p_OTUAbundancesFamily, nrow=1.5, ncol=2)
  }
}
# Plot individual abundances for all OTUs.
sapply(familyList, FUN=plotFamilyOTUAbundancesHighAbundance, 3)
sapply(familyList, FUN=plotFamilyOTUAbundancesHighAbundance, 5)

# High abundance one parent. ----------------------------------------------

plotdirHighAbundanceOneParent <- paste0(plotdir, "/trajectoryPlots-highAbundanceOneParent")
plotdirHighAbundanceOneParentScatter <- paste0(plotdir, "/scatterPlots-highAbundanceOneParent")
plotdirParentDiff <- paste0(plotdir, "/trajectoryPlots-parentDiff-noTheoretical")

# Plot OTU abundances for each OTU, faceted by combo.
# Filtered to high abundance in one parent and present in only one parent.
plotFamilyOTUAbundancesHighAbundanceOneParent <- function(family, passageNum) {
  df <- theoretical_mixtures_HighAbundanceOneParent %>% 
    filter(Family==family & passage==passageNum)
    OTUnums <- unique(df$OTUnum)
  foreach(x=OTUnums) %do% {
    p_OTUAbundancesFamily <- df %>%
      filter(OTUnum==x) %>%
      ggplot() +
      geom_line(aes(x=ratio, y=logAbundance, color=line, group=line), alpha=0.6) +
      geom_point(aes(x=ratio, y=logAbundance, color=line, shape=pointType), alpha=0.8) +
      geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
      xlab("Ratio") +
      ylab("Relative abundance (log)") +
      ylim(c(-3, 0)) +
      ggtitle(paste0(family, "-", x, "-p", passageNum)) +
      scale_color_manual(name="Mixture type",
                         labels=c("theoretical-1"="Theoretical 1",
                                  "theoretical-2"="Theoretical 2",
                                  "actual-1"="Actual 1",
                                  "actual-2"="Actual 2",
                                  "actual-3"="Actual 3"), 
                         values=c("theoretical-1"="#313695", 
                                  "theoretical-2"="#4575B4", 
                                  "actual-1"="#D73027", 
                                  "actual-2"="#F46D43",
                                  "actual-3"="#FDAE61")) +
      scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
      DEFAULTS.THEME_PRES +
      facet_wrap(~combo, scales="free") +
      theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
            strip.text.y=element_text(angle=0))
    plotName <- paste("OTU_abundances", family, x, sep="_")
    save_plot(paste0(plotdirHighAbundanceOneParent, "/", plotName, "-P", passageNum, ".png"),
              p_OTUAbundancesFamily, nrow=1.5, ncol=2)
  }
}
sapply(familyList, FUN=plotFamilyOTUAbundancesHighAbundanceOneParent, 3)
sapply(familyList, FUN=plotFamilyOTUAbundancesHighAbundanceOneParent, 5)

# Plot OTU abundances for each combo, faceted by OTU.

categories <- c("dose-dependent","overcolonizing","undercolonizing","noisy","flat","DD-under","DD-over")

# Use this function to plot ASVs in each category separately.
plotComboOTUAbundancesHighAbundanceOneParent <- function(currentCombo, passageNum, currCategory) {
  df <- theoretical_mixtures_df_categorized %>%
    filter(category==currCategory) %>% 
    filter(combo==currentCombo & passage==passageNum) %>% 
    mutate(FamilyOTU=paste(Family, OTUnum, sep="-"))
  if (nrow(df)>0) {
  p_OTUAbundancesCombo <- df %>%
    ggplot() +
    geom_line(aes(x=ratio, y=logAbundance, color=line, group=line), alpha=0.6) +
    geom_point(aes(x=ratio, y=logAbundance, color=line, shape=pointType), alpha=0.8) +
    geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
    xlab("Ratio") +
    ylab("Relative abundance (log)") +
    ylim(c(-3, 0)) +
    ggtitle(paste0(currentCombo, "-p", passageNum)) +
    scale_color_manual(name="Mixture type",
                       labels=c("theoretical-1"="Theoretical 1",
                                "theoretical-2"="Theoretical 2",
                                "actual-1"="Actual 1",
                                "actual-2"="Actual 2",
                                "actual-3"="Actual 3"), 
                       values=c("theoretical-1"="#313695", 
                                "theoretical-2"="#4575B4", 
                                "actual-1"="#D73027", 
                                "actual-2"="#F46D43",
                                "actual-3"="#FDAE61")) +
    scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
    DEFAULTS.THEME_PRES +
    facet_wrap(~FamilyOTU, scales="free") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
          strip.text.y=element_text(angle=0))
    plotName <- paste("OTU_abundances", currentCombo, sep="_")
    save_plot(paste0(plotdirParentDiff, "/", plotName, "-P", passageNum, "_", currCategory, ".png"),
              p_OTUAbundancesCombo, nrow=7, ncol=7, limitsize=FALSE)
  }
}
#sapply(combos, FUN=plotComboOTUAbundancesHighAbundanceOneParent, 3)
#sapply(combos, FUN=plotComboOTUAbundancesHighAbundanceOneParent, 5)

foreach(x=categories) %do% {
  sapply(combos, FUN=plotComboOTUAbundancesHighAbundanceOneParent, 5, x)
}

# Plot family and ASV-level abundances side-by-side to look for neutral ASVs.
familyComboPassages <- theoretical_mixtures_HighAbundanceOneParent %>% 
  mutate(familyComboPassage=paste0(Family,".",combo,":",passage))
familyComboPassages <- unique(familyComboPassages$familyComboPassage)

foreach(x=familyComboPassages) %do% {
  family=sub("\\..*","",x)
  currCombo=sub(".*\\.","",sub(":.*","",x))
  passageNum=sub(".*:","",x)
  dfTemp <- theoretical_mixtures_HighAbundanceOneParent %>% 
    filter(Family==family & combo==currCombo & passage==passageNum) %>% 
    mutate(FamilyOTU=paste(Family, OTUnum, sep="-"),
           abundance=logAbundance,
           statType=FamilyOTU) %>%
    select(!FamilyOTU) %>% 
    rbind(theoretical_mixtures_HighAbundanceOneParent %>% 
            filter(Family==family & combo==currCombo & passage==passageNum) %>% 
            mutate(abundance=logFamilyAbundance,
                   statType=Family))
  p_FamilyOTUAbundancesSplit <- dfTemp %>% 
    ggplot() +
    geom_line(aes(x=ratio, y=abundance, color=line, group=line), alpha=0.6) +
    geom_point(aes(x=ratio, y=abundance, color=line, shape=pointType), alpha=0.8) +
    geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
    xlab("Ratio") +
    ylab("Relative abundance (log)") +
    ylim(c(-3, 0)) +
    ggtitle(paste0(family, "-", currCombo, "-p", passageNum)) +
    scale_color_manual(name="Mixture type",
                       labels=c("theoretical-1"="Theoretical 1",
                                "theoretical-2"="Theoretical 2",
                                "actual-1"="Actual 1",
                                "actual-2"="Actual 2",
                                "actual-3"="Actual 3"), 
                       values=c("theoretical-1"="#313695", 
                                "theoretical-2"="#4575B4", 
                                "actual-1"="#D73027", 
                                "actual-2"="#F46D43",
                                "actual-3"="#FDAE61")) +
    scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
    DEFAULTS.THEME_PRES +
    facet_wrap(~statType) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
          strip.text.y=element_text(angle=0))
  plotName <- paste0("FamilyOTUabundances_",family,"-",currCombo,"-p",passageNum)
  save_plot(paste0(plotdirHighAbundanceOneParent, "/", plotName, ".png"),
            p_FamilyOTUAbundancesSplit, nrow=1.75, ncol=2)
}

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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~plot) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none")
p_OTUAbundancesExampleFlat
#save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleFlat.png"),
#          p_OTUAbundancesExampleFlat, base_width=1.6, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/OTU_abundance_ExampleFlat.pdf"),
          p_OTUAbundancesExampleFlat, base_width=1.7, base_height=1.6)

trajectoryLegend <- get_legend(p_OTUAbundancesExampleFlat)
save_plot(paste0(plotdirParentDiff, "/OTUtrajectoriesLegend.pdf"), 
          trajectoryLegend, base_width=0.8, base_height=0.7)

p_OTUAbundancesExampleFlat2 <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(passage==5 & combo=="XFA-XFB" & Family=="Bacteroidaceae" & OTUnum==10) %>% 
  mutate(plot="Flat\nXFA-XFB\nBacteroidaceae-10") %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=line, group=line), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=line, shape=pointType), alpha=0.8, size=0.5) +
  geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  ylab("Relative abundance (log10)") +
  ylim(c(-3, 0)) +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("theoretical-1"="#313695", 
                              "theoretical-2"="#4575B4", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=5.5)) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  facet_grid(~plot)
p_OTUAbundancesExampleFlat2
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleFlat2.pdf"),
          p_OTUAbundancesExampleFlat2, base_width=3, base_height=1.75)

p_OTUAbundancesExampleNoisy <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XBA-XCA" & Family=="Tannerellaceae" & OTUnum==5) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=line, group=line), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=line, shape=pointType), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("theoretical-1"="#313695", 
                              "theoretical-2"="#4575B4", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  #DEFAULTS.THEME_PRINT +
  DEFAULTS.THEME_PRES +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
p_OTUAbundancesExampleNoisy
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleNoisy.png"),
          p_OTUAbundancesExampleNoisy, base_width=4, base_height=2.5)

p_OTUAbundancesExampleNoisy2 <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XBA-XBB" & Family=="Rikenellaceae" & OTUnum==12) %>% 
  mutate(plot="Noisy\nXBA-XBB\nRikenellaceae-12") %>%
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
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~plot)
p_OTUAbundancesExampleNoisy2
#save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleNoisy2.png"),
#          p_OTUAbundancesExampleNoisy2, base_width=1.6, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/OTU_abundance_ExampleNoisy.pdf"),
          p_OTUAbundancesExampleNoisy2, base_width=1.7, base_height=2)

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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~plot)
p_OTUAbundancesExampleNoisy3
#save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleNoisy2.png"),
#          p_OTUAbundancesExampleNoisy2, base_width=1.6, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/OTU_abundance_ExampleNoisy.pdf"),
          p_OTUAbundancesExampleNoisy3, base_width=1.7, base_height=1.6)

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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~plot)
p_OTUAbundancesExampleUndercolonizing
#save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleUndercolonizing.png"),
#          p_OTUAbundancesExampleUndercolonizing, base_width=1.6, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/OTU_abundance_ExampleUndercolonizing.pdf"),
          p_OTUAbundancesExampleUndercolonizing, base_width=1.7, base_height=1.6)

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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~plot)
p_OTUAbundancesExampleOvercolonizing
#save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleOvercolonizingTheoretical.png"),
#          p_OTUAbundancesExampleOvercolonizingTheoretical, base_width=2.5, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/OTU_abundance_ExampleOvercolonizing.pdf"),
          p_OTUAbundancesExampleOvercolonizing, base_width=1.7, base_height=1.6)

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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~plot)
p_OTUAbundancesExampleDoseDependent
#save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleDoseDependent.png"),
#          p_OTUAbundancesExampleDoseDependent, base_width=2.5, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/OTU_abundance_ExampleDoseDependent.pdf"),
          p_OTUAbundancesExampleDoseDependent, base_width=1.7, base_height=1.6)

p_OTUAbundancesExampleCombined <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 &
           ((combo=="XDA-XDB" & Family=="Bacteroidaceae" & OTUnum==10) |
              (combo=="XBA-XBB" & Family=="Rikenellaceae" & OTUnum==12) |
              (combo=="XFA-XBA" & Family=="Bacteroidaceae" & OTUnum==9) |
              (combo=="XFA-XBA" & Family=="Lachnospiraceae" & OTUnum==115) |
              (combo=="XCA-XDA" & Family=="Bacteroidaceae" & OTUnum==9))) %>% 
  ungroup() %>% 
  mutate(plot=case_when(
    combo=="XCA-XDA" & Family=="Bacteroidaceae" & OTUnum==9 ~ 
      paste0("Dose-dependent\n",combo,"\n",Family,"-",OTUnum),
    combo=="XFA-XBA" & Family=="Lachnospiraceae" & OTUnum==115 ~ 
      paste0("Overcolonizing\n",combo,"\n",Family,"-",OTUnum),
    combo=="XFA-XBA" & Family=="Bacteroidaceae" & OTUnum==9 ~ 
      paste0("Undercolonizing\n",combo,"\n",Family,"-",OTUnum),
    combo=="XDA-XDB" & Family=="Bacteroidaceae" & OTUnum==10 ~ 
      paste0("Flat\n",combo,"\n",Family,"-",OTUnum),
    combo=="XBA-XBB" & Family=="Rikenellaceae" & OTUnum==12 ~ 
      paste0("Noisy\n",combo,"\n",Family,"-",OTUnum)),
    plot=fct_relevel(plot,
                     "Dose-dependent\nXCA-XDA\nBacteroidaceae-9",
                     "Overcolonizing\nXFA-XBA\nLachnospiraceae-115",
                     "Undercolonizing\nXFA-XBA\nBacteroidaceae-9",
                     "Flat\nXDA-XDB\nBacteroidaceae-10",
                     "Noisy\nXBA-XBB\nRikenellaceae-12")) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=line, group=line), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=line, shape=pointType), alpha=0.8, size=0.5) +
  geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  ylab("Relative abundance (log10)") +
  ylim(c(-3, 0)) +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("theoretical-1"="#313695", 
                              "theoretical-2"="#4575B4", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  DEFAULTS.THEME_PRINT +
  #DEFAULTS.THEME_PRES +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  facet_grid(~plot)
p_OTUAbundancesExampleCombined
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleCombined.pdf"),
          p_OTUAbundancesExampleCombined, base_width=7, base_height=2)
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleCombined.png"),
          p_OTUAbundancesExampleCombined, nrow=0.75, ncol=1.5)

# Neutral examples. -------------------------------------------------------

# Combined example with Bacteroidaceae and Streptococcaceae.
p_OTUAbundancesExampleNeutral <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & ((combo=="XFA-XFB" & Family=="Streptococcaceae" & OTUnum %in% c(7,17)) | 
              (combo=="XCA-XDA" & Family=="Bacteroidaceae" & OTUnum %in% c(7,9)))) %>% 
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
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2")) %>% 
  ungroup() %>% 
  mutate(plot = case_when(OTUID=="Strep-7" ~ "XFA-XFB\nStrep-7",
                          OTUID=="Strep-17" ~ "XFA-XFB\nStrep-17",
                          OTUID=="Bacte-7" ~ "XCA-XDA\nBacte-7",
                          OTUID=="Bacte-9" ~ "XCA-XDA\nBacte-9"),
         plot = fct_relevel(plot, c("XCA-XDA\nBacte-7","XCA-XDA\nBacte-9",
                                    "XFA-XFB\nStrep-7","XFA-XFB\nStrep-17"))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~plot)
p_OTUAbundancesExampleNeutral
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleNeutral.pdf"),
          p_OTUAbundancesExampleNeutral, base_width=3, base_height=2.4)

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
         plot = fct_relevel(plot, c("L. garvieae","L. lactis",
                                    "E. faecalis","E. casseliflavus"))) %>% 
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
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
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleNeutralStrepEnteC-p3.png"),
          p_OTUAbundancesExampleNeutralStrepEnteC, base_width=3, base_height=2.4)

# EnteC2 and EnteC3 example.
p_OTUAbundancesExampleEnteC <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XFA-XFB" & Family=="Enterococcaceae" & OTUnum %in% c(2,3)) %>% 
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
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2")) %>% 
  ungroup() %>% 
  mutate(plot = case_when(OTUID=="EnteC-2" ~ "D1/D2\nEnteC-2",
                          OTUID=="EnteC-3" ~ "D1/D2\nEnteC-3")) %>%  
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~plot)
p_OTUAbundancesExampleEnteC
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleEnteC.pdf"),
          p_OTUAbundancesExampleEnteC, base_width=3, base_height=1.8)

# ASV patterns over mixtures examples. ------------------------------------

# Bacteroidaceae 10.
p_OTUAbundancesExampleBacte10 <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo %in% c("XCA-XCB","XDA-XFA") & Family=="Bacteroidaceae" & OTUnum==10) %>% 
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
         plot=ifelse(
           combo=="XCA-XCB", paste0(OTUID,"\nXCA-XCB\nDose-dependent"),
           paste0(OTUID,"\nXDA-XFA\nOvercolonizing"))) %>% 
  ungroup() %>%
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~category)
p_OTUAbundancesExampleBacte10
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleBacte10.png"),
          p_OTUAbundancesExampleBacte10, base_width=3.5, base_height=2)

# Lachnospiraceae 219.
p_OTUAbundancesExampleLachn219 <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo %in% c("XBA-XCA","XCA-XCB") & Family=="Lachnospiraceae" & OTUnum==219) %>% 
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
         plot=ifelse(
           combo=="XBA-XCA", paste0(OTUID,"\nXBA-XCA\nDose-dependent"),
           paste0(OTUID,"\nXCA-XCB\nUndercolonizing")),
         plot = fct_relevel(plot, c("Lachn-219\nXCA-XCB\nUndercolonizing", "Lachn219\nXBA-XCA\nDose-dependent"))) %>% 
  ungroup() %>%
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~plot)
p_OTUAbundancesExampleLachn219
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleLachn219.png"),
          p_OTUAbundancesExampleLachn219, base_width=3.5, base_height=2)

# Example Streptococcaceae plot. ------------------------------------------
p_OTUAbundancesExampleStrep <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XFA-XFB" & Family=="Streptococcaceae" & OTUnum %in% c(7,17)) %>% 
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
                                 ratio=="0:1" & inoculationReplicate==2 ~ "rightParent2")) %>% 
  ungroup() %>% 
  mutate(OTUID = fct_relevel(OTUID, c("Strep-7","Strep-17"))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~OTUID)
p_OTUAbundancesExampleStrep
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleStrep.png"),
          p_OTUAbundancesExampleStrep, base_width=3.5, base_height=1.6)

# Streptococcaceae complementation with each other.
p_OTUAbundancesExampleStrep7Strep17 <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XFA-XFB" & Family=="Streptococcaceae" & OTUnum %in% c(7, 17)) %>% 
  group_by(line, ratio) %>% 
  mutate(sum_rel_abundance = ifelse(rel_abundance <= 10^-3, 0, rel_abundance),
         sum_rel_abundance = sum(rel_abundance),
         sumLogAbundance = log10(sum_rel_abundance)) %>%  
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
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=sumLogAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=sumLogAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#313695", 
                              "theoretical-2"="#4575B4", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  ggtitle("Strep7-Strep17") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT
p_OTUAbundancesExampleStrep7Strep17
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleStrep7Strep17.png"),
          p_OTUAbundancesExampleStrep7Strep17, base_width=4.5, base_height=2.8)

# Streptococcaceae complementation with Entero ASV 3.
p_OTUAbundancesExampleStrep17EnteC3 <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XFA-XFB" & 
           ((Family=="Streptococcaceae" & OTUnum==17) | (Family=="Enterococcaceae" & OTUnum==3))) %>%  
  group_by(line, ratio) %>% 
  mutate(sum_rel_abundance = ifelse(rel_abundance <= 10^-3, 0, rel_abundance),
         sum_rel_abundance = sum(rel_abundance),
         sumLogAbundance = log10(sum_rel_abundance)) %>%  
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
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=sumLogAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=sumLogAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#313695", 
                              "theoretical-2"="#4575B4", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  ggtitle("Strep17-EnteC3") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT
p_OTUAbundancesExampleStrep17EnteC3
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleStrep17EnteC3.png"),
          p_OTUAbundancesExampleStrep17EnteC3, base_width=4.5, base_height=2.8)

# Streptococcaceae complementation with Entero ASVs 2 and 3.
p_OTUAbundancesExampleStrep17EnteC23 <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XFA-XFB" & 
           ((Family=="Streptococcaceae" & OTUnum==17) | (Family=="Enterococcaceae" & OTUnum %in% c(2,3)))) %>%  
  group_by(line, ratio) %>% 
  mutate(sum_rel_abundance = ifelse(rel_abundance <= 10^-3, 0, rel_abundance),
         sum_rel_abundance = sum(rel_abundance),
         sumLogAbundance = log10(sum_rel_abundance)) %>%  
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
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=sumLogAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=sumLogAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#313695", 
                              "theoretical-2"="#4575B4", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  ggtitle("Strep17-EnteC2/3") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT
p_OTUAbundancesExampleStrep17EnteC23
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleStrep17EnteC23.png"),
          p_OTUAbundancesExampleStrep17EnteC23, base_width=4.5, base_height=2.8)

# Streptococcaceae complementation with all Entero ASVs.
p_OTUAbundancesExampleStrep17EnteCAll <- theoretical_mixtures_df_categorized %>%
  filter(passage==5 & combo=="XFA-XFB" & 
           ((Family=="Streptococcaceae" & OTUnum==17) | (Family=="Enterococcaceae" & category!="lowAbundance"))) %>%  
  group_by(line, ratio) %>% 
  mutate(sum_rel_abundance = ifelse(rel_abundance <= 10^-3, 0, rel_abundance),
         sum_rel_abundance = sum(rel_abundance),
         sumLogAbundance = log10(sum_rel_abundance)) %>%  
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
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=sumLogAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=sumLogAbundance, color=lineGroup, shape=pointType), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#313695", 
                              "theoretical-2"="#4575B4", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  ggtitle("Strep17-EnteCAll") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT
p_OTUAbundancesExampleStrep17EnteCAll
save_plot(paste0(plotdirParentDiff, "/OTU_abundances_ExampleStrep17EnteCAll.png"),
          p_OTUAbundancesExampleStrep17EnteCAll, base_width=4.5, base_height=2.8)


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
  mutate(plot = ifelse(passage==3,"Passage 3\nDose-dependent", 
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
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
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
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
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
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
  scale_shape_manual(values=c("circle"=16,"openCircle"=1), guide="none") +
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

# Plot scatterplots for high abundance one parent ASVs. -------------------

# Create abundance dataframe with theoretical mean in one column and actual mean in another.
ASV_sum_abundances <- full_join(
  t_sum %>% 
    filter(ratio %in% ratios) %>% 
    group_by(OTU, passage, combo, ratio) %>% 
    slice(1) %>% 
    select(OTU, OTUnum, Family, passage, combo, ratio, avgAbundance, rangeAbundance) %>% 
    rename(t_avgAbundance=avgAbundance, t_rangeAbundance=rangeAbundance),
  a_sum %>%
    group_by(OTU, passage, combo, ratio) %>% 
    slice(1) %>% 
    select(OTU, passage, combo, ratio, avgAbundance, rangeAbundance) %>% 
    rename(a_avgAbundance=avgAbundance, a_rangeAbundance=rangeAbundance),
  by=c("OTU", "passage", "combo", "ratio"))

# Function to plot OTU scatterplots for each family, faceted by OTU.
plotFamilyOTUScatterplotsHighAbundanceOneParent <- function(family, passageNum) {
  df <- ASV_sum_abundances %>% 
    filter(Family==family & passage==passageNum)
  # Annotate OTUs for plotting in separate plots.
  OTUnums <- unique(df$OTUnum)
  foreach(x=OTUnums) %do% {
    p_OTUAbundancesFamilyScatterplot <- df %>%
      filter(OTUnum==x) %>%
      ggplot() +
      geom_point(aes(x=t_avgAbundance, y=a_avgAbundance, color=ratio)) +
      geom_errorbarh(aes(y=a_avgAbundance, color=ratio, height=0,
                         xmin=t_avgAbundance-t_rangeAbundance/2, xmax=t_avgAbundance+t_rangeAbundance/2)) +
      geom_errorbar(aes(x=t_avgAbundance, color=ratio,
                        ymin=a_avgAbundance-a_rangeAbundance/2, ymax=a_avgAbundance+a_rangeAbundance/2)) +
      scale_color_brewer(palette="Reds", name="Ratio") +
      geom_abline(aes(slope=1, intercept=0)) +
      geom_hline(aes(yintercept=-4), linetype="dotted") +
      geom_vline(aes(xintercept=-4), linetype="dotted") +
      ylab("Average log actual relative abundance") +
      xlab("Average log theoretical relative abundance") +
      ylim(c(-3,0)) +
      xlim(c(-3,0)) +
      ggtitle(paste0(family, "-", x, "-p", passageNum)) +
      DEFAULTS.THEME_PRES +
      facet_wrap(~combo, scales="free")
    plotName <- paste("OTU_scatterplot", family, x, sep="_")
    save_plot(paste0(plotdirHighAbundanceOneParentScatter, "/", plotName, "-P", passageNum, ".png"),
              p_OTUAbundancesFamilyScatterplot, nrow=1.5, ncol=2)
  }
}
# Plot individual abundance scatterplots for all OTUs.
sapply(familyList, FUN=plotFamilyOTUScatterplotsHighAbundanceOneParent, 3)
sapply(familyList, FUN=plotFamilyOTUScatterplotsHighAbundanceOneParent, 5)

# Function to plot OTU scatterplots for each combo, faceted by OTU.
plotComboOTUScatterplotsHighAbundanceOneParent <- function(currentCombo, passageNum) {
  df <- ASV_sum_abundances %>% 
    filter(combo==currentCombo & passage==passageNum) %>% 
    mutate(FamilyOTU=paste(Family, OTUnum, sep="-"))
  p_OTUAbundancesComboScatterplot <- df %>%
    ggplot() +
    geom_point(aes(x=t_avgAbundance, y=a_avgAbundance, color=ratio)) +
    geom_errorbarh(aes(y=a_avgAbundance, color=ratio, height=0,
                       xmin=t_avgAbundance-t_rangeAbundance/2, xmax=t_avgAbundance+t_rangeAbundance/2)) +
    geom_errorbar(aes(x=t_avgAbundance, color=ratio,
                      ymin=a_avgAbundance-a_rangeAbundance/2, ymax=a_avgAbundance+a_rangeAbundance/2)) +
    scale_color_brewer(palette="Reds", name="Ratio") +
    geom_abline(aes(slope=1, intercept=0)) +
    geom_hline(aes(yintercept=-4), linetype="dotted") +
    geom_vline(aes(xintercept=-4), linetype="dotted") +
    ylab("Average log actual relative abundance") +
    xlab("Average log theoretical relative abundance") +
    ylim(c(-3,0)) +
    xlim(c(-3,0)) +
    ggtitle(paste0(currentCombo, "-p", passageNum)) +
    DEFAULTS.THEME_PRES +
    facet_wrap(~FamilyOTU, scales="free")
  plotName <- paste("OTU_scatterplot", currentCombo, sep="_")
  save_plot(paste0(plotdirHighAbundanceOneParentScatter, "/", plotName, "-P", passageNum, ".png"),
            p_OTUAbundancesComboScatterplot, nrow=1.5, ncol=2)
}
# Plot individual abundance scatterplots for all OTUs.
sapply(combos, FUN=plotComboOTUScatterplotsHighAbundanceOneParent, 3)
sapply(combos, FUN=plotComboOTUScatterplotsHighAbundanceOneParent, 5)


# Plot all OTU scatterplots together for each combo.
plotOTUScatterplotsHighAbundanceOneParent <- function(currentCombo, passageNum) {
  df <- ASV_sum_abundances %>% 
    filter(passage==passageNum)
  p_OTUAbundancesScatterplot <- df %>% 
    ggplot() +
    geom_point(aes(x=t_avgAbundance, y=a_avgAbundance, color=ratio)) +
    geom_errorbarh(aes(y=a_avgAbundance, color=ratio, height=0,
                       xmin=t_avgAbundance-t_rangeAbundance/2, xmax=t_avgAbundance+t_rangeAbundance/2),
                   alpha=0.3) +
    geom_errorbar(aes(x=t_avgAbundance, color=ratio,
                      ymin=a_avgAbundance-a_rangeAbundance/2, ymax=a_avgAbundance+a_rangeAbundance/2),
                  alpha=0.3) +
    scale_color_brewer(palette="Reds", name="Ratio") +
    geom_abline(aes(slope=1, intercept=0)) +
    geom_hline(aes(yintercept=-4), linetype="dotted") +
    geom_vline(aes(xintercept=-4), linetype="dotted") +
    ylab("Average log actual relative abundance") +
    xlab("Average log theoretical relative abundance") +
    ylim(c(-3,0)) +
    xlim(c(-3,0)) +
    ggtitle(paste0(currentCombo, "-p", passageNum)) +
    DEFAULTS.THEME_PRES
  plotName <- paste("OTU_scatterplot_unfaceted", currentCombo, sep="_")
  save_plot(paste0(plotdirHighAbundanceOneParentScatter, "/", plotName, "-P", passageNum, ".png"),
            p_OTUAbundancesScatterplot, nrow=1.5, ncol=2)
}
sapply(combos, FUN=plotOTUScatterplotsHighAbundanceOneParent, 3)
sapply(combos, FUN=plotOTUScatterplotsHighAbundanceOneParent, 5)
