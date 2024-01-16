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

setwd("../../")
source("analysis/generateMixturesDataframe.R")
plotdir <- "analysis/individualFamilyPlots"
plotdirScatter <- paste0(plotdir, "/scatterplots")
plotdirScatterUnfiltered <- paste0(plotdir, "/scatterplotsUnfiltered")

# Family individual abundance plots. ---------------------------------------

plotdirTrajectory <- paste0(plotdir, "/trajectoryPlots")

# Generate list of all unique families in dataframe to plot over.
theoretical_mixtures_df_statsFiltered_P3 <- theoretical_mixtures_df_statsFiltered %>%
  filter(passage==3)
theoretical_mixtures_df_statsFiltered_P5 <- theoretical_mixtures_df_statsFiltered %>%
  filter(passage==5)
familyListP3 <- na.exclude(unique(theoretical_mixtures_df_statsFiltered_P3$Family))
familyListP5 <- na.exclude(unique(theoretical_mixtures_df_statsFiltered_P5$Family))

# Plot OTu abundances for each family, faceted by combo.
plotIndividualFamilyAbundances <- function(family, passageNum) {
  p_individualFamilyPlot <- theoretical_mixtures_df_statsFiltered %>% 
    filter(Family==family & passage==passageNum) %>%
    # Group individual OTU abundances into family abundances.
    select(Family, OTU, ratio, combo, inoculationReplicate, mixtureType, label, line, rel_abundance) %>%
    group_by(combo, line, ratio) %>%
    # Remove artificial 10^-4 points in theoretical mixtures.
    filter(mixtureType=="actual" | rel_abundance>10^-4) %>% 
    mutate(family_rel_abundance=sum(rel_abundance, na.rm=TRUE), 
           logFamilyAbundance=log10(family_rel_abundance)) %>%
    ungroup() %>% 
    group_by(combo, line) %>% 
    # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
    complete(Family, 
             ratio=ratiosFull,
             fill=list(rel_abundance=NA)) %>% 
    ungroup() %>%
    mutate(ratio=fct_relevel(ratio, ratiosFull)) %>% 
    ggplot() +
    geom_line(aes(x=ratio, y=logFamilyAbundance, color=line, group=line)) +
    geom_point(aes(x=ratio, y=logFamilyAbundance, color=line), alpha=0.5) +
    geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
    xlab("Ratio") +
    ylab("Relative abundance (log)") +
    ggtitle(paste0(family, "-p", passageNum)) +
    ylim(c(-4.1, 0)) +
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
          strip.text.y=element_blank())
  plotName <- paste("individualFamilyPlot", family, sep="_")
  save_plot(paste0(plotdirTrajectory, "/", plotName, "-P", passageNum, ".png"), 
            p_individualFamilyPlot, nrow=1.5, ncol=2)
}
sapply(familyListP3, FUN=plotIndividualFamilyAbundances, 3)
sapply(familyListP5, FUN=plotIndividualFamilyAbundances, 5)

# Family individual abundance plots, unfiltered. --------------------------

plotdirTrajectoryUnfiltered <- paste0(plotdir, "/trajectoryPlotsUnfiltered")

# Generate list of all unique families from unfiltered dataframe.
theoretical_mixtures_df_stats_P3 <- theoretical_mixtures_df_stats %>%
  filter(passage==3)
theoretical_mixtures_df_stats_P5 <- theoretical_mixtures_df_stats %>%
  filter(passage==5)
familyListUnfilteredP3 <- na.exclude(unique(theoretical_mixtures_df_stats_P3$Family))
familyListUnfilteredP5 <- na.exclude(unique(theoretical_mixtures_df_stats_P5$Family))

# Function to plot abundances for each family and combo.
plotIndividualFamilyAbundancesUnfiltered <- function(family, passageNum) {
  p_individualFamilyPlot <- theoretical_mixtures_df_stats %>% 
    filter(Family==family & passage==passageNum) %>%
    # Group individual OTU abundances into family abundances.
    select(Family, OTU, ratio, combo, inoculationReplicate, mixtureType, label, line, rel_abundance) %>%
    group_by(combo, line, ratio) %>%
    # Remove artificial 10^-4 points in theoretical mixtures.
    filter(mixtureType=="actual" | rel_abundance>10^-4) %>%
    mutate(family_rel_abundance=sum(rel_abundance, na.rm=TRUE), 
           logFamilyAbundance=log10(family_rel_abundance)) %>%
    ungroup() %>% 
    group_by(combo, line) %>% 
    # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
    complete(Family, 
             ratio=ratiosFull,
             fill=list(rel_abundance=NA)) %>% 
    ungroup() %>%
    mutate(ratio=fct_relevel(ratio, ratiosFull)) %>% 
    ggplot() +
    geom_line(aes(x=ratio, y=logFamilyAbundance, color=line, group=line)) +
    geom_point(aes(x=ratio, y=logFamilyAbundance, color=line), alpha=0.5) +
    geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
    xlab("Ratio") +
    ylab("Relative abundance (log)") +
    ggtitle(paste0(family, "-p", passageNum)) +
    ylim(c(-4.1, 0)) +
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
          strip.text.y=element_blank())
  plotName <- paste("individualFamilyPlot", family, sep="_")
  save_plot(paste0(plotdirTrajectoryUnfiltered, "/", plotName, "-P", passageNum, ".png"), 
            p_individualFamilyPlot, nrow=1.5, ncol=2)
}
sapply(familyListUnfilteredP3, FUN=plotIndividualFamilyAbundancesUnfiltered, 3)
sapply(familyListUnfilteredP5, FUN=plotIndividualFamilyAbundancesUnfiltered, 5)

# Family individual abundance plots, filtered to 10^-3. --------------------------

plotdirTrajectory3 <- paste0(plotdir, "/trajectoryPlots3")

# Generate list of all unique families from unfiltered dataframe.
theoretical_mixtures_df_3_P3 <- theoretical_mixtures_df_filtered3 %>%
  filter(passage==3)
theoretical_mixtures_df_3_P5 <- theoretical_mixtures_df_filtered3 %>%
  filter(passage==5)
familyList3P3 <- na.exclude(unique(theoretical_mixtures_df_3_P3$Family))
familyList3P5 <- na.exclude(unique(theoretical_mixtures_df_3_P5$Family))

# Function to plot abundances for each family and combo.
plotIndividualFamilyAbundances3 <- function(family, passageNum) {
  df <- theoretical_mixtures_df_filtered3 %>% 
    filter(Family==family & passage==passageNum) %>%
    # Group individual OTU abundances into family abundances.
    select(Family, OTU, ratio, combo, inoculationReplicate, mixtureType, label, line, rel_abundance) %>%
    group_by(combo, line, ratio) %>%
    mutate(family_rel_abundance=sum(rel_abundance, na.rm=TRUE), 
           logFamilyAbundance=log10(family_rel_abundance)) %>%
    ungroup() %>% 
    group_by(combo, line) %>% 
    # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
    complete(Family, 
             ratio=ratiosFull,
             fill=list(rel_abundance=NA)) %>% 
    ungroup() %>%
    mutate(ratio=fct_relevel(ratio, ratiosFull))
  p_individualFamilyPlot <- df %>% 
    ggplot() +
    geom_line(aes(x=ratio, y=logFamilyAbundance, color=line, group=line)) +
    geom_point(aes(x=ratio, y=logFamilyAbundance, color=line), alpha=0.5) +
    geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
    xlab("Ratio") +
    ylab("Relative abundance (log)") +
    ggtitle(paste0(family, "-p", passageNum)) +
    ylim(c(-3.1, 0)) +
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
          strip.text.y=element_blank())
  nCombo <- n_distinct(df$combo)
  plotName <- paste("individualFamilyPlot", family, sep="_")
  save_plot(paste0(plotdirTrajectory3, "/", plotName, "-P", passageNum, ".png"), 
            p_individualFamilyPlot, nrow=1.5, ncol=2)
}
sapply(familyList3P3, FUN=plotIndividualFamilyAbundances3, 3)
sapply(familyList3P5, FUN=plotIndividualFamilyAbundances3, 5)


# Plot specific family figures for presentation. --------------------------

# Plot XFA-XFB Streptococaceae alongside two other combos to make sizing consistent.
p_Streptococcaceae_XFAXFB <- theoretical_mixtures_df_stats %>% 
  filter(Family=="Streptococcaceae" & passage==3 & combo %in% c("XBA-XBB", "XDA-XDB", "XFA-XFB")) %>%
  # Group individual OTU abundances into family abundances.
  select(Family, OTU, ratio, combo, inoculationReplicate, mixtureType, label, line, rel_abundance) %>%
  group_by(combo, line, ratio) %>%
  # Remove artificial 10^-4 points in theoretical mixtures.
  filter(mixtureType=="actual" | rel_abundance>10^-4) %>%
  mutate(family_rel_abundance=sum(rel_abundance, na.rm=TRUE), 
         logFamilyAbundance=log10(family_rel_abundance)) %>%
  ungroup() %>% 
  group_by(combo, line) %>% 
  # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
  complete(Family, 
           ratio=ratiosFull,
           fill=list(rel_abundance=NA)) %>% 
  ungroup() %>%
  mutate(ratio=fct_relevel(ratio, ratiosFull)) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logFamilyAbundance, color=line, group=line)) +
  geom_point(aes(x=ratio, y=logFamilyAbundance, color=line), alpha=0.5) +
  geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
  xlab("Ratio") +
  ylab("Relative abundance (log)") +
  ggtitle("Streptococaceae-p3") +
  ylim(c(-4.1, 0)) +
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
        strip.text.y=element_blank())
save_plot(paste0(plotdir, "/Streptococcaceae-P3_XFA-XFB.png"), 
          p_Streptococcaceae_XFAXFB, nrow=1.5, ncol=2)

# Plot XCA-XDA Eubacteriaceae alongside two other combos to make sizing consistent.
p_Eubacteriaceae_XCAXDA <- theoretical_mixtures_df_stats %>% 
  filter(Family=="Eubacteriaceae" & passage==3 & combo %in% c("XBA-XBB", "XCA-XCB", "XCA-XDA")) %>%
  # Group individual OTU abundances into family abundances.
  select(Family, OTU, ratio, combo, inoculationReplicate, mixtureType, label, line, rel_abundance) %>%
  group_by(combo, line, ratio) %>%
  # Remove artificial 10^-4 points in theoretical mixtures.
  filter(mixtureType=="actual" | rel_abundance>10^-4) %>%
  mutate(family_rel_abundance=sum(rel_abundance, na.rm=TRUE), 
         logFamilyAbundance=log10(family_rel_abundance)) %>%
  ungroup() %>% 
  group_by(combo, line) %>% 
  # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
  complete(Family, 
           ratio=ratiosFull,
           fill=list(rel_abundance=NA)) %>% 
  ungroup() %>%
  mutate(ratio=fct_relevel(ratio, ratiosFull)) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logFamilyAbundance, color=line, group=line)) +
  geom_point(aes(x=ratio, y=logFamilyAbundance, color=line), alpha=0.5) +
  geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
  xlab("Ratio") +
  ylab("Relative abundance (log)") +
  ggtitle("Eubacteriaceae-p3") +
  ylim(c(-4.1, 0)) +
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
        strip.text.y=element_blank())
save_plot(paste0(plotdir, "/Eubacteriaceae-P3_XCA-XDA.png"), 
          p_Eubacteriaceae_XCAXDA, nrow=1.5, ncol=2)

# Plot XBA-XBB Rikenellaceae alongside two other combos to make sizing consistent.
p_Rikenellaceae_XBAXBB <- theoretical_mixtures_df_stats %>% 
  filter(Family=="Rikenellaceae" & passage==3 & combo %in% c("XBA-XBB", "XCA-XCB", "XDA-XDB")) %>%
  # Group individual OTU abundances into family abundances.
  select(Family, OTU, ratio, combo, inoculationReplicate, mixtureType, label, line, rel_abundance) %>%
  group_by(combo, line, ratio) %>%
  # Remove artificial 10^-4 points in theoretical mixtures.
  filter(mixtureType=="actual" | rel_abundance>10^-4) %>%
  mutate(family_rel_abundance=sum(rel_abundance, na.rm=TRUE), 
         logFamilyAbundance=log10(family_rel_abundance)) %>%
  ungroup() %>% 
  group_by(combo, line) %>% 
  # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
  complete(Family, 
           ratio=ratiosFull,
           fill=list(rel_abundance=NA)) %>% 
  ungroup() %>%
  mutate(ratio=fct_relevel(ratio, ratiosFull)) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logFamilyAbundance, color=line, group=line)) +
  geom_point(aes(x=ratio, y=logFamilyAbundance, color=line), alpha=0.5) +
  geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
  xlab("Ratio") +
  ylab("Relative abundance (log)") +
  ggtitle("Rikenellaceae-p3") +
  ylim(c(-4.1, 0)) +
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
        strip.text.y=element_blank())
save_plot(paste0(plotdir, "/Rikenellaceae-P3_XBA-XBB.png"), 
          p_Rikenellaceae_XBAXBB, nrow=1.5, ncol=2)

# Plot XDA-XDB Butyricicoccaceae alongside two other combos to make sizing consistent.
p_Butyricicoccaceae_XDAXDB <- theoretical_mixtures_df_stats %>% 
  filter(Family=="Butyricicoccaceae" & passage==3 & combo %in% c("XBA-XBB", "XCA-XCB", "XDA-XDB")) %>%
  # Group individual OTU abundances into family abundances.
  select(Family, OTU, ratio, combo, inoculationReplicate, mixtureType, label, line, rel_abundance) %>%
  group_by(combo, line, ratio) %>%
  # Remove artificial 10^-4 points in theoretical mixtures.
  filter(mixtureType=="actual" | rel_abundance>10^-4) %>%
  mutate(family_rel_abundance=sum(rel_abundance, na.rm=TRUE), 
         logFamilyAbundance=log10(family_rel_abundance)) %>%
  ungroup() %>% 
  group_by(combo, line) %>% 
  # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
  complete(Family, 
           ratio=ratiosFull,
           fill=list(rel_abundance=NA)) %>% 
  ungroup() %>%
  mutate(ratio=fct_relevel(ratio, ratiosFull)) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logFamilyAbundance, color=line, group=line)) +
  geom_point(aes(x=ratio, y=logFamilyAbundance, color=line), alpha=0.5) +
  geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
  xlab("Ratio") +
  ylab("Relative abundance (log)") +
  ggtitle("Butyricicoccaceae-p3") +
  ylim(c(-4.1, 0)) +
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
        strip.text.y=element_blank())
save_plot(paste0(plotdir, "/Butyricicoccaceae-P3_XDA-XDB.png"), 
          p_Butyricicoccaceae_XDAXDB, nrow=1.5, ncol=2)

# Plot XFA-XFB Family XI alongside two other combos to make sizing consistent.
p_FamilyXI_XFAXFB <- theoretical_mixtures_df_stats %>% 
  filter(Family=="Family_XI" & passage==5 & combo %in% c("XBA-XBB", "XCA-XCB", "XFA-XFB")) %>%
  # Group individual OTU abundances into family abundances.
  select(Family, OTU, ratio, combo, inoculationReplicate, mixtureType, label, line, rel_abundance) %>%
  group_by(combo, line, ratio) %>%
  # Remove artificial 10^-4 points in theoretical mixtures.
  filter(mixtureType=="actual" | rel_abundance>10^-4) %>%
  mutate(family_rel_abundance=sum(rel_abundance, na.rm=TRUE), 
         logFamilyAbundance=log10(family_rel_abundance)) %>%
  ungroup() %>% 
  group_by(combo, line) %>% 
  # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
  complete(Family, 
           ratio=ratiosFull,
           fill=list(rel_abundance=NA)) %>% 
  ungroup() %>%
  mutate(ratio=fct_relevel(ratio, ratiosFull)) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logFamilyAbundance, color=line, group=line)) +
  geom_point(aes(x=ratio, y=logFamilyAbundance, color=line), alpha=0.5) +
  geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
  xlab("Ratio") +
  ylab("Relative abundance (log)") +
  ggtitle("Family_XI-p5") +
  ylim(c(-4.1, 0)) +
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
        strip.text.y=element_blank())
save_plot(paste0(plotdir, "/FamilyXI-P5_XFA-XFB.png"), 
          p_FamilyXI_XFAXFB, nrow=1.5, ncol=2)

# Plot XDA-XFA Marinifilaceae alongside two other combos to make sizing consistent.
p_Marinifilaceae_XDAXFA <- theoretical_mixtures_df_stats %>% 
  filter(Family=="Marinifilaceae" & passage==3 & combo %in% c("XCA-XCB", "XFA-XFB", "XDA-XFA")) %>%
  # Group individual OTU abundances into family abundances.
  select(Family, OTU, ratio, combo, inoculationReplicate, mixtureType, label, line, rel_abundance) %>%
  group_by(combo, line, ratio) %>%
  # Remove artificial 10^-4 points in theoretical mixtures.
  filter(mixtureType=="actual" | rel_abundance>10^-4) %>%
  mutate(family_rel_abundance=sum(rel_abundance, na.rm=TRUE), 
         logFamilyAbundance=log10(family_rel_abundance)) %>%
  ungroup() %>% 
  group_by(combo, line) %>% 
  # Add NA to ratios where OTU does not appear, to avoid plot artifacts.
  complete(Family, 
           ratio=ratiosFull,
           fill=list(rel_abundance=NA)) %>% 
  ungroup() %>%
  mutate(ratio=fct_relevel(ratio, ratiosFull)) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logFamilyAbundance, color=line, group=line)) +
  geom_point(aes(x=ratio, y=logFamilyAbundance, color=line), alpha=0.5) +
  geom_hline(yintercept=c(0, -1, -2, -3, -4), linetype="dotted") +
  xlab("Ratio") +
  ylab("Relative abundance (log)") +
  ggtitle("Marinifilaceae-p3") +
  ylim(c(-4.1, 0)) +
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
        strip.text.y=element_blank())
save_plot(paste0(plotdir, "/Marinifilaceae-P3_XDA-XFA.png"), 
          p_Marinifilaceae_XDAXFA, nrow=1.5, ncol=2)


# Family trajectory plots, new filtering. ---------------------------------

plotdirHighAbundanceOneParent <- paste0(plotdir, "/trajectoryPlots-highAbundanceOneParent")

# Get list of families in each passage.
familyPassages <- theoretical_mixtures_HighAbundanceOneParent %>% 
  mutate(familyPassage=paste0(Family,".",passage))
familyPassages <- unique(familyPassages$familyPassage)

foreach(x=familyPassages) %do% {
  family=sub("\\..*","",x)
  passageNum=sub(".*\\.","",x)
  p_individualFamilyPlot <- theoretical_mixtures_HighAbundanceOneParent %>% 
    filter(Family==family & passage==passageNum) %>%
    ggplot() +
    geom_line(aes(x=ratio, y=logFamilyAbundance, color=line, group=line), alpha=0.6) +
    geom_point(aes(x=ratio, y=logFamilyAbundance, color=line, shape=familyPointType), alpha=0.8) +
    geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
    xlab("Ratio") +
    ylab("Relative abundance (log)") +
    ggtitle(paste0(family, "-p", passageNum)) +
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
    DEFAULTS.THEME_PRES +
    facet_wrap(~combo, scales="free") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
          strip.text.y=element_blank())
  plotName <- paste("family_abundances", family, sep="_")
  save_plot(paste0(plotdirHighAbundanceOneParent, "/", plotName, "-P", passageNum, ".png"), 
            p_individualFamilyPlot, nrow=1.5, ncol=2)
}

# Plot family abundances for each combo, faceted by OTU.
plotComboFamilyAbundancesHighAbundanceOneParent <- function(currentCombo, passageNum) {
  p_individualFamilyPlot <- theoretical_mixtures_HighAbundanceOneParent %>% 
    filter(combo==currentCombo & passage==passageNum) %>% 
    ggplot() +
    geom_line(aes(x=ratio, y=logFamilyAbundance, color=line, group=line), alpha=0.6) +
    geom_point(aes(x=ratio, y=logFamilyAbundance, color=line, shape=pointType), alpha=0.8) +
    geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
    xlab("Ratio") +
    ylab("Relative abundance (log)") +
    ggtitle(paste0(currentCombo, "-p", passageNum)) +
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
    DEFAULTS.THEME_PRES +
    facet_wrap(~Family, scales="free") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10),
          strip.text.y=element_blank())
  plotName <- paste("family_abundances", currentCombo, sep="_")
  save_plot(paste0(plotdirHighAbundanceOneParent, "/", plotName, "-P", passageNum, ".png"), 
            p_individualFamilyPlot, nrow=2, ncol=2.5)
}
sapply(combos, FUN=plotComboFamilyAbundancesHighAbundanceOneParent, 3)
sapply(combos, FUN=plotComboFamilyAbundancesHighAbundanceOneParent, 5)

# Family scatterplots, new filtering. ----------------------------------------------------

plotdirHighAbundanceOneParentScatter <- paste0(plotdir, "/scatterPlots-highAbundanceOneParent")

# Create abundance dataframe with theoretical mean in one column and actual mean in another.
family_sum_abundances <- full_join(
  t_sum %>% 
    filter(ratio %in% ratios) %>% 
    group_by(Family, passage, combo, ratio) %>% 
    slice(1) %>% 
    select(Family, passage, combo, ratio, familyAvgAbundance, familyRangeAbundance) %>% 
    rename(t_familyAvgAbundance=familyAvgAbundance, t_familyRangeAbundance=familyRangeAbundance),
  a_sum %>%
    group_by(Family, passage, combo, ratio) %>% 
    slice(1) %>% 
    select(Family, passage, combo, ratio, familyAvgAbundance, familyRangeAbundance) %>% 
    rename(a_familyAvgAbundance=familyAvgAbundance, a_familyRangeAbundance=familyRangeAbundance),
  by=c("Family", "passage", "combo", "ratio"))

# Function to plot family scatterplots for each family, faceted by combo.
plotFamilyScatterplotsHighAbundanceOneParent <- function(family, passageNum) {
  p_familyScatterplot <- family_sum_abundances %>%
    filter(Family==family & passage==passageNum) %>% 
    ggplot() +
    geom_point(aes(x=t_familyAvgAbundance, y=a_familyAvgAbundance, color=ratio)) +
    geom_errorbarh(aes(y=a_familyAvgAbundance, color=ratio, height=0,
                       xmin=t_familyAvgAbundance-t_familyRangeAbundance/2, 
                       xmax=t_familyAvgAbundance+t_familyRangeAbundance/2)) +
    geom_errorbar(aes(x=t_familyAvgAbundance, color=ratio,
                      ymin=a_familyAvgAbundance-a_familyRangeAbundance/2, 
                      ymax=a_familyAvgAbundance+a_familyRangeAbundance/2)) +
    scale_color_brewer(palette="Reds", name="Ratio") +
    geom_abline(aes(slope=1, intercept=0)) +
    geom_hline(aes(yintercept=-4), linetype="dotted") +
    geom_vline(aes(xintercept=-4), linetype="dotted") +
    ylab("Avg log actual relative abundance") +
    xlab("Avg log theoretical relative abundance") +
    ylim(c(-3,0)) +
    xlim(c(-3,0)) +
    ggtitle(paste0(family, "-p", passageNum)) +
    DEFAULTS.THEME_PRES +
    facet_wrap(~combo, scales="free")
  plotName <- paste("family_scatterplot", family, sep="_")
  save_plot(paste0(plotdirHighAbundanceOneParentScatter, "/", plotName, "-P", passageNum, ".png"),
            p_familyScatterplot, nrow=1.5, ncol=2)
}
# Plot individual abundance scatterplots for all OTUs.
sapply(familyList, FUN=plotFamilyScatterplotsHighAbundanceOneParent, 3)
sapply(familyList, FUN=plotFamilyScatterplotsHighAbundanceOneParent, 5)

# Function to plot family scatterplots for each combo, faceted by family.
plotComboFamilyScatterplotsHighAbundanceOneParent <- function(currentCombo, passageNum) {
  p_familyComboScatterplot <- family_sum_abundances %>%
    filter(combo==currentCombo & passage==passageNum) %>% 
    ggplot() +
    geom_point(aes(x=t_familyAvgAbundance, y=a_familyAvgAbundance, color=ratio)) +
    geom_errorbarh(aes(y=a_familyAvgAbundance, color=ratio, height=0,
                       xmin=t_familyAvgAbundance-t_familyRangeAbundance/2, 
                       xmax=t_familyAvgAbundance+t_familyRangeAbundance/2)) +
    geom_errorbar(aes(x=t_familyAvgAbundance, color=ratio,
                      ymin=a_familyAvgAbundance-a_familyRangeAbundance/2, 
                      ymax=a_familyAvgAbundance+a_familyRangeAbundance/2)) +
    scale_color_brewer(palette="Reds", name="Ratio") +
    geom_abline(aes(slope=1, intercept=0)) +
    geom_hline(aes(yintercept=-4), linetype="dotted") +
    geom_vline(aes(xintercept=-4), linetype="dotted") +
    ylab("Avg log actual relative abundance") +
    xlab("Avg log theoretical relative abundance") +
    ylim(c(-3,0)) +
    xlim(c(-3,0)) +
    ggtitle(paste0(currentCombo, "-p", passageNum)) +
    DEFAULTS.THEME_PRES +
    facet_wrap(~Family, scales="free")
  plotName <- paste("family_scatterplot", currentCombo, sep="_")
  save_plot(paste0(plotdirHighAbundanceOneParentScatter, "/", plotName, "-P", passageNum, ".png"),
            p_familyComboScatterplot, nrow=1.5, ncol=2)
}
# Plot individual abundance scatterplots for all OTUs.
sapply(combos, FUN=plotComboFamilyScatterplotsHighAbundanceOneParent, 3)
sapply(combos, FUN=plotComboFamilyScatterplotsHighAbundanceOneParent, 5)
