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

source("analysis/generateMixturesDataframe.R")
plotdir <- "analysis/summaryStatPlots"
plotdirUnfiltered <- paste0(plotdir, "/unfiltered")
plotdirHighAbundanceOneParent <- paste0(plotdir, "/highAbundanceOneParent")
plotdirParentDiff <- paste0(plotdir, "/parentDiff-noTheoretical")

# Reorder taxa by full taxonomy for correct color ordering.
theoretical_mixtures_df_categorized <- theoretical_mixtures_df_categorized %>% 
  left_join(taxaPalette %>% select(Family, taxa), by="Family")

# Difference and distance statistics. -------------------------------------------------

# Plot histograms of difference.
p_diffHist <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=avgdifference), binwidth=0.025) +
  ylab("Count") +
  xlab("Signed difference") +
  xlim(c(-2.5, 2.5)) +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "diffHist.png"), p_diffHist)

p_diffHistFacet <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=avgdifference), binwidth=0.025) +
  ylab("Count") +
  xlab("Signed difference") +
  xlim(c(-2.5, 2.5)) +
  facet_wrap(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "diffHistFacet.png"), p_diffHistFacet)

p_diffHistFiltered <- ASV_sumFiltered %>%
  ggplot() + 
  geom_histogram(aes(x=avgdifference), binwidth=0.025) +
  ylab("Count") +
  xlab("Signed diff") +
  xlim(c(-2.5, 2.5)) +
  ylim(c(0,30)) +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "diffHistFiltered.png"), p_diffHistFiltered)

p_diffHistFilteredFacet <- ASV_sumFiltered %>%
  ggplot() + 
  geom_histogram(aes(x=avgdifference), binwidth=0.025) +
  ylab("Count") +
  xlab("Signed difference") +
  xlim(c(-2.5, 2.5)) +
  facet_wrap(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "diffHistFilteredFacet.png"), p_diffHistFilteredFacet)

# Plot histogram of deviation.
p_distHist <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=avgdistance), binwidth=0.025) +
  ylab("Count") +
  xlab("Magnitude of difference") +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "distHist.png"), p_distHist)

p_distHistFacet <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=avgdistance), binwidth=0.025) +
  ylab("Count") +
  xlab("Magnitude of difference") +
  facet_wrap(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "distHistFacet.png"), p_distHistFacet)

p_distHistFiltered <- ASV_sumFiltered %>%
  ggplot() + 
  geom_histogram(aes(x=avgdistance), binwidth=0.025) +
  ylab("Count") +
  xlab("Magnitude") +
  ylim(c(0,30)) +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "distHistFiltered.png"), p_distHistFiltered)

p_distHistFilteredFacet <- ASV_sumFiltered %>%
  ggplot() + 
  geom_histogram(aes(x=avgdistance), binwidth=0.025) +
  ylab("Count") +
  xlab("Magnitude of difference") +
  facet_wrap(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "distHistFilteredFacet.png"), p_distHistFilteredFacet)

# Scatterplot of difference and deviation, passage 3.
p_diffDistScatterP3 <- ASV_sum %>%
  filter(passage==3) %>%
  ggplot() +
  geom_point(aes(x=avgdifference, y=avgdistance, color=Family), alpha=0.5) +
  xlab("Signed diff") +
  ylab("Magnitude") +
  xlim(c(-2.5, 2.5)) +
  ylim(c(0, 2.5)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "diffDistScatterP3.png"), p_diffDistScatterP3, nrow=1.5, ncol=1.75)

# Scatterplot of difference and deviation, passage 5.
p_diffDistScatterP5 <- ASV_sum %>%
  filter(passage==5) %>%
  ggplot() +
  geom_point(aes(x=avgdifference, y=avgdistance, color=Family), alpha=0.5) +
  xlab("Signed diff") +
  ylab("Magnitude") +
  xlim(c(-2.5, 2.5)) +
  ylim(c(0, 2.5)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "diffDistScatterP5.png"), p_diffDistScatterP5, nrow=1.5, ncol=1.75)

# Filtered scatterplot of difference and deviation, passage 3.
p_diffDistScatterFilteredP3 <- ASV_sumFiltered %>%
  filter(passage==3) %>%
  ggplot() +
  geom_point(aes(x=avgdifference, y=avgdistance, color=Family), alpha=0.5) +
  xlab("Signed diff") +
  ylab("Magnitude") +
  xlim(c(-2.5, 2.5)) +
  ylim(c(0, 2.5)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "diffDistScatterFilteredP3.png"), p_diffDistScatterFilteredP3, nrow=1.5, ncol=1.75)

# Filtered scatterplot of difference and deviation, passage 5.
p_diffDistScatterFilteredP5 <- ASV_sumFiltered %>%
  filter(passage==5) %>%
  ggplot() +
  geom_point(aes(x=avgdifference, y=avgdistance, color=Family), alpha=0.5) +
  xlab("Signed diff") +
  ylab("Magnitude") +
  xlim(c(-2.5, 2.5)) +
  ylim(c(0, 2.5)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "diffDistScatterFilteredP5.png"), p_diffDistScatterFilteredP5, nrow=1.5, ncol=1.75)

# Example scatterplot for one combo
p_diffDistScatterFilteredExample <- ASV_sumFiltered %>%
  filter(combo=="XDA-XDB") %>%
  ggplot() +
  geom_point(aes(x=avgdifference, y=avgdistance, color=Family), alpha=0.5) +
  xlab("Signed diff") +
  ylab("Magnitude") +
  xlim(c(-2.5, 2.5)) +
  ylim(c(0, 2.5)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "diffDistScatterFilteredExample.png"), 
          p_diffDistScatterFilteredExample, nrow=1, ncol=1)

# Example scatterplot 2 for one combo
p_diffDistScatterFilteredExample2 <- ASV_sumFiltered %>%
  filter(combo=="XDA-XFA") %>%
  ggplot() +
  geom_point(aes(x=avgdifference, y=avgdistance, color=Family), alpha=0.5) +
  xlab("Signed diff") +
  ylab("Magnitude") +
  xlim(c(-2.5, 2.5)) +
  ylim(c(0, 2.5)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "diffDistScatterFilteredExample2.png"), 
          p_diffDistScatterFilteredExample2, nrow=1, ncol=1)

# Boxplots.
p_sdBoxplot <- ASV_sum %>%
  pivot_longer(cols=c(t_sdAvg, a_sdAvg, avgdistance, avgdifference)) %>%
  rename(stat=name, val=value) %>%
  filter(stat %in% c("a_sdAvg", "t_sdAvg")) %>%
  mutate(stat=gsub("t_sdAvg", "Theoretical", stat),
         stat=gsub("a_sdAvg", "Actual", stat)) %>%
  ggplot() +
  geom_boxplot(aes(x=stat, y=val)) +
  geom_jitter(aes(x=stat, y=val, fill=factor(stat)), pch=21, alpha=0.1, width=0.1) +
  geom_hline(yintercept=c(0, 0.25, 0.5, 0.75, 1, 1.25), linetype="dotted") +
  scale_fill_manual(values=c("blue", "green")) +
  theme(legend.position="none") +
  ylab("Standard deviation") +
  xlab("Mixture type") +
  scale_y_continuous(limits=c(0, 1.25),
                     breaks=c(0, 0.25, 0.5, 0.75, 1, 1.25)) +
  facet_grid(~passage)
save_plot(paste0(plotdir, "/", "sdBoxplot.png"), p_sdBoxplot)

p_diffDistBoxplot <- ASV_sumFiltered %>%
  pivot_longer(cols=c(t_sdAvg, a_sdAvg, avgdistance, avgdifference)) %>%
  rename(stat=name, val=value) %>%
  filter(stat %in% c("avgdistance", "avgdifference")) %>%
  mutate(stat=gsub("avgdistance", "Magnitude", stat),
         stat=gsub("avgdifference", "Signed diff", stat),
         stat=fct_relevel(stat, levels=c("Signed diff", "Magnitude"))) %>%
  ggplot() +
  geom_boxplot(aes(x=stat, y=val)) +
  geom_jitter(aes(x=stat, y=val), fill="blue", pch=21, alpha=0.1, width=0.1) +
  geom_hline(yintercept=c(-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5), linetype="dotted") +
  ylab("Avg value") +
  xlab("Statistic") +
  scale_y_continuous(limits=c(-2.5, 2.5),
                     breaks=c(-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5)) +
  facet_grid(~passage)
save_plot(paste0(plotdir, "/", "diffDistBoxplot.png"), p_diffDistBoxplot)

# Plot range statistics. --------------------------------------------------

# Replace 0s is caused by only one replicate present.
ASV_sum_plotting <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(category %in% c("dose-dependent","overcolonizing","undercolonizing")) %>% 
  group_by(OTU, passage, combo) %>% 
  summarize(t_rangeAvg, a_rangeAvg) %>% 
  mutate(t_rangeAvg=replace(t_rangeAvg,t_rangeAvg==0, NA),
         a_rangeAvg=replace(a_rangeAvg,a_rangeAvg==0, NA)) %>% 
  unique()

# Plot theoretical range.
p_trangeHist <- ASV_sum_plotting %>%
  ggplot() + 
  geom_histogram(aes(x=t_rangeAvg), binwidth=0.025) +
  ylab("Count") +
  xlab("Theoretical range") +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "trangeHist.png"), p_trangeHist)

p_trangeHistFacet <- ASV_sum_plotting %>%
  ggplot() + 
  geom_histogram(aes(x=t_rangeAvg), binwidth=0.025) +
  ylab("Count") +
  xlab("Theoretical range") +
  facet_grid(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "trangeHistFacet.png"), p_trangeHistFacet)

# Plot actual range.
p_arangeHist <- ASV_sum_plotting %>%
  ggplot() + 
  geom_histogram(aes(x=a_rangeAvg), binwidth=0.025) +
  ylab("Count") +
  xlab("Actual range") +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "arangeHist.png"), p_arangeHist)

p_arangeHistFacet <- ASV_sum_plotting %>%
  ggplot() + 
  geom_histogram(aes(x=a_rangeAvg), binwidth=0.025) +
  ylab("Count") +
  xlab("Actual range") +
  facet_wrap(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdir, "/", "arangeHistFacet.png"), p_arangeHistFacet)

# Plot boxplot of theoretical and actual ranges.
p_rangeBoxplot <- ASV_sum_plotting %>%
  pivot_longer(cols=c(t_rangeAvg, a_rangeAvg)) %>% 
  rename(stat=name, val=value) %>%
  mutate(stat=gsub("t_rangeAvg", "Theoretical", stat),
         stat=gsub("a_rangeAvg", "Actual", stat)) %>%
  ggplot() +
  geom_boxplot(aes(x=stat, y=val)) +
  geom_jitter(aes(x=stat, y=val, fill=factor(stat)), pch=21, alpha=0.1, width=0.1) +
  geom_hline(yintercept=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5), linetype="dotted") +
  scale_fill_manual(values=c("blue", "green")) +
  theme(legend.position="none") +
  ylab("Range") +
  xlab("Mixture type") +
  facet_grid(~passage)
p_rangeBoxplot
save_plot(paste0(plotdir, "/", "rangeBoxplot.png"), p_rangeBoxplot)

# Difference and distance statistics with new filtering. -------------------

# Scatterplot of difference and magnitude, passage 3.
p_diffDistScatterP3 <- ASV_sum %>%
  filter(passage==3) %>%
  ggplot() +
  geom_point(aes(x=avgdifference, y=avgdistance, color=Family), alpha=0.5) +
  xlab("Signed diff") +
  ylab("Magnitude") +
  xlim(c(-1.5, 1.5)) +
  ylim(c(0, 1.5)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "diffDistScatterP3.png"), 
          p_diffDistScatterP3, nrow=1.5, ncol=1.75)

# Scatterplot of difference and magnitude, passage 5.
p_diffDistScatterP5 <- ASV_sum %>%
  filter(passage==5) %>%
  ggplot() +
  geom_point(aes(x=avgdifference, y=avgdistance, color=Family), alpha=0.5) +
  xlab("Signed diff") +
  ylab("Magnitude") +
  xlim(c(-1.5, 1.5)) +
  ylim(c(0, 1.5)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "diffDistScatterP5.png"), 
          p_diffDistScatterP5, nrow=1.5, ncol=1.75)

# Plot histograms of difference.
p_diffHist <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=avgdifference), binwidth=0.025) +
  ylab("Count") +
  xlab("Signed difference") +
  xlim(c(-1.5, 1.5)) +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "diffHist.png"), p_diffHist)

p_diffHistFacet <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=avgdifference), binwidth=0.025) +
  ylab("Count") +
  xlab("Signed difference") +
  xlim(c(-1.5, 1.5)) +
  facet_wrap(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "diffHistFacet.png"), p_diffHistFacet)

# Plot histograms of magnitude.
p_distHist <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=avgdistance), binwidth=0.025) +
  ylab("Count") +
  xlab("Magnitude of difference") +
  xlim(c(0,1.5)) +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "distHist.png"), p_distHist)

p_distHistFacet <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=avgdistance), binwidth=0.025) +
  ylab("Count") +
  xlab("Magnitude of difference") +
  xlim(c(0,1.5)) +
  facet_wrap(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "distHistFacet.png"), p_distHistFacet)

# Dose distance and difference statistics with new filtering. ------------

# Scatterplot of dose difference and dose distance, passage 3.
p_doseDiffDistScatterP3 <- ASV_sum %>%
  filter(passage==3) %>%
  ggplot() +
  geom_point(aes(x=a_doseDifference, y=a_doseDistance, color=Family), alpha=0.5) +
  xlab("Signed diff between doses") +
  ylab("Magnitude between doses") +
  xlim(c(-0.5, 0.5)) +
  ylim(c(0, 1.2)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "doseDiffDistScatterP3.png"), 
          p_doseDiffDistScatterP3, nrow=1.5, ncol=1.75)

# Scatterplot of dose difference and dose distance, passage 5.
p_doseDiffDistScatterP5 <- ASV_sum %>%
  filter(passage==5) %>%
  ggplot() +
  geom_point(aes(x=a_doseDifference, y=a_doseDistance, color=Family), alpha=0.5) +
  xlab("Signed diff between doses") +
  ylab("Magnitude between doses") +
  xlim(c(-0.5, 0.5)) +
  ylim(c(0, 1.2)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "doseDiffDistScatterP5.png"), 
          p_doseDiffDistScatterP5, nrow=1.5, ncol=1.75)

# Plot histograms of dose difference.
p_doseDiffHist <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=a_doseDifference), binwidth=0.025) +
  ylab("Count") +
  xlab("Signed diff between doses") +
  xlim(c(-0.5, 0.5)) +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "doseDiffHist.png"), p_doseDiffHist)

p_doseDiffHistFacet <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=a_doseDifference), binwidth=0.025) +
  ylab("Count") +
  xlab("Signed diff between doses") +
  xlim(c(-0.5, 0.5)) +
  facet_wrap(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "doseDiffHistFacet.png"), p_doseDiffHistFacet)

# Plot histogram of deviation.
p_doseDistHist <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=a_doseDistance), binwidth=0.025) +
  ylab("Count") +
  xlab("Magnitude between doses") +
  xlim(c(0,1.2)) +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "doseDistHist.png"), p_doseDistHist)

p_doseDistHistFacet <- ASV_sum %>%
  ggplot() + 
  geom_histogram(aes(x=a_doseDistance), binwidth=0.025) +
  ylab("Count") +
  xlab("Magnitude between doses") +
  xlim(c(0,1.2)) +
  facet_wrap(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "doseDistHistFacet.png"), p_doseDistHistFacet)

# Difference and dose difference statistics with new filtering. -----------------------------------------

# Scatterplot of difference and dose difference, passage 3.
p_diffDoseDiffScatterP3 <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(passage==3) %>%
  ggplot() +
  geom_point(aes(x=avgdifference, y=a_doseDifference, color=taxa)) +
  scale_x_continuous(name="Difference from theoretical", 
                     breaks=c(-1,-0.5,0,0.5,1), limits=c(-1.25,1.25)) +
  scale_y_continuous(name="Signed change across mixture ratios",
                     breaks=c(-0.5,-0.25,0,0.25,0.5), limits=c(-0.5,0.5)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "diffDoseDiffScatterP3.png"), 
          p_diffDoseDiffScatterP3, nrow=1.5, ncol=1.75)

# Scatterplot of difference and dose difference, passage 5.
p_diffDoseDiffScatterP5 <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(category %in% c("overcolonizing","undercolonizing","dose-dependent") & passage==5) %>%
  ungroup() %>% 
  mutate(category=fct_relevel(category, c("overcolonizing","undercolonizing","dose-dependent"))) %>% 
  ggplot() +
  geom_vline(aes(xintercept=-0.2), linetype="dashed", size=0.5) +
  geom_vline(aes(xintercept=0.2), linetype="dashed", size=0.5) +
  geom_hline(aes(yintercept=-0.155), linetype="dashed", size=0.5) +
  geom_hline(aes(yintercept=0.155), linetype="dashed", size=0.5) +
  geom_point(aes(x=avgdifference, y=a_doseDifference, color=category), size=0.5) +
  geom_segment(aes(x=-1.5, xend=-0.21, y=-0.145, yend=-0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.21, y=0.145, yend=0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=-0.21, xend=-0.21, y=-0.145, yend=0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=0.21, xend=1.5, y=-0.145, yend=-0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=0.21, xend=1.5, y=0.145, yend=0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=0.21, xend=0.21, y=-0.145, yend=0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.19, y=-0.165, yend=-0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.19, y=0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=-0.19, xend=-0.19, y=-0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=1.5, y=-0.165, yend=-0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=1.5, y=0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=0.19, y=-0.165, yend=0.165), color="#A63A43", size=0.25) +
  scale_x_continuous(name="Difference from theoretical expectation", 
                     breaks=c(-1.5,-0.75,0,0.75,1.5), limits=c(-1.5,1.5), 
                     labels=c("-1.5","-0.75","0","0.75","1.5")) +
  scale_y_continuous(name="Signed change across mixture ratios",
                     breaks=c(-0.5,-0.25,0,0.25,0.5), limits=c(-0.5,0.5), 
                     labels=c("-0.5","-0.25","0","0.25","0.5")) +
  scale_color_manual(name="Category", 
                     labels=c("overcolonizing"="Overcolonizing",
                              "undercolonizing"="Undercolonizing",
                              "dose-dependent"="Dose-dependent"),
                     values=c("overcolonizing"="#D5A021",
                              "undercolonizing"="#3C91E6",
                              "dose-dependent"="#A63A43")) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRINT
p_diffDoseDiffScatterP5
save_plot(paste0(plotdirParentDiff, "/", "diffDoseDiffScatterP5.pdf"), 
          p_diffDoseDiffScatterP5, base_width=7, base_height=3.5)

# Unfaceted scatterplot of difference and dose difference, passage 5.
p_diffDoseDiffScatterP5Combined <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(category %in% c("dose-dependent","undercolonizing", "overcolonizing") & passage==5) %>% 
  ggplot() +
  geom_vline(aes(xintercept=-0.2), linetype="dashed", size=0.5) +
  geom_vline(aes(xintercept=0.2), linetype="dashed", size=0.5) +
  geom_hline(aes(yintercept=-0.155), linetype="dashed", size=0.5) +
  geom_hline(aes(yintercept=0.155), linetype="dashed", size=0.5) +
  geom_point(aes(x=avgdifference, y=a_doseDifference, color=taxa), size=0.5) +
  geom_segment(aes(x=-1.5, xend=-0.21, y=-0.145, yend=-0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.21, y=0.145, yend=0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=-0.21, xend=-0.21, y=-0.145, yend=0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=0.21, xend=1.5, y=-0.145, yend=-0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=0.21, xend=1.5, y=0.145, yend=0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=0.21, xend=0.21, y=-0.145, yend=0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.19, y=-0.165, yend=-0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.19, y=0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=-0.19, xend=-0.19, y=-0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=1.5, y=-0.165, yend=-0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=1.5, y=0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=0.19, y=-0.165, yend=0.165), color="#A63A43", size=0.25) +
  scale_x_continuous(name="Difference from theoretical expectation", 
                     breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), limits=c(-1.5,1.5)) +
  scale_y_continuous(name="Signed change across \nmixture ratios",
                     breaks=c(-0.5,-0.25,0,0.25,0.5), limits=c(-0.5,0.5)) +
  scale_color_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_diffDoseDiffScatterP5Combined
save_plot(paste0(plotdirParentDiff, "/", "diffDoseDiffScatterP5Combined.pdf"), 
          p_diffDoseDiffScatterP5Combined, base_width=2.75, base_height=1.75)

# Scatterplot of difference and dose difference for all categories, passage 5.
p_diffDoseDiffScatterP5All <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(passage==5) %>%
  group_by(combo, OTU) %>% 
  summarize(avgdifference, a_doseDifference, category) %>% 
  unique() %>% 
  ungroup() %>% 
  mutate(category=fct_relevel(category, c("noisy","flat","overcolonizing","undercolonizing","dose-dependent"))) %>% 
  arrange(category) %>% 
  ggplot() +
  geom_vline(aes(xintercept=-0.2), linetype="dashed", size=0.5) +
  geom_vline(aes(xintercept=0.2), linetype="dashed", size=0.5) +
  geom_hline(aes(yintercept=-0.155), linetype="dashed", size=0.5) +
  geom_hline(aes(yintercept=0.155), linetype="dashed", size=0.5) +
  geom_point(aes(x=avgdifference, y=a_doseDifference, color=category), size=0.5, alpha=0.75) +
  scale_x_continuous(name="Difference from theoretical", 
                     breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), limits=c(-1.5,1.5)) +
  scale_y_continuous(name="Signed change across \nmixture ratios",
                     breaks=c(-0.5,-0.25,0,0.25,0.5), limits=c(-0.5,0.5)) +
  scale_color_manual(name="Category", 
                    labels=c("noisy"="Noisy",
                             "flat"="Flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("noisy"="#CECECE",
                             "flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#3C91E6",
                             "dose-dependent"="#A63A43")) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  DEFAULTS.THEME_PRINT
p_diffDoseDiffScatterP5All
save_plot(paste0(plotdirParentDiff, "/diffDoseDiffScatterP5All.pdf"), 
          p_diffDoseDiffScatterP5All, base_width=3.5, base_height=2.5)

# Scatterplot of difference and dose difference for flat ASVs, passage 5.
p_diffDoseDiffScatterP5flat <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(passage==5 & category=="flat") %>%
  mutate(category="Flat") %>% 
  group_by(combo, OTU) %>% 
  summarize(avgdifference, a_doseDifference, category) %>% 
  unique() %>% 
  ggplot() +
  geom_vline(aes(xintercept=-0.2), linetype="dashed", size=0.5) +
  geom_vline(aes(xintercept=0.2), linetype="dashed", size=0.5) +
  geom_hline(aes(yintercept=-0.155), linetype="dashed", size=0.5) +
  geom_hline(aes(yintercept=0.155), linetype="dashed", size=0.5) +
  geom_point(aes(x=avgdifference, y=a_doseDifference), color="#616667", size=0.5, alpha=0.75) +
  geom_segment(aes(x=-1.5, xend=-0.21, y=-0.145, yend=-0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.21, y=0.145, yend=0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=-0.21, xend=-0.21, y=-0.145, yend=0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=0.21, xend=1.5, y=-0.145, yend=-0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=0.21, xend=1.5, y=0.145, yend=0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=0.21, xend=0.21, y=-0.145, yend=0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.19, y=-0.165, yend=-0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.19, y=0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=-0.19, xend=-0.19, y=-0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=1.5, y=-0.165, yend=-0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=1.5, y=0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=0.19, y=-0.165, yend=0.165), color="#A63A43", size=0.25) +
  scale_x_continuous(name="Difference from theoretical expectation", 
                     breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), limits=c(-1.5,1.5)) +
  scale_y_continuous(name="Signed change across \nmixture ratios",
                     breaks=c(-0.5,-0.25,0,0.25,0.5), limits=c(-0.5,0.5)) +
  DEFAULTS.THEME_PRINT
p_diffDoseDiffScatterP5flat
save_plot(paste0(plotdirParentDiff, "/diffDoseDiffScatterP5flat.pdf"), 
          p_diffDoseDiffScatterP5flat, base_width=2.75, base_height=1.75)

# Scatterplot of difference and dose difference for noisy ASVs, passage 5.
p_diffDoseDiffScatterP5noisy <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(passage==5 & category=="noisy") %>%
  mutate(category="Noisy") %>% 
  group_by(combo, OTU) %>% 
  summarize(avgdifference, a_doseDifference, category) %>% 
  unique() %>% 
  ggplot() +
  geom_vline(aes(xintercept=-0.2), linetype="dashed", size=0.5) +
  geom_vline(aes(xintercept=0.2), linetype="dashed", size=0.5) +
  geom_hline(aes(yintercept=-0.155), linetype="dashed", size=0.5) +
  geom_hline(aes(yintercept=0.155), linetype="dashed", size=0.5) +
  geom_point(aes(x=avgdifference, y=a_doseDifference), color="#999999", size=0.5, alpha=0.75) +
  geom_segment(aes(x=-1.5, xend=-0.21, y=-0.145, yend=-0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.21, y=0.145, yend=0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=-0.21, xend=-0.21, y=-0.145, yend=0.145), color="#3C91E6", size=0.25) +
  geom_segment(aes(x=0.21, xend=1.5, y=-0.145, yend=-0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=0.21, xend=1.5, y=0.145, yend=0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=0.21, xend=0.21, y=-0.145, yend=0.145), color="#D5A021", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.19, y=-0.165, yend=-0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=-1.5, xend=-0.19, y=0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=-0.19, xend=-0.19, y=-0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=1.5, y=-0.165, yend=-0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=1.5, y=0.165, yend=0.165), color="#A63A43", size=0.25) +
  geom_segment(aes(x=0.19, xend=0.19, y=-0.165, yend=0.165), color="#A63A43", size=0.25) +
  scale_x_continuous(name="Difference from theoretical expectation", 
                     breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), limits=c(-1.5,1.5)) +
  scale_y_continuous(name="Signed change across \nmixture ratios",
                     breaks=c(-0.5,-0.25,0,0.25,0.5), limits=c(-0.5,0.5)) +
  DEFAULTS.THEME_PRINT
p_diffDoseDiffScatterP5noisy
save_plot(paste0(plotdirParentDiff, "/diffDoseDiffScatterP5noisy.pdf"), 
          p_diffDoseDiffScatterP5noisy, base_width=2.75, base_height=1.75)

# Plot for paper figure.
save_plot(paste0(plotdirParentDiff, "/diffDoseDiffScatter_all.pdf"),
          p_diffDoseDiffScatterP5Combined + p_diffDoseDiffScatterP5All,
          base_width=6, base_height=2)

# Scatterplot of difference and dose difference for all categories, each combo.
foreach(x=combos) %do% {
  p_diffDoseDiffScatterP5All_combo <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
    filter(passage==5) %>%
    group_by(combo, OTU) %>% 
    summarize(avgdifference, a_doseDifference, category) %>% 
    unique() %>% 
    ungroup() %>% 
    mutate(category=fct_relevel(category, c("noisy","flat","overcolonizing","undercolonizing","dose-dependent"))) %>% 
    filter(combo==x) %>% 
    arrange(category) %>% 
    ggplot() +
    geom_vline(aes(xintercept=-0.2), linetype="dashed", size=0.5) +
    geom_vline(aes(xintercept=0.2), linetype="dashed", size=0.5) +
    geom_hline(aes(yintercept=-0.155), linetype="dashed", size=0.5) +
    geom_hline(aes(yintercept=0.155), linetype="dashed", size=0.5) +
    geom_point(aes(x=avgdifference, y=a_doseDifference, color=category), size=0.5, alpha=0.75) +
    geom_segment(aes(x=-1.5, xend=-0.21, y=-0.145, yend=-0.145), color="#3C91E6", size=0.25) +
    geom_segment(aes(x=-1.5, xend=-0.21, y=0.145, yend=0.145), color="#3C91E6", size=0.25) +
    geom_segment(aes(x=-0.21, xend=-0.21, y=-0.145, yend=0.145), color="#3C91E6", size=0.25) +
    geom_segment(aes(x=0.21, xend=1.5, y=-0.145, yend=-0.145), color="#D5A021", size=0.25) +
    geom_segment(aes(x=0.21, xend=1.5, y=0.145, yend=0.145), color="#D5A021", size=0.25) +
    geom_segment(aes(x=0.21, xend=0.21, y=-0.145, yend=0.145), color="#D5A021", size=0.25) +
    geom_segment(aes(x=-1.5, xend=-0.19, y=-0.165, yend=-0.165), color="#A63A43", size=0.25) +
    geom_segment(aes(x=-1.5, xend=-0.19, y=0.165, yend=0.165), color="#A63A43", size=0.25) +
    geom_segment(aes(x=-0.19, xend=-0.19, y=-0.165, yend=0.165), color="#A63A43", size=0.25) +
    geom_segment(aes(x=0.19, xend=1.5, y=-0.165, yend=-0.165), color="#A63A43", size=0.25) +
    geom_segment(aes(x=0.19, xend=1.5, y=0.165, yend=0.165), color="#A63A43", size=0.25) +
    geom_segment(aes(x=0.19, xend=0.19, y=-0.165, yend=0.165), color="#A63A43", size=0.25) +
    scale_x_continuous(name="Difference from theoretical", 
                       breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), limits=c(-1.5,1.5)) +
    scale_y_continuous(name="Signed change across \nmixture ratios",
                       breaks=c(-0.5,-0.25,0,0.25,0.5), limits=c(-0.5,0.5)) +
    scale_color_manual(name="Category", 
                       labels=c("noisy"="Noisy",
                                "flat"="Flat",
                                "overcolonizing"="Overcolonizing",
                                "undercolonizing"="Undercolonizing",
                                "dose-dependent"="Dose-dependent"),
                       values=c("noisy"="#999999",
                                "flat"="#616667",
                                "overcolonizing"="#D5A021",
                                "undercolonizing"="#3C91E6",
                                "dose-dependent"="#A63A43")) +
    theme(legend.title=element_blank(),
          legend.key.size=unit(0.5,"lines")) +
    DEFAULTS.THEME_PRINT
  save_plot(paste0(plotdirParentDiff, "/diffDoseDiffScatterP5All_",x,".png"),
            p_diffDoseDiffScatterP5All_combo, base_width=3.5, base_height=2.5)
}

# Scatterplot of family difference and dose difference, passage 3.
p_diffDoseDiffScatterP3_family <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(passage==3) %>%
  ggplot() +
  geom_point(aes(x=familyAvgDifference, y=a_familyDoseDifference, color=taxa)) +
  scale_x_continuous(name="Difference from theoretical", 
                     breaks=c(-1,-0.5,0,0.5,1), limits=c(-1.25,1.25)) +
  scale_y_continuous(name="Signed change across mixture ratios",
                     breaks=c(-0.5,-0.25,0,0.25,0.5), limits=c(-0.5,0.5)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "diffDoseDiffScatterP3_family.png"), 
          p_diffDoseDiffScatterP3_family, nrow=1.5, ncol=1.75)

# Scatterplot of difference and dose difference, passage 5.
p_diffDoseDiffScatterP5_family <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(passage==5) %>%
  ggplot() +
  geom_point(aes(x=familyAvgDifference, y=a_familyDoseDifference, color=taxa)) +
  scale_x_continuous(name="Difference from theoretical", 
                     breaks=c(-1,-0.5,0,0.5,1), limits=c(-1.25,1.25)) +
  scale_y_continuous(name="Signed change across mixture ratios",
                     breaks=c(-0.5,-0.25,0,0.25,0.5), limits=c(-0.5,0.5)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "diffDoseDiffScatterP5_family.png"), 
          p_diffDoseDiffScatterP5_family, nrow=1.5, ncol=1.75)

# Scatterplot of difference and dose magnitude, passage 3.
p_diffDoseDistScatterP3 <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(passage==3) %>%
  ggplot() +
  geom_point(aes(x=avgdifference, y=a_doseDistance, color=taxa)) +
  scale_x_continuous(name="Difference from theoretical", 
                     breaks=c(-1,-0.5,0,0.5,1), limits=c(-1.25,1.25)) +
  scale_y_continuous(name="Magnitude change across mixture ratios",
                     breaks=c(0.25,0.5,0.75,1), limits=c(0,1.1)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "diffDoseDistScatterP3.png"), 
          p_diffDoseDistScatterP3, nrow=1.5, ncol=1.75)

# Scatterplot of difference and dose magnitude, passage 5.
p_diffDoseDistScatterP5 <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(passage==5) %>%
  ggplot() +
  geom_point(aes(x=avgdifference, y=a_doseDistance, color=taxa)) +
  scale_x_continuous(name="Difference from theoretical", 
                     breaks=c(-1,-0.5,0,0.5,1), limits=c(-1.25,1.25)) +
  scale_y_continuous(name="Magnitude change across mixture ratios",
                     breaks=c(0.25,0.5,0.75,1), limits=c(0,1.1)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "diffDoseDistScatterP5.png"), 
          p_diffDoseDistScatterP5, nrow=1.5, ncol=1.75)

# Scatterplot of family difference and dose magnitude, passage 3.
p_diffDoseDistScatterP3_family <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(passage==3) %>%
  ggplot() +
  geom_point(aes(x=familyAvgDifference, y=a_familyDoseDistance, color=taxa)) +
  scale_x_continuous(name="Difference from theoretical", 
                     breaks=c(-1,-0.5,0,0.5,1), limits=c(-1.25,1.25)) +
  scale_y_continuous(name="Magnitude change across mixture ratios",
                     breaks=c(0.25,0.5,0.75,1), limits=c(0,1.1)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "diffDoseDistScatterP3_family.png"), 
          p_diffDoseDistScatterP3_family, nrow=1.5, ncol=1.75)

# Scatterplot of difference and dose magnitude, passage 5.
p_diffDoseDistScatterP5_family <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  filter(passage==5) %>%
  ggplot() +
  geom_point(aes(x=familyAvgDifference, y=a_familyDoseDistance, color=taxa)) +
  scale_x_continuous(name="Difference from theoretical", 
                     breaks=c(-1,-0.5,0,0.5,1), limits=c(-1.25,1.25)) +
  scale_y_continuous(name="Magnitude change across mixture ratios",
                     breaks=c(0.25,0.5,0.75,1), limits=c(0,1.1)) +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=2, ncol=4, scales="free") +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "diffDoseDistScatterP5_family.png"), 
          p_diffDoseDistScatterP5_family, nrow=1.5, ncol=1.75)

# Histograms of difference and dose difference. ---------------------------

p_diffHistFlat <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
  filter(passage==5 & category=="flat") %>%
  mutate(category="Flat") %>% 
  group_by(combo, Family, OTUnum) %>% 
  summarize(category, avgdifference) %>% 
  unique() %>% 
  ggplot() +
  geom_histogram(aes(x=avgdifference), bins=80, fill="#616667") +
  scale_x_continuous(name="Difference from theoretical", 
                     breaks=c(-1.5,-0.75,0,0.75,1.5), limits=c(-1.5,1.5),
                     labels=c("-1.5","-0.75","0","0.75","1.5")) +
  ylab("Count") +
  ylim(c(0, 25)) +
  facet_wrap(~category) +
  DEFAULTS.THEME_PRINT
p_diffHistFlat

p_doseDiffHistFlat <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
  filter(passage==5 & category=="flat") %>% 
  mutate(category="Flat") %>%
  group_by(combo, Family, OTUnum) %>% 
  summarize(category, a_doseDifference) %>% 
  unique() %>% 
  ggplot() +
  geom_histogram(aes(x=a_doseDifference), bins=80, fill="#616667") +
  scale_x_continuous(name="Signed change \nacross mixture ratios",
                     breaks=c(-0.5,-0.25,0,0.25,0.5), limits=c(-0.5,0.5),
                     labels=c("-0.5","-0.25","0","0.25","0.5")) +
  ylab("Count") +
  ylim(c(0, 25)) +
  facet_wrap(~category) +
  DEFAULTS.THEME_PRINT
p_doseDiffHistFlat

save_plot(paste0(plotdirParentDiff, "/diffDoseDiffHistsFlat.pdf"),
          p_diffHistFlat + p_doseDiffHistFlat, base_width=3, base_height=1.75)

p_diffHistNoisy <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
  filter(passage==5 & category=="noisy") %>% 
  mutate(category="Noisy") %>% 
  group_by(combo, Family, OTUnum) %>% 
  summarize(category, avgdifference) %>% 
  unique() %>% 
  ggplot() +
  geom_histogram(aes(x=avgdifference), bins=80, fill="#999999") +
  scale_x_continuous(name="Difference from theoretical", 
                     breaks=c(-1.5,-0.75,0,0.75,1.5), limits=c(-1.5,1.5),
                     labels=c("-1.5","-0.75","0","0.75","1.5")) +
  ylab("Count") +
  ylim(c(0, 25)) +
  facet_wrap(~category) +
  DEFAULTS.THEME_PRINT
p_diffHistNoisy

p_doseDiffHistNoisy <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
  filter(passage==5 & category=="noisy") %>% 
  mutate(category="Noisy") %>% 
  group_by(combo, Family, OTUnum) %>% 
  summarize(category, a_doseDifference) %>% 
  unique() %>% 
  ggplot() +
  geom_histogram(aes(x=a_doseDifference), bins=80, fill="#999999") +
  scale_x_continuous(name="Signed change \nacross mixture ratios",
                     breaks=c(-0.5,-0.25,0,0.25,0.5), limits=c(-0.5,0.5),
                     labels=c("-0.5","-0.25","0","0.25","0.5")) +
  ylab("Count") +
  ylim(c(0, 25)) +
  facet_wrap(~category) +
  DEFAULTS.THEME_PRINT
p_doseDiffHistNoisy

save_plot(paste0(plotdirParentDiff, "/diffDoseDiffHistsNoisy.pdf"),
          p_diffHistNoisy + p_doseDiffHistNoisy, base_width=3, base_height=1.75)

# Range statistics with new filtering. -----------------------------------

# Replace 0s is caused by only one replicate present.
theoretical_mixtures_HighAbundanceOneParent_summaryFiltered_plotting <- 
  theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
  mutate(t_rangeAvg=replace(t_rangeAvg, t_rangeAvg==0, NA),
         a_rangeAvg=replace(a_rangeAvg, a_rangeAvg==0, NA)) %>% 
  group_by(passage, combo, Family, OTUnum) %>% 
  summarize(t_rangeAvg, a_rangeAvg) %>% 
  unique()

# Plot theoretical range.
p_trangeHist <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
  ggplot() + 
  geom_histogram(aes(x=t_rangeAvg), binwidth=0.05) +
  ylab("Count") +
  xlab("Theoretical range") +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
p_trangeHist
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "trangeHist.png"), p_trangeHist)

p_trangeHistFacet <- ASV_sum_plotting %>%
  ggplot() + 
  geom_histogram(aes(x=t_rangeAvg), binwidth=0.05) +
  ylab("Count") +
  xlab("Theoretical range") +
  facet_grid(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "trangeHistFacet.png"), p_trangeHistFacet)

# Plot actual range.
p_arangeHist <- ASV_sum_plotting %>%
  ggplot() + 
  geom_histogram(aes(x=a_rangeAvg), binwidth=0.05) +
  ylab("Count") +
  xlab("Actual range") +
  facet_grid(~passage) +
  DEFAULTS.THEME_PRES
p_arangeHist
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "arangeHist.png"), p_arangeHist)

p_arangeHistFacet <- ASV_sum_plotting %>%
  ggplot() + 
  geom_histogram(aes(x=a_rangeAvg), binwidth=0.05) +
  ylab("Count") +
  xlab("Actual range") +
  facet_wrap(~passage~combo) +
  DEFAULTS.THEME_PRES
save_plot(paste0(plotdirHighAbundanceOneParent, "/", "arangeHistFacet.png"), p_arangeHistFacet)

# Plot boxplot of theoretical and actual ranges.
p_rangeBoxplot <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered_plotting %>%
  pivot_longer(cols=c(t_rangeAvg, a_rangeAvg)) %>% 
  rename(stat=name, val=value) %>%
  mutate(stat=gsub("t_rangeAvg", "Theoretical", stat),
         stat=gsub("a_rangeAvg", "Actual", stat)) %>%
  ggplot() +
  geom_jitter(aes(x=stat, y=val, fill=factor(stat)), pch=21, width=0.1) +
  geom_boxplot(aes(x=stat, y=val), alpha=0) +
  geom_hline(yintercept=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5), linetype="dotted") +
  scale_fill_manual(values=c("orange", "blue")) +
  theme(legend.position="none") +
  ylab("Range") +
  xlab("Mixture type") +
  facet_grid(~passage) + 
  DEFAULTS.THEME_PRES
p_rangeBoxplot
save_plot(paste0(plotdirHighAbundanceOneParent, "/rangeBoxplot.png"), p_rangeBoxplot)

# Plot ASVs sorted by summary stats. --------------------------------------

# Calculate number of ASVs in each combo in each category.
ASVcategoriesNumber <- theoretical_mixtures_df_categorized %>%
  #filter(passage==5) %>% 
  group_by(combo, OTU, passage) %>% 
  summarize(category) %>% 
  unique() %>% 
  group_by(combo, category, passage) %>% 
  filter(!is.na(category)) %>%
  summarize(nOTU=n()) %>% 
  ungroup() %>% 
  mutate(category=fct_relevel(category, c("lowAbundance", "noisy","flat","overcolonizing","undercolonizing","dose-dependent")))

# Calculate number of ASVs in each category total.
ASVcategoriesNumberTotal <- ASVcategoriesNumber %>% 
  group_by(category, passage) %>% 
  summarize(nOTU=sum(nOTU))

# Plot number of ASVs in each combo in each category, passage 5.
p_ASVcategoriesNumberRBY <- ASVcategoriesNumber %>%
  filter(passage==5) %>% 
  #filter(category=="dose-dependent") %>% 
  filter(category %in% c("overcolonizing","undercolonizing","dose-dependent")) %>% 
  #filter(category %in% c("noisy","flat","overcolonizing","undercolonizing","dose-dependent")) %>%
  #mutate(combo = case_when(
  #  combo=="XBA-XBB" ~ "1",
  #  combo=="XCA-XCB" ~ "2",
  #  combo=="XDA-XDB" ~ "3",
  #  combo=="XFA-XFB" ~ "4",
  #  combo=="XBA-XCA" ~ "5",
  #  combo=="XCA-XDA" ~ "6",
  #  combo=="XDA-XFA" ~ "7",
  #  combo=="XFA-XBA" ~ "8",
  #)) %>% 
  #mutate(combo=fct_relevel(combo, c("1","2","3","4","5","6","7","8"))) %>% 
  mutate(combo = case_when(
   combo=="XBA-XBB" ~ "A1/A2",
   combo=="XCA-XCB" ~ "B1/B2",
   combo=="XDA-XDB" ~ "C1/C2",
   combo=="XFA-XFB" ~ "D1/D2",
   combo=="XBA-XCA" ~ "A1/B1",
   combo=="XCA-XDA" ~ "B1/C1",
   combo=="XDA-XFA" ~ "C1/D1",
   combo=="XFA-XBA" ~ "D1/A1",
  )) %>%
  mutate(combo=fct_relevel(combo, c("A1/A2","B1/B2","C1/C2","D1/D2","A1/B1","B1/C1","C1/D1","D1/A1"))) %>%
  # mutate(combo = case_when(
  #  combo=="XBA-XBB" ~ "XBA\nXBB",
  #  combo=="XCA-XCB" ~ "XCA\nXCB",
  #  combo=="XDA-XDB" ~ "XDA\nXDB",
  #  combo=="XFA-XFB" ~ "XFA\nXFB",
  #  combo=="XBA-XCA" ~ "XBA\nXCA",
  #  combo=="XCA-XDA" ~ "XCA\nXDA",
  #  combo=="XDA-XFA" ~ "XDA\nXFA",
  #  combo=="XFA-XBA" ~ "XFA\nXBA",
  # )) %>%
  # mutate(combo=fct_relevel(combo, c("XBA\nXBB","XCA\nXCB","XDA\nXDB","XFA\nXFB","XBA\nXCA","XCA\nXDA","XDA\nXFA","XFA\nXBA"))) %>% 
  ggplot() +
  #geom_bar(aes(x=nOTU, y=combo, fill=category), stat="identity", color="black", size=0.01) +
  geom_bar(aes(x=combo, y=nOTU, fill=category), stat="identity", color="black", size=0.1) +
  scale_fill_manual(name="Category",
                    labels=c(#"noisy"="Noisy",
                             #"flat"="Flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c(#"noisy"="#999999",
                             #"flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#81B8EF",
                             "dose-dependent"="#820402")) +
  #scale_y_discrete(name="Mixture") +
  #scale_x_continuous(name="# ASVs", breaks=c(0,20,40,60), limits=c(0,65)) +
  #xlab("# ASVs") +
  xlab("Mixture") +
  ylab("# ASVs (1:1)") +
  ylim(0,15) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
        #axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
  #DEFAULTS.THEME_PRES
p_ASVcategoriesNumberRBY
save_plot(paste0(plotdirParentDiff, "/ASVcategoriesNumberRBY.pdf"), p_ASVcategoriesNumberRBY, 
          base_width=3.4, base_height=1.2)
#save_plot(paste0(plotdirParentDiff, "/ASVcategoriesNumberRBY.png"), p_ASVcategoriesNumberRBY, 
#          base_width=3, base_height=1.75)

# Plot number of ASVs in each combo in each category, passage 3.
p_ASVcategoriesNumberRBYp3 <- ASVcategoriesNumber %>%
  filter(passage==3) %>%  
  filter(category %in% c("overcolonizing","undercolonizing","dose-dependent")) %>% 
  mutate(combo = case_when(
    combo=="XBA-XBB" ~ "A1/A2",
    combo=="XCA-XCB" ~ "B1/B2",
    combo=="XDA-XDB" ~ "C1/C2",
    combo=="XFA-XFB" ~ "D1/D2",
    combo=="XBA-XCA" ~ "A1/B1",
    combo=="XCA-XDA" ~ "B1/C1",
    combo=="XDA-XFA" ~ "C1/D1",
    combo=="XFA-XBA" ~ "D1/A1",
  )) %>%
    mutate(combo=fct_relevel(combo, c("A1/A2","B1/B2","C1/C2","D1/D2","A1/B1","B1/C1","C1/D1","D1/A1"))) %>%
  ggplot() +
    geom_bar(aes(x=combo, y=nOTU, fill=category), stat="identity", color="black", size=0.1) +
    scale_fill_manual(name="Category",
                      labels=c(#"noisy"="Noisy",
                        #"flat"="Flat",
                        "overcolonizing"="Overcolonizing",
                        "undercolonizing"="Undercolonizing",
                        "dose-dependent"="Dose-dependent"),
                      values=c(#"noisy"="#999999",
                        #"flat"="#616667",
                        "overcolonizing"="#D5A021",
                        "undercolonizing"="#81B8EF",
                        "dose-dependent"="#820402")) +
    xlab("Mixture") +
    ylab("# ASVs (1:1)") +
    #ylim(0,15) +
    theme(legend.title=element_blank(),
          legend.key.size=unit(0.5,"lines")) +
    theme(legend.position="none") +
    DEFAULTS.THEME_PRINT
p_ASVcategoriesNumberRBYp3
save_plot(paste0(plotdirParentDiff, "/ASVcategoriesNumberRBYp3.pdf"), p_ASVcategoriesNumberRBYp3, 
          base_width=2.75, base_height=1.75)
save_plot(paste0(plotdirParentDiff, "/ASVcategoriesNumberRBYp3.png"), p_ASVcategoriesNumberRBYp3, 
          base_width=3, base_height=1.75)

# Calculate relative abundance of ASVs in each combo in each category.
ASVcategoriesRelAb <- theoretical_mixtures_df_categorized %>% 
  filter((ratio %in% c("1:0","0:1") | (ratio=="1:1" & mixtureType=="actual")) & !is.na(rel_abundance_old)) %>%
  #filter(passage==5 & ratio %in% c("1000:1","1:1","1:1000") & mixtureType=="actual" & !is.na(rel_abundance_old)) %>% 
  group_by(passage, Family, OTUnum, combo, ratio) %>% 
  summarize(category, avgAbundance=ifelse(ratio=="1:1", sum(rel_abundance_old/3) , sum(rel_abundance_old/2))) %>%
  #summarize(category, avgAbundance = sum(rel_abundance_old/3)) %>% 
  unique() %>% 
  ungroup() %>%
  mutate(category=fct_relevel(category, c("lowAbundance","noisy","flat","overcolonizing","undercolonizing","dose-dependent"))) %>% 
  group_by(passage, category, combo, ratio) %>% 
  arrange(-avgAbundance)

# Plot relative abundance of ASVs in each combo in each category, p5.
p_ASVcategoriesRelAb <- ASVcategoriesRelAb %>% 
  filter(passage==5) %>% 
  unique() %>% 
  ungroup() %>% 
  # mutate(combo=case_when(
  #   combo=="XBA-XBB" ~ "XBA\nXBB",
  #   combo=="XCA-XCB" ~ "XCA\nXCB",
  #   combo=="XDA-XDB" ~ "XDA\nXDB",
  #   combo=="XFA-XFB" ~ "XFA\nXFB",
  #   combo=="XBA-XCA" ~ "XBA\nXCA",
  #   combo=="XCA-XDA" ~ "XCA\nXDA",
  #   combo=="XDA-XFA" ~ "XDA\nXFA",
  #   combo=="XFA-XBA" ~ "XFA\nXBA"
  # ),
  #combo = fct_relevel(combo, c("XBA\nXBB","XCA\nXCB","XDA\nXDB","XFA\nXFB","XBA\nXCA","XCA\nXDA","XDA\nXFA","XFA\nXBA"))) %>% 
  # combo = fct_relevel(combo, c("XFA\nXBA","XDA\nXFA","XCA\nXDA","XBA\nXCA","XFA\nXFB","XDA\nXDB","XCA\nXCB","XBA\nXBB"))) %>%
  mutate(combo = case_when(
    combo=="XBA-XBB" ~ "A1/A2",
    combo=="XCA-XCB" ~ "B1/B2",
    combo=="XDA-XDB" ~ "C1/C2",
    combo=="XFA-XFB" ~ "D1/D2",
    combo=="XBA-XCA" ~ "A1/B1",
    combo=="XCA-XDA" ~ "B1/C1",
    combo=="XDA-XFA" ~ "C1/D1",
    combo=="XFA-XBA" ~ "D1/A1",
  )) %>%
  mutate(combo=fct_relevel(combo, c("A1/A2","B1/B2","C1/C2","D1/D2","A1/B1","B1/C1","C1/D1","D1/A1"))) %>%
  filter(ratio=="1:1") %>% 
  #mutate(category=fct_relevel(category, c("overcolonizing","lowAbundance","undercolonizing",
  #                                        "noisy","dose-dependent","flat"))) %>% 
  ggplot() +
  geom_bar(aes(x=combo, y=avgAbundance, fill=category), stat="identity", color="black", size=0.1) +
  #geom_bar(aes(x=avgAbundance, y=combo, fill=category), stat="identity", color="black", size=0.1) +
  scale_fill_manual(name="Category", 
                    labels=c("lowAbundance"="Low abundance",
                             "noisy"="Noisy",
                             "flat"="Resident",
                             "overcolonizing"="Strong colonizer (S)",
                             "undercolonizing"="Weak colonizer (W)",
                             "dose-dependent"="Dose-dependent (DD)"),
                    values=c("lowAbundance"="#1F2020",
                             "noisy"="#999999",
                             "flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#81B8EF",
                             "dose-dependent"="#820402")) +
  #xlab("Community") +
  #xlab("Relative abundance (1:1)") +
  xlab("Mixture") +
  #ylab("Mixture") +
  ylab("Relative abundance") +
  DEFAULTS.THEME_PRINT +
  #DEFAULTS.THEME_PRES +
  #facet_grid(~combo, scales="free") +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  guides(fill = guide_legend(nrow=2)) +
  theme(legend.position="none")
p_ASVcategoriesRelAb
#save_plot(paste0(plotdirParentDiff, "/ASVcategoriesRelAbSideways.pdf"), p_ASVcategoriesRelAb,
#          base_width=4.5, base_height=2.5)
save_plot(paste0(plotdirParentDiff, "/ASVcategoriesRelAb.pdf"), p_ASVcategoriesRelAb,
          base_width=3.4, base_height=1.2)

categoriesLegend <- get_legend(p_ASVcategoriesRelAb)
save_plot(paste0(plotdirParentDiff, "/categoriesLegend.pdf"), categoriesLegend, base_width=3, base_height=0.4)

# Plot relative abundance of ASVs in each combo in each category, p3.
p_ASVcategoriesRelAbp3 <- ASVcategoriesRelAb %>% 
  filter(passage==3) %>% 
  unique() %>% 
  ungroup() %>%
  mutate(combo = case_when(
    combo=="XBA-XBB" ~ "A1/A2",
    combo=="XCA-XCB" ~ "B1/B2",
    combo=="XDA-XDB" ~ "C1/C2",
    combo=="XFA-XFB" ~ "D1/D2",
    combo=="XBA-XCA" ~ "A1/B1",
    combo=="XCA-XDA" ~ "B1/C1",
    combo=="XDA-XFA" ~ "C1/D1",
    combo=="XFA-XBA" ~ "D1/A1",
  )) %>%
    mutate(combo=fct_relevel(combo, c("A1/A2","B1/B2","C1/C2","D1/D2","A1/B1","B1/C1","C1/D1","D1/A1"))) %>%
    filter(ratio=="1:1") %>% 
    ggplot() +
    geom_bar(aes(x=combo, y=avgAbundance, fill=category), stat="identity", color="black", size=0.1) +
    #geom_bar(aes(x=avgAbundance, y=combo, fill=category), stat="identity", color="black", size=0.1) +
    scale_fill_manual(name="Category", 
                      labels=c("lowAbundance"="Low abundance",
                               "noisy"="Noisy",
                               "flat"="Flat",
                               "overcolonizing"="Overcolonizing",
                               "undercolonizing"="Undercolonizing",
                               "dose-dependent"="Dose-dependent"),
                      values=c("lowAbundance"="#1F2020",
                               "noisy"="#999999",
                               "flat"="#616667",
                               "overcolonizing"="#D5A021",
                               "undercolonizing"="#81B8EF",
                               "dose-dependent"="#820402")) +
    #xlab("Relative abundance (1:1)") +
    xlab("Mixture") +
    #ylab("Mixture") +
    ylab("Relative abundance") +
    DEFAULTS.THEME_PRINT +
    #DEFAULTS.THEME_PRES +
    #facet_grid(~combo, scales="free") +
    theme(legend.title=element_blank(),
          legend.key.size=unit(0.5,"lines")) +
    theme(legend.position="none")
p_ASVcategoriesRelAbp3
save_plot(paste0(plotdirParentDiff, "/ASVcategoriesRelAbp3.pdf"), p_ASVcategoriesRelAbp3,
          base_width=2.75, base_height=3)
save_plot(paste0(plotdirParentDiff, "/ASVcategoriesRelAbp3.png"), p_ASVcategoriesRelAbp3,
          base_width=2.75, base_height=3)

#save_plot(paste0(plotdirParentDiff, "/ASVcategoriesCombined.pdf"), 
#          p_ASVcategoriesNumber + p_ASVcategoriesRelAb + plot_layout(widths=c(1,3)), base_width=7, base_height=2)


# Plot number and rel ab comparisons between passages. --------------------

# Number of ASVs in each category.
p_passageCompNumber <- ASVcategoriesNumberTotal %>% 
  filter(category!="lowAbundance") %>% 
  mutate(category=fct_relevel(category, c("dose-dependent","undercolonizing","overcolonizing","flat","noisy"))) %>% 
  ggplot() +
  geom_bar(aes(x=as.character(passage), y=nOTU, fill=category), position="dodge", stat="identity") +
  ylab("# ASVs") +
  xlab("Passage") +
  scale_fill_manual(name="Category", 
                    labels=c("lowAbundance"="Low abundance",
                             "noisy"="Noisy",
                             "flat"="Flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("lowAbundance"="#1F2020",
                             "noisy"="#999999",
                             "flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#81B8EF",
                             "dose-dependent"="#820402")) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_passageCompNumber
save_plot(paste0(plotdirParentDiff, "/passageCompNumber.pdf"), p_passageCompNumber, base_width=2, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/passageCompNumber.png"), p_passageCompNumber, base_width=2, base_height=1.6)

# Number of ASVs in each category, line plot.
p_passageCompNumberLines <- ASVcategoriesNumberTotal %>% 
  filter(category!="lowAbundance") %>% 
  ggplot() +
  geom_line(aes(x=as.character(passage), y=nOTU, group=category, color=category), size=0.4) +
  geom_point(aes(x=as.character(passage), y=nOTU, color=category), size=0.5) +
  xlab("Passage") +
  ylab("# ASVs") +
  ylim(0,225) +
  scale_color_manual(name="Category", 
                    labels=c("lowAbundance"="Low abundance",
                             "noisy"="Noisy",
                             "flat"="Flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("lowAbundance"="#1F2020",
                             "noisy"="#999999",
                             "flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#81B8EF",
                             "dose-dependent"="#820402")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_passageCompNumberLines
save_plot(paste0(plotdirParentDiff, "/passageCompNumberLines.pdf"), p_passageCompNumberLines, base_width=1.75, base_height=1.75)

# Plot average dose-dependence of each ASV in community, line plot.
passageCompDD <- theoretical_mixtures_df_categorized %>%
  filter(category!="lowAbundance") %>% 
  group_by(passage, combo, OTU) %>% 
  summarize(avgDD = mean(a_doseDifference)) %>% 
  ungroup() %>% 
  group_by(passage) %>% 
  summarize(avgPassageDD = abs(mean(avgDD)))

# Total relative abundance in each category.
p_passageCompRelab <- ASVcategoriesRelAb %>% 
  filter(category!="lowAbundance" & ratio=="1:1") %>%
  group_by(passage, category) %>% 
  summarize(sumRelAb = sum(avgAbundance)) %>% 
  mutate(category=fct_relevel(category, c("dose-dependent","undercolonizing","overcolonizing","flat","noisy"))) %>% 
  ggplot() +
  geom_bar(aes(x=as.character(passage), y=sumRelAb, fill=category), position="dodge", stat="identity") +
  ylab("Total relative abundance (1:1)") +
  xlab("Passage") +
  scale_fill_manual(name="Category", 
                    labels=c("lowAbundance"="Low abundance",
                             "noisy"="Noisy",
                             "flat"="Flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("lowAbundance"="#1F2020",
                             "noisy"="#999999",
                             "flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#81B8EF",
                             "dose-dependent"="#820402")) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_passageCompRelab
save_plot(paste0(plotdirParentDiff, "/passageCompRelab.pdf"), p_passageCompRelab, base_width=2, base_height=1.6)
save_plot(paste0(plotdirParentDiff, "/passageCompRelab.png"), p_passageCompRelab, base_width=2, base_height=1.6)

# Proportion of ASVs that go from DD to not vs not to DD.
changingDDASVs <- theoretical_mixtures_df_categorized %>% 
  group_by(OTU, combo) %>% 
  filter("dose-dependent" %in% category & 
           ("lowAbundance" %in% category | "noisy" %in% category | "flat" %in% category |
              "overcolonizing" %in% category | "undercolonizing" %in% category)) %>% 
  mutate(p3 = passage==3 & category=="dose-dependent",
         p5 = passage==5 & category=="dose-dependent") %>% 
  summarize(direction = ifelse(TRUE %in% p3, "lost", "gained"))
  
# Plot # ASVs that lose vs gain dose-dependence.
p_changingDDASVs <- changingDDASVs %>% 
  group_by(direction) %>% 
  summarize(nOTUs = n()) %>%
  mutate(direction = ifelse(direction == "lost", "Lost", "Gained"),
         direction = fct_relevel(direction, c("Lost","Gained"))) %>% 
  ggplot() +
  geom_bar(aes(x=direction, y=nOTUs), stat="identity") +
  ylab("# ASVs") +
  xlab("Change in dose dependence") +
  DEFAULTS.THEME_PRINT
p_changingDDASVs
save_plot(paste0(plotdirParentDiff, "/changingDDASVs.pdf"), p_changingDDASVs, base_width=1.4, base_height=1.4)
save_plot(paste0(plotdirParentDiff, "/changingDDASVs.png"), p_changingDDASVs, base_width=1.4, base_height=1.4)

# Proportion of ASVs that go from OC to not vs not to OC.
changingOCASVs <- theoretical_mixtures_df_categorized %>% 
  group_by(OTU, combo) %>% 
  filter("overcolonizing" %in% category & 
           ("lowAbundance" %in% category | "noisy" %in% category | "flat" %in% category |
              "dose-dependent" %in% category | "undercolonizing" %in% category)) %>% 
  mutate(p3 = passage==3 & category=="overcolonizing",
         p5 = passage==5 & category=="overcolonizing") %>% 
  summarize(direction = ifelse(TRUE %in% p3, "lost", "gained"))

# Plot # ASVs that lose vs gain overcolonization.
p_changingOCASVs <- changingOCASVs %>% 
  group_by(direction) %>% 
  summarize(nOTUs = n()) %>%
  mutate(direction = ifelse(direction == "lost", "Lost", "Gained"),
         direction = fct_relevel(direction, c("Lost","Gained"))) %>% 
  ggplot() +
  geom_bar(aes(x=direction, y=nOTUs), stat="identity") +
  ylab("# ASVs") +
  xlab("Change in strong colonizers") +
  DEFAULTS.THEME_PRINT
p_changingOCASVs
save_plot(paste0(plotdirParentDiff, "/changingOCASVs.pdf"), p_changingOCASVs, base_width=1.4, base_height=1.4)
save_plot(paste0(plotdirParentDiff, "/changingOCASVs.png"), p_changingOCASVs, base_width=1.4, base_height=1.4)

# Proportion of ASVs that go from UC to not vs not to UC.
changingUCASVs <- theoretical_mixtures_df_categorized %>% 
  group_by(OTU, combo) %>% 
  filter("undercolonizing" %in% category & 
           ("lowAbundance" %in% category | "noisy" %in% category | "flat" %in% category |
              "dose-dependent" %in% category | "overcolonizing" %in% category)) %>% 
  mutate(p3 = passage==3 & category=="undercolonizing",
         p5 = passage==5 & category=="undercolonizing") %>% 
  summarize(direction = ifelse(TRUE %in% p3, "lost", "gained"))

# Plot # ASVs that lose vs gain undercolonization.
p_changingUCASVs <- changingUCASVs %>% 
  group_by(direction) %>% 
  summarize(nOTUs = n()) %>%
  mutate(direction = ifelse(direction == "lost", "Lost", "Gained"),
         direction = fct_relevel(direction, c("Lost","Gained"))) %>% 
  ggplot() +
  geom_bar(aes(x=direction, y=nOTUs), stat="identity") +
  ylab("# ASVs") +
  xlab("Change in weak colonizers") +
  DEFAULTS.THEME_PRINT
p_changingUCASVs
save_plot(paste0(plotdirParentDiff, "/changingUCASVs.pdf"), p_changingUCASVs, base_width=1.4, base_height=1.4)
save_plot(paste0(plotdirParentDiff, "/changingUCASVs.png"), p_changingUCASVs, base_width=1.4, base_height=1.4)

# Proportion of ASVs that go from flat to not vs not to flat.
changingflatASVs <- theoretical_mixtures_df_categorized %>% 
  group_by(OTU, combo) %>% 
  filter("flat" %in% category & 
           ("lowAbundance" %in% category | "noisy" %in% category | "undercolonizing" %in% category |
              "dose-dependent" %in% category | "overcolonizing" %in% category)) %>% 
  mutate(p3 = passage==3 & category=="flat",
         p5 = passage==5 & category=="flat") %>% 
  summarize(direction = ifelse(TRUE %in% p3, "lost", "gained"))

# Plot # ASVs that lose vs gain flat colonization.
p_changingflatASVs <- changingflatASVs %>% 
  group_by(direction) %>% 
  summarize(nOTUs = n()) %>%
  mutate(direction = ifelse(direction == "lost", "Lost", "Gained"),
         direction = fct_relevel(direction, c("Lost","Gained"))) %>% 
  ggplot() +
  geom_bar(aes(x=direction, y=nOTUs), stat="identity") +
  ylab("# ASVs") +
  xlab("Change in residents") +
  DEFAULTS.THEME_PRINT
p_changingflatASVs
save_plot(paste0(plotdirParentDiff, "/changingflatASVs.pdf"), p_changingflatASVs, base_width=1.4, base_height=1.4)
save_plot(paste0(plotdirParentDiff, "/changingflatASVs.png"), p_changingflatASVs, base_width=1.4, base_height=1.4)

# Proportion of ASVs that go from noisy to not vs not to noisy.
changingnoisyASVs <- theoretical_mixtures_df_categorized %>% 
  group_by(OTU, combo) %>% 
  filter("noisy" %in% category & 
           ("lowAbundance" %in% category | "flat" %in% category | "undercolonizing" %in% category |
              "dose-dependent" %in% category | "overcolonizing" %in% category)) %>% 
  mutate(p3 = passage==3 & category=="noisy",
         p5 = passage==5 & category=="noisy") %>% 
  summarize(direction = ifelse(TRUE %in% p3, "lost", "gained"))

# Plot # ASVs that lose vs gain noisy colonization.
p_changingnoisyASVs <- changingnoisyASVs %>% 
  group_by(direction) %>% 
  summarize(nOTUs = n()) %>%
  mutate(direction = ifelse(direction == "lost", "Lost", "Gained"),
         direction = fct_relevel(direction, c("Lost","Gained"))) %>% 
  ggplot() +
  geom_bar(aes(x=direction, y=nOTUs), stat="identity") +
  ylab("# ASVs") +
  xlab("Change in noisy ASVs") +
  DEFAULTS.THEME_PRINT
p_changingnoisyASVs
save_plot(paste0(plotdirParentDiff, "/changingnoisyASVs.pdf"), p_changingnoisyASVs, base_width=1.4, base_height=1.4)
save_plot(paste0(plotdirParentDiff, "/changingnoisyASVs.png"), p_changingnoisyASVs, base_width=1.4, base_height=1.4)

# Plot all changing ASVs on same plot.
changingASVsAll <- rbind(
  changingDDASVs %>% mutate(category = "DD"),
  changingOCASVs %>% mutate(category = "Strong"),
  changingUCASVs %>% mutate(category = "Weak"),
  changingflatASVs %>% mutate(category = "Resident"),
  changingnoisyASVs %>% mutate(category = "Noisy")
  ) %>% 
  group_by(category, direction) %>% 
  summarize(count = n()) %>% 
  mutate(count = ifelse(direction=="lost", count*(-1), count)) %>% 
  ungroup()

# Calculate changing ASVs as a proportion of p3 ASVs in each category.
changingASVsAll <- changingASVsAll %>% 
  mutate(proportion = case_when(
    category == "DD" ~ count / ASVcategoriesNumberTotal %>%
      filter(passage==3 & category=="dose-dependent") %>% pull(nOTU),
    category == "Strong" ~ count / ASVcategoriesNumberTotal %>%
      filter(passage==3 & category=="overcolonizing") %>% pull(nOTU),
    category == "Weak" ~ count / ASVcategoriesNumberTotal %>%
      filter(passage==3 & category=="undercolonizing") %>% pull(nOTU),
    category == "Resident" ~ count / ASVcategoriesNumberTotal %>%
      filter(passage==3 & category=="flat") %>% pull(nOTU),
    category == "Noisy" ~ count / ASVcategoriesNumberTotal %>%
      filter(passage==3 & category=="noisy") %>% pull(nOTU),
  ))

p_changingASVsAll <- changingASVsAll %>% 
  mutate(category = fct_relevel(category, c("DD","Strong","Weak","Resident","Noisy"))) %>% 
  ggplot() +
  geom_bar(aes(x=category, y=count, fill=category), stat="identity") +
  geom_hline(yintercept=0) +
  xlab("Category") +
  ylab("# ASVs gained/lost from passage 3 to 5") +
  ylim(-140,140) +
  scale_fill_manual(name="Category",
                    values=c("Noisy"="#999999",
                             "Resident"="#616667",
                             "Strong"="#D5A021",
                             "Weak"="#81B8EF",
                             "DD"="#820402")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_changingASVsAll
save_plot(paste0(plotdirParentDiff, "/changingASVsAll.pdf"), p_changingASVsAll, base_width=2.5, base_height=1.6)

p_changingASVsAllProportion <- changingASVsAll %>% 
  mutate(category = fct_relevel(category, c("DD","Strong","Weak","Resident","Noisy"))) %>% 
  ggplot() +
  geom_bar(aes(x=category, y=proportion, fill=category), stat="identity") +
  geom_hline(yintercept=0) +
  xlab("Category") +
  ylab("% p3 ASVs gained/lost from p 3 to p5") +
  #ylim(-140,140) +
  scale_fill_manual(name="Category",
                    values=c("Noisy"="#999999",
                             "Resident"="#616667",
                             "Strong"="#D5A021",
                             "Weak"="#81B8EF",
                             "DD"="#820402")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_changingASVsAllProportion
save_plot(paste0(plotdirParentDiff, "/changingASVsAllProportion.pdf"), p_changingASVsAllProportion, base_width=2.5, base_height=1.6)

# Identify categories that DD ASVs change to.
p3DDASVs <- theoretical_mixtures_df_categorized %>% 
  filter(passage==3 & category=="dose-dependent") %>%
  mutate(OTUcombo = paste0(OTU, "_", combo)) %>% 
  pull(OTUcombo) %>% 
  unique()

trackDDASVs <- theoretical_mixtures_df_categorized %>% 
  mutate(OTUcombo = paste0(OTU, "_", combo)) %>% 
  filter(OTUcombo %in% p3DDASVs & passage==5) %>% 
  group_by(category) %>% 
  summarize(count = n_distinct(OTUcombo))

# Identify categories that OC ASVs change to.
p3OCASVs <- theoretical_mixtures_df_categorized %>% 
  filter(passage==3 & category=="overcolonizing") %>%
  mutate(OTUcombo = paste0(OTU, "_", combo)) %>% 
  pull(OTUcombo) %>% 
  unique()

trackOCASVs <- theoretical_mixtures_df_categorized %>% 
  mutate(OTUcombo = paste0(OTU, "_", combo)) %>% 
  filter(OTUcombo %in% p3OCASVs & passage==5) %>% 
  group_by(category) %>% 
  summarize(count = n_distinct(OTUcombo))

# Identify categories that UC ASVs change to.
p3UCASVs <- theoretical_mixtures_df_categorized %>% 
  filter(passage==3 & category=="undercolonizing") %>%
  mutate(OTUcombo = paste0(OTU, "_", combo)) %>% 
  pull(OTUcombo) %>% 
  unique()

trackUCASVs <- theoretical_mixtures_df_categorized %>% 
  mutate(OTUcombo = paste0(OTU, "_", combo)) %>% 
  filter(OTUcombo %in% p3UCASVs & passage==5) %>% 
  group_by(category) %>% 
  summarize(count = n_distinct(OTUcombo))

# Plot ASVs sorted by summary stats, new sum categories. ------------------

# Calculate number of ASVs in each combo in each category.
ASVcategoriesSumNumber <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
  filter(passage==5) %>% 
  group_by(combo, OTU) %>% 
  summarize(categorySum) %>% 
  unique() %>% 
  group_by(combo, categorySum) %>% 
  filter(!is.na(categorySum)) %>%
  summarize(nOTU=n()) %>% 
  ungroup() %>% 
  mutate(categorySum=fct_relevel(categorySum, c("lowAbundance", "noisy","flat","actualFlat","overcolonizing","undercolonizing","dose-dependent")))

# Calculate number of ASVs in each category total.
ASVcategoriesSumNumberTotal <- ASVcategoriesSumNumber %>% 
  group_by(categorySum) %>% 
  summarize(nOTU=sum(nOTU))

# Plot number of ASVs in each combo in each category.
p_ASVcategoriesSumNumber <- ASVcategoriesSumNumber %>% 
  filter(categorySum %in% c("noisy","flat","overcolonizing","undercolonizing","dose-dependent","actualFlat")) %>% 
  ggplot() +
  geom_bar(aes(x=nOTU, y=combo, fill=categorySum), stat="identity", color="black", size=0.01) +
  scale_fill_manual(name="Category",
                    labels=c("noisy"="Noisy",
                             "flat"="Flat",
                             "actualFlat"="Actual flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("noisy"="#999999",
                             "flat"="#616667",
                             "actualFlat"="#120309",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#3C91E6",
                             "dose-dependent"="#A63A43")) +
  scale_y_discrete(name="Mixture") +
  scale_x_continuous(name="# ASVs", breaks=c(0,20,40,60), limits=c(0,65)) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  DEFAULTS.THEME_PRINT
p_ASVcategoriesSumNumber

# Calculate relative abundance of ASVs in each combo in each category.
ASVcategoriesSumRelAb <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
  filter(passage==5 & 
           (ratio %in% c("1:0","0:1") | (ratio=="1:1" & mixtureType=="actual")) & !is.na(rel_abundance_old)) %>%
  group_by(Family, OTUnum, combo, ratio) %>% 
  summarize(categorySum, avgAbundance=ifelse(ratio=="1:1", sum(rel_abundance_old/3) , sum(rel_abundance_old/2))) %>%
  unique() %>% 
  ungroup() %>%
  mutate(categorySum=fct_relevel(categorySum, c("lowAbundance","noisy","flat","actualFlat","overcolonizing","undercolonizing","dose-dependent"))) %>% 
  group_by(categorySum, combo, ratio) %>% 
  arrange(-avgAbundance)

# Plot relative abundance of ASVs in each combo in each category.
p_ASVcategoriesSumRelAb <- ASVcategoriesSumRelAb %>% 
  unique() %>% 
  mutate(combo=case_when(
    combo=="XBA-XBB" ~ "XBA-\nXBB",
    combo=="XCA-XCB" ~ "XCA-\nXCB",
    combo=="XDA-XDB" ~ "XDA-\nXDB",
    combo=="XFA-XFB" ~ "XFA-\nXFB",
    combo=="XBA-XCA" ~ "XBA-\nXCA",
    combo=="XCA-XDA" ~ "XCA-\nXDA",
    combo=="XDA-XFA" ~ "XDA-\nXFA",
    combo=="XFA-XBA" ~ "XFA-\nXBA"),
    combo=fct_relevel(combo, c("XBA-\nXBB","XCA-\nXCB","XDA-\nXDB","XFA-\nXFB","XBA-\nXCA","XCA-\nXDA","XDA-\nXFA","XFA-\nXBA"))) %>% 
  ggplot() +
  geom_bar(aes(x=ratio, y=avgAbundance, fill=categorySum), stat="identity", color="black", size=0.01) +
  scale_fill_manual(name="Category",
                    labels=c("noisy"="Noisy",
                             "flat"="Flat",
                             "actualFlat"="Actual flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("noisy"="#999999",
                             "flat"="#616667",
                             "actualFlat"="#120309",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#3C91E6",
                             "dose-dependent"="#A63A43")) +
  xlab("Community") +
  ylab("Relative abundance per category") +
  DEFAULTS.THEME_PRINT +
  facet_grid(~combo, scales="free") +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        axis.text.x=element_text(angle=90, vjust=.2, hjust=-0.4))
p_ASVcategoriesSumRelAb
