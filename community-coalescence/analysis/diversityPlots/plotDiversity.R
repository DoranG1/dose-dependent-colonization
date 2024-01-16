library(dada2)
library(phyloseq)
library(tidyverse)
library(ggfittext)
library(seqinr)
library(cowplot)
library(exactRankTests)
library(forcats)
library(rlist)
library(foreach)

setwd("/oak/stanford/groups/relman/users/dorang/210827-e0017-16S")
source("analysis/generateMixturesDataframe.R")
plotdir <- "analysis/diversityPlots"

# Calculate alpha diversity. ------------------------------------------------------

# Calculate diversity metrics for each sample.
theoretical_alpha <- estimate_richness(ps_theoretical, 
                                       measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
theoretical_alpha <- theoretical_alpha %>% 
  mutate(sample=rownames(theoretical_alpha), sample=sub("\\.","-",sample),
         sample=sub("\\.","-",sample), sample=sub("\\.",":",sample),
         sample=sub("\\.","-",sample))
# Join metadata to alpha diversity.
theoretical_alpha <- theoretical_alpha %>% 
  left_join(df_sam_theoretical, by="sample")

# # Calculate diversity metrics for each sample, filtered to >10^-4.
# theoretical_alpha_relAb4 <- estimate_richness(ps_theoretical_relAb4, 
#                                        measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
# theoretical_alpha_relAb4 <- theoretical_alpha_relAb4 %>% 
#   mutate(sample=rownames(theoretical_alpha_relAb4), sample=sub("\\.","-",sample),
#          sample=sub("\\.","-",sample), sample=sub("\\.",":",sample),
#          sample=sub("\\.","-",sample))
# # Join metadata to alpha diversity.
# theoretical_alpha_relAb4 <- theoretical_alpha_relAb4 %>% 
#   left_join(df_sam_theoretical, by="sample")
# 
# # Calculate diversity metrics for each sample, filtered to >10^-3.
# theoretical_alpha_relAb3 <- estimate_richness(ps_theoretical_relAb3, 
#                                               measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
# theoretical_alpha_relAb3 <- theoretical_alpha_relAb3 %>% 
#   mutate(sample=rownames(theoretical_alpha_relAb3), sample=sub("\\.","-",sample),
#          sample=sub("\\.","-",sample), sample=sub("\\.",":",sample),
#          sample=sub("\\.","-",sample))
# # Join metadata to alpha diversity.
# theoretical_alpha_relAb3 <- theoretical_alpha_relAb3 %>% 
#   left_join(df_sam_theoretical, by="sample")

# Plot ASV counts ---------------------------------------------------------

# Plot histogram of ASVs per samples, filtered to relative abundance >10^-4.
p_ASVcountHistRelAb4 <- theoretical_mixtures_df %>%
  ggplot() +
  geom_histogram(aes(x=ASVcountRelAb4)) +
  DEFAULTS.THEME_PRES +
  facet_wrap(~passage)
p_ASVcountHistRelAb4
save_plot(paste0(plotdir, "/ASVcountHistRelAb4.png"), p_ASVcountHistRelAb4)

# Plot histogram of ASVs per samples, filtered to relative abundance >10^-3.
p_ASVcountHistRelAb3 <- theoretical_mixtures_df %>%
  ggplot() +
  geom_histogram(aes(x=ASVcountRelAb3)) +
  DEFAULTS.THEME_PRES +
  facet_wrap(~passage)
p_ASVcountHistRelAb3
save_plot(paste0(plotdir, "/ASVcountHistRelAb3.png"), p_ASVcountHistRelAb3)

# Plot ASV counts for all mixtures, passage 3, rel ab >10^-4.
p_allASVCountsP3RelAb4 <- theoretical_mixtures_df %>%
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(
    x=ratio, y=ASVcountRelAb4, 
    group=line, 
    color=line)) +
  xlab("Community") +
  ylab("# ASVs") +
  ggtitle("ASV counts, passage 3") +
  ylim(c(0, 150)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_allASVCountsP3RelAb4
save_plot(paste0(plotdir, "/allASVCountsP3RelAb4.png"), p_allASVCountsP3RelAb4, nrow=1.5, ncol=1.75)

# Plot ASV counts for all mixtures, passage 5, rel ab >10^-4
p_allASVCountsP5RelAb4 <- theoretical_mixtures_df %>%
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(
    x=ratio, y=ASVcountRelAb4, 
    group=line, 
    color=line)) +
  xlab("Community") +
  ylab("# ASVs") +
  ggtitle("ASV counts, passage 5") +
  ylim(c(0, 150)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_allASVCountsP5RelAb4
save_plot(paste0(plotdir, "/allASVCountsP5RelAb4.png"), p_allASVCountsP5RelAb4, nrow=1.5, ncol=1.75)

# Plot ASV counts for all mixtures, passage 3, rel ab >10^-3.
p_allASVCountsP3RelAb3 <- theoretical_mixtures_df %>%
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(
    x=ratio, y=ASVcountRelAb3, 
    group=line, 
    color=line)) +
  xlab("Community") +
  ylab("# ASVs") +
  ggtitle("ASV counts, passage 3") +
  ylim(c(0, 150)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_allASVCountsP3RelAb3
save_plot(paste0(plotdir, "/allASVCountsP3RelAb3.png"), p_allASVCountsP3RelAb3, nrow=1.5, ncol=1.75)

# Plot ASV counts for all mixtures, passage 5, rel ab >10^-3.
p_allASVCountsP5RelAb3 <- theoretical_mixtures_df %>%
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(
    x=ratio, y=ASVcountRelAb3, 
    group=line, 
    color=line)) +
  xlab("Community") +
  ylab("# ASVs") +
  #ggtitle("ASV counts, passage 5") +
  ylim(c(0, 80)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_allASVCountsP5RelAb3
save_plot(paste0(plotdir, "/allASVCountsP5RelAb3.png"), p_allASVCountsP5RelAb3, nrow=1.5, ncol=1.75)

# Plot peak ASV counts, passage 5, rel ab >10^-4
p_peakASVCountsRelAb4 <- theoretical_mixtures_df %>%
  filter(ratio %in% c("1:0","1:1","0:1")) %>%
  rbind(theoretical_mixtures_df %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          mutate(mixtureType="actual")) %>% 
  group_by(combo, passage, ratio, mixtureType) %>% 
  summarize(ASVcountRelAb4_avg=mean(ASVcountRelAb4),
            ASVcountRelAb4_range=max(ASVcountRelAb4)-min(ASVcountRelAb4),
            ASVcountRelAb4_min=ASVcountRelAb4_avg-ASVcountRelAb4_range/2,
            ASVcountRelAb4_max=ASVcountRelAb4_avg+ASVcountRelAb4_range/2) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=ASVcountRelAb4_avg, group=interaction(passage, combo,mixtureType), color=mixtureType)) +
  geom_errorbar(aes(x=ratio, ymin=ASVcountRelAb4_min, ymax=ASVcountRelAb4_max, color=mixtureType), 
                alpha=0.25, width=0.05) +
  xlab("Community") +
  ylab("# ASVs") +
  ggtitle("ASV counts") +
  ylim(c(0, 150)) +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical"="Theoretical",
                              "actual"="Actual"), 
                     values=c("theoretical"="#313695", 
                              "actual"="#D73027")) +
  DEFAULTS.THEME_PRES + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_peakASVCountsRelAb4
save_plot(paste0(plotdir, "/peakASVCountsRelAb4.png"), p_peakASVCountsRelAb4, nrow=1.5, ncol=1.5)

# Theoretical vs actual ASV counts scatter plot.
ASVcountComparison_df <- theoretical_mixtures_df %>% 
  filter(ratio %in% c("1:0","1:1","0:1") & mixtureType=="theoretical") %>% 
  rename(ASVcountRelAb4_T=ASVcountRelAb4) %>% 
  group_by(passage, combo, ratio) %>% 
  summarize(ASVcountRelAb4_avgT=mean(ASVcountRelAb4_T),
         ASVcountRelAb4_rangeT=max(ASVcountRelAb4_T)-min(ASVcountRelAb4_T),
         ASVcountRelAb4_minT=ASVcountRelAb4_avgT-ASVcountRelAb4_rangeT/2,
         ASVcountRelAb4_maxT=ASVcountRelAb4_avgT+ASVcountRelAb4_rangeT/2) %>%
  left_join(
    theoretical_mixtures_df %>% 
      mutate(mixtureType=replace(mixtureType, ratio %in% c("1:0","0:1"), "actual")) %>% 
      filter(ratio %in% c("1:0","1:1","0:1") & mixtureType=="actual") %>% 
      rename(ASVcountRelAb4_A=ASVcountRelAb4) %>% 
      group_by(passage, combo, ratio) %>% 
      summarize(ASVcountRelAb4_avgA=mean(ASVcountRelAb4_A),
             ASVcountRelAb4_rangeA=max(ASVcountRelAb4_A)-min(ASVcountRelAb4_A),
             ASVcountRelAb4_minA=ASVcountRelAb4_avgA-ASVcountRelAb4_rangeA/2,
             ASVcountRelAb4_maxA=ASVcountRelAb4_avgA+ASVcountRelAb4_rangeA/2)
  )

p_ASVcountTheoreticalActualScatter <- ASVcountComparison_df %>% 
  ggplot() +
  geom_point(aes(x=ASVcountRelAb4_avgT, y=ASVcountRelAb4_avgA)) +
  geom_errorbar(aes(x=ASVcountRelAb4_avgT, ymin=ASVcountRelAb4_minA, ymax=ASVcountRelAb4_maxA)) +
  geom_errorbarh(aes(y=ASVcountRelAb4_avgA, xmin=ASVcountRelAb4_minT, xmax=ASVcountRelAb4_maxT), alpha=0.25) +
  geom_abline(slope=1, linetype="dotted") +
  xlab("Theoretical ASV count") +
  ylab("Actual ASV count") +
  xlim(c(0,150)) +
  ylim(c(0,150)) +
  ggtitle("Theoretical vs actual ASV counts") +
  DEFAULTS.THEME_PRES
p_ASVcountTheoreticalActualScatter
save_plot(paste0(plotdir, "/ASVcountTheoreticalActualScatter.png"), p_ASVcountTheoreticalActualScatter, nrow=1, ncol=1)

# Plot distribution of parent and 1:1 points separately.
set.seed(1)
p_peakASVdistributionsP3 <- theoretical_mixtures_df %>% 
  ungroup() %>% 
  filter(passage==3 & ratio %in% c("1:0","1:1","0:1") | mixtureType=="actual") %>% 
  #filter(ratio %in% c("1:0","1:1","0:1")) %>% 
  select(passage, combo, ratio, mixtureType, inoculationReplicate, ASVcountRelAb3) %>% 
  unique() %>%
  mutate(type=case_when(
    ratio %in% c("1:0","0:1") ~ "parent",
    ratio=="1:1" & mixtureType=="theoretical" ~ "1:1\n",
    ratio=="1:1" & mixtureType=="actual" ~ "1:1\n ",
    ratio %in% c("10:1","1:10") ~ "1:10",
    ratio %in% c("100:1","1:100") ~ "1:100",
    ratio %in% c("1000:1","1:1000") ~ "1:1000"),
    type=fct_relevel(
      type, c("parent","1:1\n","1:1\n ","1:10","1:100","1:1000"))) %>% 
  #mutate(type=case_when(ratio %in% c("1:0","0:1") ~ "parent",
  #                      ratio=="1:1" & mixtureType=="theoretical" ~ "theoretical\n1:1",
  #                      ratio=="1:1" & mixtureType=="actual" ~ "actual\n1:1"),
  #       type=fct_relevel(type, c("parent","theoretical\n1:1","actual\n1:1"))) %>% 
  ggplot() +
  #geom_jitter(aes(x=type, y=ASVcountRelAb3, color=type), width=0.15, size=0.25) +
  #geom_jitter(aes(x=type, y=ASVcountRelAb3, color=type), width=0.15, size=0.5) +
  geom_violin(aes(x=type, y=ASVcountRelAb3, fill=type), size=0.25) +
  #geom_boxplot(aes(x=type, y=ASVcountRelAb3), alpha=0, size=0.25) +
  geom_boxplot(aes(x=type, y=ASVcountRelAb3), size=0.25, width=0.1, outlier.shape=NA) +
  scale_fill_manual(values=c("parent"="#4DAF4A",
                              "1:1\n"="#808080",
                              "1:1\n "="#F46D43",
                              "1:10"="#F46D43",
                              "1:100"="#F46D43",
                              "1:1000"="#F46D43")) +
  #scale_color_manual(values=c("parent"="#4DAF4A",
  #                            "theoretical\n1:1"="#808080",
  #                            "actual\n1:1"="#F46D43")) +
  theme(legend.position="none") +
  #ylim(c(0,150)) +
  ylim(c(0,80)) +
  xlab("Community") +
  ylab("# ASVs") +
  #ggtitle("ASV counts") +
  DEFAULTS.THEME_PRINT
  #DEFAULTS.THEME_PRES
p_peakASVdistributionsP3
save_plot(paste0(plotdir, "/peakASVdistributionsP3.pdf"), p_peakASVdistributionsP3, base_width=2.35, base_height=1.8)
#save_plot(paste0(plotdir, "/peakASVdistributions3ratios.png"), p_peakASVdistributions3ratios, 
#          base_height=1.5, base_width=2)

# Statistical significance of differences across groups.
aovdata <- theoretical_mixtures_df %>% 
  ungroup() %>%
  filter(passage==5 & (ratio %in% c("1:0","1:1","0:1") | mixtureType=="actual")) %>% 
  select(passage, combo, ratio, mixtureType, inoculationReplicate, ASVcountRelAb3) %>% 
  unique() %>%
  mutate(type=case_when(
    ratio %in% c("1:0","0:1") ~ "parent",
    ratio=="1:1" & mixtureType=="theoretical" ~ "1:1 theoretical",
    ratio=="1:1" & mixtureType=="actual" ~ "1:1 actual ",
    ratio %in% c("10:1","1:10") ~ "1:10",
    ratio %in% c("100:1","1:100") ~ "1:100",
    ratio %in% c("1000:1","1:1000") ~ "1:1000")) %>% 
  filter(!is.na(type) & !is.na(ASVcountRelAb3))
summaryAOV <- aov(ASVcountRelAb3 ~ type, data=aovdata)
summary(summaryAOV)
summaryTukey <- TukeyHSD(summaryAOV, conf.level=0.95)

# Plot comparison of ASV counts for passages.
p_ASVcountComparison <- theoretical_mixtures_df_bottomed %>%
  filter(mixtureType=="actual") %>%
  ggplot() +
  geom_boxplot(aes(x=factor(passage), y=ASVcount)) +
  geom_point(aes(x=factor(passage), y=ASVcount, fill=factor(passage)), pch=21, alpha=0.1) +
  scale_fill_manual(values=c("blue", "red")) +
  theme(legend.position="none") +
  ylab("ASV counts") +
  xlab("Passage") +
  facet_grid(~combo) +
  DEFAULTS.THEME_PRES
p_ASVcountComparison

# Plot ASV counts for all mixtures, both passages.
p_allASVCounts <- theoretical_mixtures_df_bottomed %>%
  ggplot() +
  geom_line(aes(
    x=ratio, y=ASVcount, 
    group=line, 
    color=line)) +
  xlab("Community") +
  ylab("# ASVs") +
  ggtitle("ASV counts, passage 5") +
  ylim(c(0, 150)) +
  scale_color_manual(name="Mixture type",
                     labels=c("Theoretical 1", "Theoretical 2", "Actual 1", "Actual 2", "Actual 3"), 
                     values=c("chartreuse2", "green1", "slateblue4", "steelblue4", "royalblue2")) +
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo~passage, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/allASVCounts.png"), p_allASVCounts, nrow=2.25, ncol=1.75)

# Plot ASV counts, example plot with one theoretical mixture.
p_ASVCountExampleTheoretical <- theoretical_mixtures_df %>%
  filter(combo=="XFA-XFB" & passage==5 & mixtureType=="theoretical") %>%
  mutate(lineGroup=case_when(ratio=="1:0" ~ "leftParent",
                             ratio=="0:1" ~ "rightParent",
                             ratio %in% ratios ~ "mixture")) %>% 
  ggplot() +
  geom_line(aes(
    x=ratio, y=ASVcountRelAb3, 
    group=interaction(lineGroup, line), 
    color=line)) +
  geom_point(aes(x=ratio, y=ASVcountRelAb3, color=group), size=1) +
  xlab("Community") +
  ylab("# ASVs") +
  #ggtitle("ASV counts, passage 5") +
  ylim(c(0, 80)) +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "mixture"="Mixture",
                              "control"="Parent"), 
                     values=c("theoretical-1"="#313695", 
                              "theoretical-2"="#4575B4",
                              "mixture"="#313695",
                              "control"="#4DAF4A")) +
  #DEFAULTS.THEME_PRES + 
  DEFAULTS.THEME_PRINT +
  #facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  theme(legend.position="none")
p_ASVCountExampleTheoretical
save_plot(paste0(plotdir, "/ASVCountExampleTheoreticalXFAXFB.png"), p_ASVCountExampleTheoretical, 
          base_width=2.35, base_height=1.8)

# Plot ASV counts, example plot with one theoretical and actual mixture.
p_ASVCountExampleXBAXCA <- theoretical_mixtures_df %>%
  filter(combo=="XBA-XCA" & passage==5) %>%
  filter(!(combo=="XFA-XBA" & ratio=="100:1" & inoculationReplicate==1)) %>% 
  mutate(combo="A1/B1") %>% 
  group_by(label, line) %>% 
  summarize(combo, group, ratio, inoculationReplicate, ASVcountRelAb3) %>% 
  unique() %>% 
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
  mutate(lineGroup = fct_relevel(lineGroup, c("parent","theoretical-1","theoretical-2","actual-1","actual-2","actual-3"))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=ASVcountRelAb3, color=lineGroup, group=interaction(lineGroup, parentGroup)),
            alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=ASVcountRelAb3, color=lineGroup), alpha=0.8, size=0.5) +
  #geom_line(aes(
  #  x=ratio, y=ASVcountRelAb3, 
  #  group=line, 
  #  color=line)) +
  xlab("Community") +
  ylab("# ASVs") +
  ylim(c(0, 80)) +
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
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  DEFAULTS.THEME_PRINT + 
  #DEFAULTS.THEME_PRES +
  facet_wrap(~combo) +
  #theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=5.5),
        legend.title=element_blank(),
        #legend.position=c(0.68,1),
        legend.key.size=unit(0.5, "lines")) +
  theme(legend.position="none")
p_ASVCountExampleXBAXCA
save_plot(paste0(plotdir, "/ASVCountExampleXBAXCA.pdf"), p_ASVCountExampleXBAXCA, base_width=1.75, base_height=1.5)
# save_plot(paste0(plotdir, "/ASVCountExample.png"), p_ASVCountExample, nrow=1, ncol=1)

ASVlinesLegend <- get_legend(p_ASVCountExampleXFAXFB)
save_plot(paste0(plotdir, "/ASVlinesLegend.pdf"), ASVlinesLegend, base_width=0.75, base_height=.8)

save_plot(paste0(plotdir, "/ASVcountsPaperCombined.png"), p_ASVCountExample / p_peakASVdistributions, 
                 nrow=1.25, ncol=1)

# Plot observed ASVs without filtering. -----------------------------------

# Plot observed ASVs for all mixtures, passage 3.
p_allObservedP3 <- theoretical_alpha %>% 
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Observed, group=line, color=line)) +
  xlab("Community") +
  ylab("Observed ASVs") +
  ggtitle("Observed ASVs, passage 3") +
  ylim(c(0, 200)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/observedASVsP3.png"), p_allObservedP3, nrow=1.5, ncol=1.75)

# Plot observed ASVs for all mixtures, passage 5.
p_allObservedP5 <- theoretical_alpha %>% 
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Observed, group=line, color=line)) +
  xlab("Community") +
  ylab("Observed ASVs") +
  ggtitle("Observed ASVs, passage 5") +
  ylim(c(0, 200)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/observedASVsP5.png"), p_allObservedP5, nrow=1.5, ncol=1.75)

# Plot Shannon diversity. -------------------------------------------------

# Plot effective # species from Shannon, example plot with one theoretical and actual mixture.
p_ShannonExampleXFAXFB <- theoretical_alpha %>%
  filter(combo=="XFA-XFB" & passage==5) %>%
  filter(!(combo=="XFA-XBA" & ratio=="100:1" & inoculationReplicate==1)) %>% 
  mutate(combo="D1/D2") %>% 
  group_by(label, line) %>% 
  summarize(combo, group, ratio, inoculationReplicate, Shannon) %>% 
  unique() %>% 
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
  mutate(lineGroup=fct_relevel(lineGroup, c("parent","theoretical-1","theoretical-2","actual-1","actual-2","actual-3"))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=exp(Shannon), color=lineGroup, group=interaction(lineGroup, parentGroup)),
            alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=exp(Shannon), color=lineGroup), alpha=0.8, size=0.5) +
  #geom_line(aes(
  #  x=ratio, y=ASVcountRelAb3, 
  #  group=line, 
  #  color=line)) +
  xlab("Mixture ratio") +
  ylab("Effective # species") +
  ylim(c(0, 40)) +
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
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  DEFAULTS.THEME_PRINT + 
  #DEFAULTS.THEME_PRES +
  facet_wrap(~combo) +
  #theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=5.5),
        legend.title=element_blank(),
        #legend.position=c(0.68,1),
        legend.key.size=unit(0.5, "lines"),
        strip.text=element_text(face="bold")) +
  #theme(legend.position="bottom") +
  guides(color = guide_legend(nrow=1)) +
  theme(legend.position="none")
p_ShannonExampleXFAXFB
save_plot(paste0(plotdir, "/ShannonExampleXFAXFB.pdf"), p_ShannonExampleXFAXFB, base_width=1.7, base_height=1.5)
#save_plot(paste0(plotdir, "/ShannonExampleXFAXFB.pdf"), p_ShannonExampleXFAXFB, base_width=2, base_height=1.75)

ASVlinesLegendHorizontal <- get_legend(p_ShannonExampleXFAXFB)
save_plot(paste0(plotdir, "/ASVlinesLegendHorizontal.pdf"), ASVlinesLegendHorizontal, base_width=3.5, base_height=0.5)

# Plot Shannon diversity for all mixtures, passage 3.
p_allShannonP3 <- theoretical_alpha %>% 
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Shannon, group=line, color=line)) +
  xlab("Community") +
  ylab("Shannon diversity") +
  ggtitle("Shannon diversity, passage 3") +
  ylim(c(2, 3.75)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allShannonP3
save_plot(paste0(plotdir, "/ShannonP3.png"), p_allShannonP3, nrow=1.5, ncol=1.75)

# Plot Shannon diversity for all mixtures, passage 5.
p_allShannonP5 <- theoretical_alpha %>% 
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Shannon, group=line, color=line)) +
  xlab("Community") +
  ylab("Shannon diversity") +
  ggtitle("Shannon diversity, passage 5") +
  ylim(c(2, 3.75)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allShannonP5
save_plot(paste0(plotdir, "/ShannonP5.png"), p_allShannonP5, nrow=1.5, ncol=1.75)

# Plot effective number of species from Shannon diversity, passage 3.
# Effective # of species is exp(Shannon).
p_allShannonP3Effective <- theoretical_alpha %>% 
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=exp(Shannon), group=line, color=line)) +
  xlab("Community") +
  ylab("Effective # species, Shannon") +
  ggtitle("Shannon effective species, passage 3") +
  ylim(c(0, 35)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allShannonP3Effective
save_plot(paste0(plotdir, "/ShannonP3Effective.png"), p_allShannonP3Effective, nrow=1.5, ncol=1.75)

# Plot effective number of species from Shannon diversity, passage 5.
# Effective # of species is exp(Shannon).
p_allShannonP5Effective <- theoretical_alpha %>% 
  filter(passage==5) %>%
  filter(!(combo=="XFA-XBA" & ratio=="100:1" & inoculationReplicate==1)) %>% 
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
  geom_line(aes(x=ratio, y=exp(Shannon), color=lineGroup, group=interaction(lineGroup, parentGroup)),
            alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=exp(Shannon), color=lineGroup), alpha=0.8, size=0.75) +
  xlab("Mixture ratio") +
  ylab("Effective # species") +
  ylim(c(0, 40)) +
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
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  DEFAULTS.THEME_PRINT + 
  #DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2) +
  #theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=5.5),
        legend.title=element_blank(),
        #legend.position=c(0.68,1),
        legend.key.size=unit(0.5, "lines")) +
  theme(legend.position="none") 
p_allShannonP5Effective
save_plot(paste0(plotdir, "/ShannonP5Effective.png"), p_allShannonP5Effective, nrow=1.5, ncol=1.75)

# Plot example effective number of species from Shannon diversity, passage 5.
# Effective # of species is exp(Shannon).
p_ShannonEffectiveExample <- theoretical_alpha %>% 
  filter(combo=="XDA-XDB" & passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=exp(Shannon), group=line, color=line), size=0.25) +
  xlab("Community") +
  ylab("Shannon effective # species") +
  ylim(c(0, 35)) +
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
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  DEFAULTS.THEME_PRINT + 
  facet_wrap(~combo) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=5.5),
        legend.title=element_blank(),
        legend.position=c(0.68,1),
        legend.key.size=unit(0.5, "lines")) 
p_ShannonEffectiveExample
save_plot(paste0(plotdir, "/ShannonEffectiveExample.pdf"), 
          p_ShannonEffectiveExample, base_width=2.25, base_height=1.8)

# Plot Shannon diversity for all mixtures, passage 3, filtered to >10^-4.
p_allShannonP3_relAb4 <- theoretical_alpha_relAb4 %>% 
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Shannon, group=line, color=line)) +
  xlab("Community") +
  ylab("Shannon diversity") +
  ggtitle("Shannon diversity, passage 3") +
  ylim(c(2, 3.75)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allShannonP3_relAb4
save_plot(paste0(plotdir, "/ShannonP3_relAb4.png"), p_allShannonP3_relAb4, nrow=1.5, ncol=1.75)

# Plot Shannon diversity for all mixtures, passage 5, filtered to >10^-4.
p_allShannonP5_relAb4 <- theoretical_alpha_relAb4 %>% 
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Shannon, group=line, color=line)) +
  xlab("Community") +
  ylab("Shannon diversity") +
  ggtitle("Shannon diversity, passage 5") +
  ylim(c(2, 3.75)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allShannonP5_relAb4
save_plot(paste0(plotdir, "/ShannonP5_relAb4.png"), p_allShannonP5_relAb4, nrow=1.5, ncol=1.75)

# Plot effective number of species from Shannon diversity, passage 3, filtered to >10^-4.
# Effective # of species is exp(Shannon).
p_allShannonP3Effective_relAb4 <- theoretical_alpha_relAb4 %>% 
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=exp(Shannon), group=line, color=line)) +
  xlab("Community") +
  ylab("Effective # species, Shannon") +
  ggtitle("Shannon effective species, passage 3") +
  ylim(c(0, 35)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allShannonP3Effective_relAb4
save_plot(paste0(plotdir, "/ShannonP3Effective_relAb4.png"), p_allShannonP3Effective_relAb4, nrow=1.5, ncol=1.75)

# Plot effective number of species from Shannon diversity, passage 5, filtered to >10^-4.
# Effective # of species is exp(Shannon).
p_allShannonP5Effective_relAb4 <- theoretical_alpha_relAb4 %>% 
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=exp(Shannon), group=line, color=line)) +
  xlab("Community") +
  ylab("Effective # species, Shannon") +
  ggtitle("Shannon effective species, passage 5") +
  ylim(c(0, 35)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allShannonP5Effective_relAb4
save_plot(paste0(plotdir, "/ShannonP5Effective_relAb4.png"), p_allShannonP5Effective_relAb4, nrow=1.5, ncol=1.75)

# Plot Shannon diversity for all mixtures, passage 3, filtered to >10^-3.
p_allShannonP3_relAb3 <- theoretical_alpha_relAb3 %>% 
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Shannon, group=line, color=line)) +
  xlab("Community") +
  ylab("Shannon diversity") +
  ggtitle("Shannon diversity, passage 3") +
  ylim(c(2, 3.75)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allShannonP3_relAb3
save_plot(paste0(plotdir, "/ShannonP3_relAb3.png"), p_allShannonP3_relAb3, nrow=1.5, ncol=1.75)

# Plot Shannon diversity for all mixtures, passage 5, filtered to >10^-3.
p_allShannonP5_relAb3 <- theoretical_alpha_relAb3 %>% 
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Shannon, group=line, color=line)) +
  xlab("Community") +
  ylab("Shannon diversity") +
  ggtitle("Shannon diversity, passage 5") +
  ylim(c(2, 3.75)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allShannonP5_relAb3
save_plot(paste0(plotdir, "/ShannonP5_relAb3.png"), p_allShannonP5_relAb3, nrow=1.5, ncol=1.75)

# Plot effective number of species from Shannon diversity, passage 3, filtered to >10^-3.
# Effective # of species is exp(Shannon).
p_allShannonP3Effective_relAb3 <- theoretical_alpha_relAb3 %>% 
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=exp(Shannon), group=line, color=line)) +
  xlab("Community") +
  ylab("Effective # species, Shannon") +
  ggtitle("Shannon effective species, passage 3") +
  ylim(c(0, 35)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allShannonP3Effective_relAb3
save_plot(paste0(plotdir, "/ShannonP3Effective_relAb3.png"), p_allShannonP3Effective_relAb3, nrow=1.5, ncol=1.75)

# Plot effective number of species from Shannon diversity, passage 5, filtered to >10^-3.
# Effective # of species is exp(Shannon).
p_allShannonP5Effective_relAb3 <- theoretical_alpha_relAb3 %>% 
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=exp(Shannon), group=line, color=line)) +
  xlab("Community") +
  ylab("Effective # species, Shannon") +
  ggtitle("Shannon effective species, passage 5") +
  ylim(c(0, 35)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allShannonP5Effective_relAb3
save_plot(paste0(plotdir, "/ShannonP5Effective_relAb3.png"), p_allShannonP5Effective_relAb3, nrow=1.5, ncol=1.75)

# Plot distribution of parent and 1:1 points separately for Shannon effective # species.
# Effective # of species is exp(Shannon)
p_peakShannonDistributions <- theoretical_alpha %>% 
  # Remove poorly sequenced XFA-XBA sample.
  filter(!(passage==5 & combo=="XFA-XBA" & line=="actual-1" & ratio=="100:1")) %>% 
  filter(passage==5 & ratio %in% c("1:0","1:1","0:1") | mixtureType=="actual") %>% 
  select(passage, combo, ratio, mixtureType, inoculationReplicate, Shannon) %>% 
  unique() %>%
  mutate(type=case_when(
    ratio %in% c("1:0","0:1") ~ "parent",
    ratio=="1:1" & mixtureType=="theoretical" ~ "1:1\n",
    ratio=="1:1" & mixtureType=="actual" ~ "1:1\n ",
    ratio %in% c("10:1","1:10") ~ "1:10",
    ratio %in% c("100:1","1:100") ~ "1:100",
    ratio %in% c("1000:1","1:1000") ~ "1:1000"),
    type=fct_relevel(
      type, c("parent","1:1\n","1:1\n ","1:10","1:100","1:1000"))) %>% 
  ggplot() +
  #geom_jitter(aes(x=type, y=exp(Shannon), color=type), width=0.15, size=0.25) +
  geom_violin(aes(x=type, y=exp(Shannon), fill=type), size=0.25) +
  #geom_boxplot(aes(x=type, y=exp(Shannon)), alpha=0, size=0.25) +
  geom_boxplot(aes(x=type, y=exp(Shannon)), size=0.25, width=0.1, outlier.shape=NA) +
  scale_fill_manual(values=c("parent"="#4DAF4A",
                             "1:1\n"="#808080",
                             "1:1\n "="#F46D43",
                             "1:10"="#F46D43",
                             "1:100"="#F46D43",
                             "1:1000"="#F46D43")) +
  theme(legend.position="none") +
  ylim(c(0,40)) +
  xlab("Mixture ratio") +
  ylab("Effective # species") +
  DEFAULTS.THEME_PRINT
p_peakShannonDistributions
save_plot(paste0(plotdir, "/peakShannonDistributions.pdf"), p_peakShannonDistributions, base_width=2.35, base_height=1.8)

# Statistical significance of differences across groups.
aovdataShannon <- theoretical_alpha %>% 
  filter(!(passage==5 & combo=="XFA-XBA" & line=="actual-1" & ratio=="100:1")) %>% 
  ungroup() %>%
  filter(passage==5 & (ratio %in% c("1:0","1:1","0:1") | mixtureType=="actual")) %>% 
  select(passage, combo, ratio, mixtureType, inoculationReplicate, Shannon) %>% 
  unique() %>%
  mutate(type=case_when(
    ratio %in% c("1:0","0:1") ~ "parent",
    ratio=="1:1" & mixtureType=="theoretical" ~ "1:1 theoretical",
    ratio=="1:1" & mixtureType=="actual" ~ "1:1 actual ",
    ratio %in% c("10:1","1:10") ~ "1:10",
    ratio %in% c("100:1","1:100") ~ "1:100",
    ratio %in% c("1000:1","1:1000") ~ "1:1000"))
summaryAOVShannon <- aov(exp(Shannon) ~ type, data=aovdataShannon)
summary(summaryAOVShannon)
summaryTukeyShannon <- TukeyHSD(summaryAOVShannon, conf.level=0.95)

# Statistical significance using Wilcoxon test.
shannonParentTheoreticalDiff <- wilcox.test(
  exp(theoretical_alpha %>% filter(passage==5 & ratio %in% c("1:0","0:1")) %>% pull(Shannon)),
  exp(theoretical_alpha %>% filter(passage==5 & ratio=="1:1" & mixtureType=="theoretical") %>% pull(Shannon)),
  alternative = c("less"))

shannon11TheoreticalDiff <- wilcox.test(
  exp(theoretical_alpha %>% filter(passage==5 & ratio=="1:1" & mixtureType=="actual") %>% pull(Shannon)),
  exp(theoretical_alpha %>% filter(passage==5 & ratio=="1:1" & mixtureType=="theoretical") %>% pull(Shannon)),
  alternative = c("less"))

shannon110TheoreticalDiff <- wilcox.test(
  exp(theoretical_alpha %>% filter(passage==5 & ratio %in% c("1:10","10:1") & mixtureType=="actual") %>% pull(Shannon)),
  exp(theoretical_alpha %>% filter(passage==5 & ratio=="1:1" & mixtureType=="theoretical") %>% pull(Shannon)),
  alternative = c("less"))

shannon1100TheoreticalDiff <- wilcox.test(
  exp(theoretical_alpha %>% filter(passage==5 & ratio %in% c("1:100","100:1") & mixtureType=="actual") %>% pull(Shannon)),
  exp(theoretical_alpha %>% filter(passage==5 & ratio=="1:1" & mixtureType=="theoretical") %>% pull(Shannon)),
  alternative = c("less"))

shannon11000TheoreticalDiff <- wilcox.test(
  exp(theoretical_alpha %>% filter(passage==5 & ratio %in% c("1:1000","1000:1") & mixtureType=="actual") %>% pull(Shannon)),
  exp(theoretical_alpha %>% filter(passage==5 & ratio=="1:1" & mixtureType=="theoretical") %>% pull(Shannon)),
  alternative = c("less"))

# Check normality of values.
theoretical_alpha %>% filter(passage==5 & ratio=="1:1" & mixtureType=="theoretical") %>% 
  ggplot() +
  geom_histogram(aes(x=exp(Shannon)))  +
  xlim(0,45)

p_peakASVdistributionsP3 <- theoretical_mixtures_df %>% 
  ungroup() %>% 
  filter(passage==3 & ratio %in% c("1:0","1:1","0:1") | mixtureType=="actual") %>% 
  #filter(ratio %in% c("1:0","1:1","0:1")) %>% 
  select(passage, combo, ratio, mixtureType, inoculationReplicate, ASVcountRelAb3) %>% 
  unique() %>%
  mutate(type=case_when(
    ratio %in% c("1:0","0:1") ~ "parent",
    ratio=="1:1" & mixtureType=="theoretical" ~ "1:1\n",
    ratio=="1:1" & mixtureType=="actual" ~ "1:1\n ",
    ratio %in% c("10:1","1:10") ~ "1:10",
    ratio %in% c("100:1","1:100") ~ "1:100",
    ratio %in% c("1000:1","1:1000") ~ "1:1000"),
    type=fct_relevel(
      type, c("parent","1:1\n","1:1\n ","1:10","1:100","1:1000"))) %>% 
  #mutate(type=case_when(ratio %in% c("1:0","0:1") ~ "parent",
  #                      ratio=="1:1" & mixtureType=="theoretical" ~ "theoretical\n1:1",
  #                      ratio=="1:1" & mixtureType=="actual" ~ "actual\n1:1"),
  #       type=fct_relevel(type, c("parent","theoretical\n1:1","actual\n1:1"))) %>% 
  ggplot() +
  #geom_jitter(aes(x=type, y=ASVcountRelAb3, color=type), width=0.15, size=0.25) +
  #geom_jitter(aes(x=type, y=ASVcountRelAb3, color=type), width=0.15, size=0.5) +
  geom_violin(aes(x=type, y=ASVcountRelAb3, fill=type), size=0.25) +
  #geom_boxplot(aes(x=type, y=ASVcountRelAb3), alpha=0, size=0.25) +
  geom_boxplot(aes(x=type, y=ASVcountRelAb3), size=0.25, width=0.1, outlier.shape=NA) +
  scale_fill_manual(values=c("parent"="#4DAF4A",
                             "1:1\n"="#808080",
                             "1:1\n "="#F46D43",
                             "1:10"="#F46D43",
                             "1:100"="#F46D43",
                             "1:1000"="#F46D43")) +
  #scale_color_manual(values=c("parent"="#4DAF4A",
  #                            "theoretical\n1:1"="#808080",
  #                            "actual\n1:1"="#F46D43")) +
  theme(legend.position="none") +
  #ylim(c(0,150)) +
  ylim(c(0,80)) +
  xlab("Community") +
  ylab("# ASVs") +
  #ggtitle("ASV counts") +
  DEFAULTS.THEME_PRINT
#DEFAULTS.THEME_PRES
p_peakASVdistributionsP3
save_plot(paste0(plotdir, "/peakASVdistributionsP3.pdf"), p_peakASVdistributionsP3, base_width=2.35, base_height=1.8)

# Plot Simpson diversity. -------------------------------------------------

# Plot Simpson diversity for all mixtures, passage 3.
p_allSimpsonP3 <- theoretical_alpha %>% 
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Simpson, group=line, color=line)) +
  xlab("Community") +
  ylab("Simpson diversity") +
  ggtitle("Simpson diversity, passage 3") +
  ylim(c(0.80, 0.96)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allSimpsonP3
save_plot(paste0(plotdir, "/SimpsonP3.png"), p_allSimpsonP3, nrow=1.5, ncol=1.75)

# Plot Simpson diversity for all mixtures, passage 5.
p_allSimpsonP5 <- theoretical_alpha %>% 
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Simpson, group=line, color=line)) +
  xlab("Community") +
  ylab("Simpson diversity") +
  ggtitle("Simpson diversity, passage 5") +
  ylim(c(0.80, 0.96)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allSimpsonP5
save_plot(paste0(plotdir, "/SimpsonP5.png"), p_allSimpsonP5, nrow=1.5, ncol=1.75)

# Plot effective number of species from Simpson diversity, passage 3.
# Effective # of species is InvSimpson (1/(1-Simpson)).
p_allSimpsonP3Effective <- theoretical_alpha %>% 
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=InvSimpson, group=line, color=line)) +
  xlab("Community") +
  ylab("Effective # species, Simpson") +
  ggtitle("Simpson effective species, passage 3") +
  ylim(c(0, 23)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allSimpsonP3Effective
save_plot(paste0(plotdir, "/SimpsonP3Effective.png"), p_allSimpsonP3Effective, nrow=1.5, ncol=1.75)

# Plot effective number of species from Simpson diversity, passage 5.
# Effective # of species is InvSimpson (1/(1-Simpson)).
p_allSimpsonP5Effective <- theoretical_alpha %>% 
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=InvSimpson, group=line, color=line)) +
  xlab("Community") +
  ylab("Effective # species, Simpson") +
  ggtitle("Simpson effective species, passage 5") +
  ylim(c(0, 23)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allSimpsonP5Effective
save_plot(paste0(plotdir, "/SimpsonP5Effective.png"), p_allSimpsonP5Effective, nrow=1.5, ncol=1.75)

# Plot Simpson diversity for all mixtures, passage 3, filtered to >10^-4.
p_allSimpsonP3_relAb4 <- theoretical_alpha_relAb4 %>% 
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Simpson, group=line, color=line)) +
  xlab("Community") +
  ylab("Simpson diversity") +
  ggtitle("Simpson diversity, passage 3") +
  ylim(c(0.80, 0.96)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allSimpsonP3_relAb4
save_plot(paste0(plotdir, "/SimpsonP3_relAb4.png"), p_allSimpsonP3_relAb4, nrow=1.5, ncol=1.75)

# Plot Simpson diversity for all mixtures, passage 5, filtered to >10^-4.
p_allSimpsonP5_relAb4 <- theoretical_alpha_relAb4 %>% 
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=Simpson, group=line, color=line)) +
  xlab("Community") +
  ylab("Simpson diversity") +
  ggtitle("Simpson diversity, passage 5") +
  ylim(c(0.80, 0.96)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allSimpsonP5_relAb4
save_plot(paste0(plotdir, "/SimpsonP5_relAb4.png"), p_allSimpsonP5_relAb4, nrow=1.5, ncol=1.75)

# Plot effective number of species from Simpson diversity, passage 3, filtered to >10^-4.
# Effective # of species is InvSimpson (1/(1-Simpson)).
p_allSimpsonP3Effective_relAb4 <- theoretical_alpha_relAb4 %>% 
  filter(passage==3) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=InvSimpson, group=line, color=line)) +
  xlab("Community") +
  ylab("Effective # species, Simpson") +
  ggtitle("Simpson effective species, passage 3") +
  ylim(c(0, 23)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allSimpsonP3Effective_relAb4
save_plot(paste0(plotdir, "/SimpsonP3Effective_relAb4.png"), p_allSimpsonP3Effective_relAb4, nrow=1.5, ncol=1.75)

# Plot effective number of species from Simpson diversity, passage 5, filtered to >10^-4.
# Effective # of species is InvSimpson (1/(1-Simpson)).
p_allSimpsonP5Effective_relAb4 <- theoretical_alpha_relAb4 %>% 
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=InvSimpson, group=line, color=line)) +
  xlab("Community") +
  ylab("Effective # species, Simpson") +
  ggtitle("Simpson effective species, passage 5") +
  ylim(c(0, 23)) +
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
  DEFAULTS.THEME_PRES + 
  facet_wrap(~combo, scales="free", nrow=2, ncol=4) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))  
p_allSimpsonP5Effective_relAb4
save_plot(paste0(plotdir, "/SimpsonP5Effective_relAb4.png"), p_allSimpsonP5Effective_relAb4, nrow=1.5, ncol=1.75)