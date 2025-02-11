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

# Plot ASV counts ---------------------------------------------------------

# Plot distribution of parent and 1:1 points separately.
set.seed(1)
p_peakASVdistributionsP5 <- theoretical_mixtures_df %>% 
  ungroup() %>% 
  filter(passage==5 & (ratio %in% c("1:0","1:1","0:1") | mixtureType=="actual")) %>% 
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
  # Remove duplicates of parent communities.
  mutate(parent = case_when(ratio=="1:0" ~ sub("-.*","",combo),
                            ratio=="0:1" ~ sub(".*-","",combo),
                            .default=paste0(combo, "-", ratio, "-", mixtureType, "-", inoculationReplicate))) %>% 
  group_by(parent, inoculationReplicate) %>% 
  slice(1) %>% 
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
p_peakASVdistributionsP5
save_plot(paste0(plotdir, "/peakASVdistributions.pdf"), p_peakASVdistributionsP5, base_width=2.35, base_height=1.8)
#save_plot(paste0(plotdir, "/peakASVdistributions3ratios.png"), p_peakASVdistributions3ratios, 
#          base_height=1.5, base_width=2)

# Use Wilcoxon test to calculate significance of group comparisons for ASV counts.
parentTheoreticalASVCountTest <- wilcox.test(
  theoretical_mixtures_df %>% ungroup() %>% filter(passage==5 & ratio %in% c("1:0","0:1")) %>% 
    mutate(parent = ifelse(ratio=="1:0", sub("-.*","",combo), sub(".*-","",combo))) %>% 
    select(parent, inoculationReplicate, ASVcountRelAb3) %>% unique() %>% pull(ASVcountRelAb3),
  theoretical_mixtures_df %>% ungroup() %>% filter(passage==5 & ratio=="1:1" & mixtureType=="theoretical") %>% 
    select(combo, ratio, inoculationReplicate, ASVcountRelAb3) %>% unique() %>% pull(ASVcountRelAb3),
  alternative="less")
theoretical1to1CountTest <- wilcox.test(
  theoretical_mixtures_df %>% ungroup() %>% filter(passage==5 & ratio=="1:1" & mixtureType=="theoretical") %>% 
    select(combo, ratio, inoculationReplicate, ASVcountRelAb3) %>% unique() %>% pull(ASVcountRelAb3),
  theoretical_mixtures_df %>% ungroup() %>% filter(passage==5 & ratio=="1:1" & mixtureType=="actual") %>% 
    select(combo, ratio, inoculationReplicate, ASVcountRelAb3) %>% unique() %>% pull(ASVcountRelAb3),
  alternative="greater")
theoretical1to10CountTest <- wilcox.test(
  theoretical_mixtures_df %>% ungroup() %>% filter(passage==5 & ratio=="1:1" & mixtureType=="theoretical") %>% 
    select(combo, ratio, inoculationReplicate, ASVcountRelAb3) %>% unique() %>% pull(ASVcountRelAb3),
  theoretical_mixtures_df %>% ungroup() %>% filter(passage==5 & ratio %in% c("10:1","1:10") & mixtureType=="actual") %>% 
    select(combo, ratio, inoculationReplicate, ASVcountRelAb3) %>% unique() %>% pull(ASVcountRelAb3),
  alternative="greater")
theoretical1to100CountTest <- wilcox.test(
  theoretical_mixtures_df %>% ungroup() %>% filter(passage==5 & ratio=="1:1" & mixtureType=="theoretical") %>% 
    select(combo, ratio, inoculationReplicate, ASVcountRelAb3) %>% unique() %>% pull(ASVcountRelAb3),
  theoretical_mixtures_df %>% ungroup() %>% filter(passage==5 & ratio %in% c("100:1","1:100") & mixtureType=="actual") %>% 
    select(combo, ratio, inoculationReplicate, ASVcountRelAb3) %>% unique() %>% pull(ASVcountRelAb3),
  alternative="greater")
theoretical1to1000CountTest <- wilcox.test(
  theoretical_mixtures_df %>% ungroup() %>% filter(passage==5 & ratio=="1:1" & mixtureType=="theoretical") %>% 
    select(combo, ratio, inoculationReplicate, ASVcountRelAb3) %>% unique() %>% pull(ASVcountRelAb3),
  theoretical_mixtures_df %>% ungroup() %>% filter(passage==5 & ratio %in% c("1000:1","1:1000") & mixtureType=="actual") %>% 
    select(combo, ratio, inoculationReplicate, ASVcountRelAb3) %>% unique() %>% pull(ASVcountRelAb3),
  alternative="greater")

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

# Plot Shannon diversity. -------------------------------------------------

# Plot effective # species from Shannon, example plot with one theoretical and actual mixture.
p_ShannonExampleXFAXFBtheoretical <- theoretical_alpha %>%
  filter(combo=="XFA-XFB" & passage==5) %>%
  filter(mixtureType=="theoretical") %>% 
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
p_ShannonExampleXFAXFBtheoretical
save_plot(paste0(plotdir, "/ShannonExampleXFAXFB.pdf"), p_ShannonExampleXFAXFB, base_width=1.7, base_height=1.5)
#save_plot(paste0(plotdir, "/ShannonExampleXFAXFB.pdf"), p_ShannonExampleXFAXFB, base_width=2, base_height=1.75)
save_plot(paste0(plotdir, "/ShannonExampleXFAXFBtheoretical.png"), p_ShannonExampleXFAXFBtheoretical, base_width=2, base_height=1.75)

ASVlinesLegendHorizontal <- get_legend(p_ShannonExampleXFAXFB)
save_plot(paste0(plotdir, "/ASVlinesLegendHorizontal.pdf"), ASVlinesLegendHorizontal, base_width=3.5, base_height=0.5)

# Plot distribution of parent and 1:1 points separately for Shannon effective # species.
# Effective # of species is exp(Shannon)
p_peakShannonDistributions <- theoretical_alpha %>% 
  # Remove poorly sequenced XFA-XBA sample.
  filter(!(passage==5 & combo=="XFA-XBA" & line=="actual-1" & ratio=="100:1")) %>% 
  filter(passage==5 & (ratio %in% c("1:0","1:1","0:1") | mixtureType=="actual")) %>% 
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
  # Remove duplicates of parent communities.
  mutate(parent = case_when(ratio=="1:0" ~ sub("-.*","",combo),
                            ratio=="0:1" ~ sub(".*-","",combo),
                            .default=paste0(combo, "-", ratio, "-", mixtureType, "-", inoculationReplicate))) %>% 
  group_by(parent, inoculationReplicate) %>% 
  slice(1) %>% 
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

# Statistical significance using Wilcoxon test.
shannonParentTheoreticalDiff <- wilcox.test(
  exp(theoretical_alpha %>% filter(passage==5 & ratio %in% c("1:0","0:1")) %>% 
        mutate(parent = ifelse(ratio=="1:0", sub("-.*","",combo), sub(".*-","",combo))) %>% 
        select(parent, inoculationReplicate, Shannon) %>% unique() %>% pull(Shannon)),
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
