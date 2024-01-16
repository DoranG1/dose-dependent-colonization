library(dada2)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggfittext)
library(seqinr)
library(cowplot)
library(exactRankTests)
library(forcats)
library(rlist)
library(foreach)

setwd("/oak/stanford/groups/relman/users/dorang/210827-e0017-16S")
source("analysis/generateMixturesDataframe.R")
plotdir <- "analysis/paintedASVPlots"

# Paint ASVs by community of origin ---------------------------------------
# Example plot with one theoretical mixture.

p_abundancesPaintedTheoreticalExample <- theoretical_mixtures_df_bottomed %>%
  filter(combo=="XFA-XBA" & mixtureType=="theoretical" & inoculationReplicate==2) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parent)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/abundancesPaintedTheoreticalExample.png"), 
          p_abundancesPaintedTheoreticalExample, nrow=1, ncol=1.5) 

# Example plot with one theoretical and actual mixture.
p_abundancesPaintedExample <- theoretical_mixtures_df_bottomed %>%
  filter(combo=="XFA-XBA" & inoculationReplicate==2) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parent)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/abundancesPaintedExample.png"), p_abundancesPaintedExample, nrow=1, ncol=1.5)

# Example plot with one theoretical and actual mixture, black outlines.
p_abundancesPaintedExampleOutlines <- theoretical_mixtures_df_bottomed %>%
  filter(combo=="XFA-XBA" & inoculationReplicate==2) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parent)), stat="identity", color="black", size=0.01) +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/abundancesPaintedExampleOutlines.png"),
          p_abundancesPaintedExampleOutlines, nrow=1, ncol=1.5)

# Example plot 2 with one theoretical and actual mixture.
p_abundancesPaintedExample2 <- theoretical_mixtures_df_bottomed %>%
  filter(combo=="XCA-XDA" & inoculationReplicate==2) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parent)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/abundancesPaintedExample2.png"), p_abundancesPaintedExample2, nrow=1, ncol=1.5)

# Example plot with just actual mixtures from XDA-XDB.
p_abundancesPaintedExample3 <- theoretical_mixtures_df_bottomed %>%
  filter(combo=="XDA-XDB" & mixtureType=="actual") %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parent)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/abundancesPaintedExample3.png"), p_abundancesPaintedExample3, nrow=1, ncol=1.25)

# All theoretical mixtures painted.
p_abundancesPaintedTheoretical <- theoretical_mixtures_df_bottomed %>%
  filter(mixtureType=="theoretical") %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parent)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/abundancesPaintedTheoretical.png"), p_abundancesPaintedTheoretical, nrow=1.5, ncol=2)

# All painted ASV abundances, passage 3.
p_abundancesPaintedAllP3 <- theoretical_mixtures_df_bottomed %>%
  filter(passage==3) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parent)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/abundancesPaintedAllP3.png"), p_abundancesPaintedAllP3, nrow=1.5, ncol=2)  

# All painted ASV abundances, passage 5.
p_abundancesPaintedAllP5 <- theoretical_mixtures_df_bottomed %>%
  filter(passage==5) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parent)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/abundancesPaintedAllP5.png"), p_abundancesPaintedAllP5, nrow=1.5, ncol=2)  

# Paint all abundances with black outlines, passage 3.
p_abundancesPaintedBlackP3 <- theoretical_mixtures_df_bottomed %>%
  filter(passage==3) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parent)), stat="identity", color="black", size=0.001) +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/abundancesPaintedBlackP3.png"), p_abundancesPaintedBlackP3, nrow=1.5, ncol=2)

# Paint all abundances with black outlines, passage 3.
p_abundancesPaintedBlackP5 <- theoretical_mixtures_df_bottomed %>%
  filter(passage==5) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parent)), stat="identity", color="black", size=0.001) +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/abundancesPaintedBlackP5.png"), p_abundancesPaintedBlackP5, nrow=1.5, ncol=2)

# Redo painted ASV plots. -------------------------------------------------

# All painted ASV abundances, passage 3, replicate 1, rel ab 10^-4.
p_abundancesPaintedAllP3_rep1_relAb4 <- theoretical_mixtures_df %>%
  filter(passage==3) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parentOriginRep1_relAb4)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_abundancesPaintedAllP3_rep1_relAb4
save_plot(paste0(plotdir, "/abundancesPaintedAllP3_rep1Rel4.png"), 
          p_abundancesPaintedAllP3_rep1_relAb4, nrow=1.5, ncol=2) 

# All painted ASV abundances, passage 3, replicate 2, rel ab 10^-4.
p_abundancesPaintedAllP3_rep2_relAb4 <- theoretical_mixtures_df %>%
  filter(passage==3) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parentOriginRep2_relAb4)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_abundancesPaintedAllP3_rep2_relAb4
save_plot(paste0(plotdir, "/abundancesPaintedAllP3_rep2Rel4.png"), 
          p_abundancesPaintedAllP3_rep2_relAb4, nrow=1.5, ncol=2)

# All painted ASV abundances, passage 3, replicate 1, rel ab 10^-3.
p_abundancesPaintedAllP3_rep1_relAb3 <- theoretical_mixtures_df %>%
  filter(passage==3) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parentOriginRep1_relAb3)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_abundancesPaintedAllP3_rep1_relAb3
save_plot(paste0(plotdir, "/abundancesPaintedAllP3_rep1Rel3.png"), 
          p_abundancesPaintedAllP3_rep1_relAb3, nrow=1.5, ncol=2)

# All painted ASV abundances, passage 3, replicate 2, rel ab 10^-3.
p_abundancesPaintedAllP3_rep2_relAb3 <- theoretical_mixtures_df %>%
  filter(passage==3) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parentOriginRep2_relAb3)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_abundancesPaintedAllP3_rep2_relAb3
save_plot(paste0(plotdir, "/abundancesPaintedAllP3_rep2Rel3.png"),
          p_abundancesPaintedAllP3_rep2_relAb3, nrow=1.5, ncol=2)

# All painted ASV abundances, passage 5, replicate 1, rel ab 10^-4.
p_abundancesPaintedAllP5_rep1_relAb4 <- theoretical_mixtures_df %>%
  filter(passage==5) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parentOriginRep1_relAb4)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_abundancesPaintedAllP5_rep1_relAb4
save_plot(paste0(plotdir, "/abundancesPaintedAllP5_rep1Rel4.png"), 
          p_abundancesPaintedAllP5_rep1_relAb4, nrow=1.5, ncol=2) 

# All painted ASV abundances, passage 5, replicate 2, rel ab 10^-4.
p_abundancesPaintedAllP5_rep2_relAb4 <- theoretical_mixtures_df %>%
  filter(passage==5) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parentOriginRep2_relAb4)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_abundancesPaintedAllP5_rep2_relAb4
save_plot(paste0(plotdir, "/abundancesPaintedAllP5_rep2Rel4.png"), 
          p_abundancesPaintedAllP5_rep2_relAb4, nrow=1.5, ncol=2)

# All painted ASV abundances, passage 5, replicate 1, rel ab 10^-3.
p_abundancesPaintedAllP5_rep1_relAb3 <- theoretical_mixtures_df %>%
  filter(passage==5) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parentOriginRep1_relAb3)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_abundancesPaintedAllP5_rep1_relAb3
save_plot(paste0(plotdir, "/abundancesPaintedAllP5_rep1Rel3.png"), 
          p_abundancesPaintedAllP5_rep1_relAb3, nrow=1.5, ncol=2)

# All painted ASV abundances, passage 5, replicate 2, rel ab 10^-3.
p_abundancesPaintedAllP5_rep2_relAb3 <- theoretical_mixtures_df %>%
  filter(passage==5) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(parentOriginRep2_relAb3)), stat="identity") +
  xlab("Community") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", values=c("dark green", "light blue", "orange", "grey")) +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
p_abundancesPaintedAllP5_rep2_relAb3
save_plot(paste0(plotdir, "/abundancesPaintedAllP5_rep2Rel3.png"), 
          p_abundancesPaintedAllP5_rep2_relAb3, nrow=1.5, ncol=2)

# Line plots of ASV origin. -----------------------------------------------

# Combine replicates of parent origin into one column, calculate averages between lines for plotting.
parentOrigins_df <- rbind(
  theoretical_mixtures_df %>% 
    rename(parentOrigin_relAb3=parentOriginRep1_relAb3,
           parentOrigin_relAb4=parentOriginRep1_relAb4) %>% 
    mutate(parentOrigin_relAb3Rep=1, parentOrigin_relAb4Rep=1) %>% 
    select(!c(parentOriginRep2_relAb3, parentOriginRep2_relAb4)),
  theoretical_mixtures_df %>% 
    rename(parentOrigin_relAb3=parentOriginRep2_relAb3,
           parentOrigin_relAb4=parentOriginRep2_relAb4) %>% 
    mutate(parentOrigin_relAb3Rep=2, parentOrigin_relAb4Rep=2) %>% 
    select(!c(parentOriginRep1_relAb3, parentOriginRep1_relAb4))
        ) %>%
  filter(mixtureType=="actual") %>% 
  group_by(passage, combo, ratio, line, parentOrigin_relAb3, parentOrigin_relAb3Rep) %>% 
  mutate(relAbSum_relAb3=sum(rel_abundance)) %>% 
  ungroup() %>% 
  group_by(passage, combo, ratio, parentOrigin_relAb3, parentOrigin_relAb3Rep) %>% 
  mutate(relAbSumAvg_relAb3=mean(relAbSum_relAb3),
         relAbSumRange_relAb3=max(relAbSum_relAb3)-min(relAbSum_relAb3)) %>% 
  ungroup() %>% 
  group_by(passage, combo, ratio, line, parentOrigin_relAb4, parentOrigin_relAb4Rep) %>% 
  mutate(relAbSum_relAb4=sum(rel_abundance)) %>% 
  ungroup() %>% 
  group_by(passage, combo, ratio, parentOrigin_relAb4, parentOrigin_relAb4Rep) %>% 
  mutate(relAbSumAvg_relAb4=mean(relAbSum_relAb4),
         relAbSumRange_relAb4=max(relAbSum_relAb4)-min(relAbSum_relAb4)) %>% 
  ungroup()

p_ASVoriginP5_relAb3 <- parentOrigins_df %>% 
  filter(passage==5) %>%
  ggplot() +
  geom_line(aes(x=ratio, y=relAbSumAvg_relAb3, group=interaction(parentOrigin_relAb3,parentOrigin_relAb3Rep), 
                color=parentOrigin_relAb3)) +
  #geom_ribbon(aes(x=ratio, ymax=relAbSumAvg_relAb3+(relAbSumRange_relAb4/2), 
  #                  ymin=relAbSumAvg_relAb3-(relAbSumRange_relAb4/2), 
  #                  group=interaction(parentOrigin_relAb3,parentOrigin_relAb3Rep),
  #                  fill=parentOrigin_relAb3), color="black") +
  scale_color_brewer(palette="Dark2") +
  #scale_fill_brewer(palette="Dark2") +
  facet_wrap(~combo, nrow=2)
p_ASVoriginP5_relAb3
save_plot(paste0(plotdir, "/ASVoriginP5_relAb3.png"), 
          p_ASVoriginP5_relAb3, nrow=1.5, ncol=2)

# Plot # and rel ab of ASVs by origin. -------------------------------------------------------

# Summarize number of ASVs by each origin.
ASVoriginSummaryNumber <- theoretical_mixtures_df_categorized %>% 
  filter(passage==5 & category!="lowAbundance") %>% 
  group_by(combo, origin) %>% 
  summarize(nASVs = n_distinct(OTU)) %>% 
  ungroup()

# Plot # ASVs from each origin in each community.
p_ASVoriginNumber <- ASVoriginSummaryNumber %>% 
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
  geom_bar(aes(x=combo, y=nASVs, fill=origin), stat="identity") +
  xlab("Mixture") +
  ylab("# ASVs (1:1)") +
  scale_fill_manual(name="Origin", 
                    labels = c("A"="1", "B"="2", "both"="Both"),
                    values = c("A"="#BF2026", "B"="#ED7D2F", "both"="#62626A")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_ASVoriginNumber
save_plot(paste0(plotdir, "/ASVoriginNumber.pdf"), p_ASVoriginNumber, base_width=2.75, base_height=1.2)

# Summarize relative abundance of ASVs at 1:1 by each origin.
ASVoriginSummaryRelAb <- theoretical_mixtures_df_categorized %>% 
  filter(passage==5 & ratio=="1:1" & mixtureType=="actual" & !is.na(rel_abundance_old)) %>% 
  group_by(combo, OTU) %>% 
  summarize(origin, avgAbundance = sum(rel_abundance_old/3)) %>% 
  unique() %>% 
  arrange(-avgAbundance) %>% 
  ungroup()

# Plot rel ab of ASVs from each origin in each community.
p_ASVoriginRelAb <- ASVoriginSummaryRelAb %>% 
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
  geom_bar(aes(x=combo, y=avgAbundance, fill=origin), stat="identity", color="black", size=0.01) +
  xlab("Mixture") +
  ylab("Relative abundance") +
  scale_fill_manual(name="Origin", 
                    labels = c("A"="1", "B"="2", "both"="Both", "neither"="Neither"),
                    values = c("A"="#BF2026", "B"="#ED7D2F", "both"="#62626A", "neither"="#000000")) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_ASVoriginRelAb
save_plot(paste0(plotdir, "/ASVoriginRelAb.pdf"), p_ASVoriginRelAb, base_width=2.75, base_height=1.2)

ASVoriginLegend <- get_legend(p_ASVoriginRelAb)
save_plot(paste0(plotdir, "/ASVoriginLegend.pdf"), ASVoriginLegend, base_width=0.5, base_height=0.6)
