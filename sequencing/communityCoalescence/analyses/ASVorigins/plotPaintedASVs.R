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
