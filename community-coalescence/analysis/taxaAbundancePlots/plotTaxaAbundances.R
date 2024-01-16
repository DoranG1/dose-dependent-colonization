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
library(data.table)

setwd("/oak/stanford/groups/relman/users/dorang/210827-e0017-16S")
source("analysis/generateMixturesDataframe.R")
plotdir <- "analysis/taxaAbundancePlots"

# Import family codes.
familyCoding <- fread("config/familyCodingFull.txt")
theoretical_mixtures_df_categorized <- theoretical_mixtures_df_categorized %>% 
  left_join(familyCoding %>% rename(Family=family))

# Create version of palette to use with family codes.
paletteCodes <- left_join(palette, familyCoding %>% rename(taxaShort=family))
paletteVectorCodes <- paletteCodes$hex
names(paletteVectorCodes) <- paletteCodes$code

# Family abundance plots ---------------------------------------------------------

# Plot abundances for all mixtures, passage 3.
p_allFamilyAbundancesP3 <- theoretical_mixtures_df_bottomed %>%
  filter(passage==3) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(taxa)), stat="identity") +
  scale_fill_manual(values=taxaPaletteList) +
  xlab("Community") +
  ylab("Family relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/allFamilyAbundancesP3.png"), p_allFamilyAbundancesP3, nrow=1.5, ncol=1.75)

# Plot abundances for all mixtures, passage 5.
p_allFamilyAbundancesP5 <- theoretical_mixtures_df_bottomed %>%
  filter(passage==5) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(taxa)), stat="identity") +
  scale_fill_manual(values=taxaPaletteList) +
  xlab("Community") +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/allFamilyAbundancesP5.png"), p_allFamilyAbundancesP5, nrow=1.5, ncol=1.75)

# Plot abundances for theoretical mixtures.
p_theoreticalFamilyAbundances <- theoretical_mixtures_df_bottomed %>%
  filter(mixtureType=="theoretical") %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(taxa)), stat="identity") +
  scale_fill_manual(values=taxaPaletteList) +
  xlab("Community") +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/theoreticalFamilyAbundances.png"),
          p_theoreticalFamilyAbundances, nrow=1, ncol=1.75)

# Plot family abundances, example plot with one theoretical mixture.
p_theoreticalFamilyAbundancesExample <- theoretical_mixtures_df_categorized %>%
  filter(mixtureType=="theoretical" & combo=="XFA-XBA" & inoculationReplicate==1 & passage==5) %>%
  #mutate(combo="Theoretical\nXFA-XBA") %>% 
  mutate(combo="Theoretical\nD1/A1") %>% 
  # Filter relative abundance for legend sizing.
  filter(rel_abundance_old > 0.005) %>% 
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance_old, fill=Family), stat="identity", color="black", size=0.05) +
  scale_fill_manual(values=paletteVector) +
  xlab("Mixture ratio") +
  ylab("Relative abundance") +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  guides(fill=guide_legend(ncol=1)) +
  #theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_grid(~combo) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=5.5))
p_theoreticalFamilyAbundancesExample
save_plot(paste0(plotdir, "/theoreticalFamilyAbundancesExample.pdf"),
          p_theoreticalFamilyAbundancesExample, base_width=2.8, base_height=1.8)

theoreticalExampleLegend <- get_legend(p_theoreticalFamilyAbundancesExample)
save_plot(paste0(plotdir, "/theoreticalExampleLegend.pdf"), theoreticalExampleLegend, base_width=1.1, base_height=1.75)

# Plot family abundances, example plot with one theoretical and actual mixture.
p_familyAbundancesExample <- theoretical_mixtures_df_bottomed %>%
  filter(combo=="XDA-XDB" & inoculationReplicate==2) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(taxa)), stat="identity") +
  scale_fill_manual(values=taxaPaletteList) +
  xlab("Community") +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/familyAbundancesExample.png"),
          p_familyAbundancesExample, nrow=1, ncol=1)

# Plot abundances for all actual mixtures.
p_allActualAbundances <- theoretical_mixtures_df_categorized %>% 
  filter(mixtureType=="actual" & passage==5) %>% 
  rbind(theoretical_mixtures_df_categorized %>% 
          filter(passage==5 & ratio %in% c("1:0","0:1"))) %>%  
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
  ggplot() +
  geom_bar(aes(x=inoculationReplicate, y=rel_abundance_old, fill=Family), stat="identity", color="black", size=0.01) +
  scale_fill_manual(values=paletteVector) +
  xlab("Replicate") +
  ylab("Relative abundance") +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        legend.position="bottom") +
  theme(legend.position="none") +
  facet_grid(~combo~ratio, scales="free", switch="y") +
  DEFAULTS.THEME_PRINT +
  theme(strip.placement="outside",
        strip.text.y.left=element_text(angle=0))
p_allActualAbundances
save_plot(paste0(plotdir, "/allActualAbundances.pdf"), p_allActualAbundances, 
          base_width=7, base_height=6.5)

e0017Legend <- get_legend(p_allActualAbundances)
save_plot(paste0(plotdir, "/e0017Legend.pdf"), e0017Legend, base_width=5, base_height=0.75)

# Plot abundances for all actual mixtures.
p_allActualAbundancesXDAXDB <- theoretical_mixtures_df_categorized %>% 
  filter(mixtureType=="actual" & passage==5) %>%
  rbind(theoretical_mixtures_df_categorized %>% 
          filter(passage==5 & ratio %in% c("1:0","0:1"))) %>%  
  filter(combo=="XDA-XDB") %>% 
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
  ggplot() +
  geom_bar(aes(x=inoculationReplicate, y=rel_abundance_old, fill=Family), stat="identity", color="black", size=0.01) +
  scale_fill_manual(values=paletteVector) +
  xlab("Mixture ratio") +
  ylab("Relative abundance") +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        legend.position="bottom") +
  theme(legend.position="none") +
  facet_grid(~combo~ratio, scales="free", switch="x") +
  ggtitle("C1/C2") +
  theme(plot.title=element_text(size=5.5, face="plain", hjust=0.5)) +
  DEFAULTS.THEME_PRINT +
  theme(strip.placement="outside",
        strip.text.y.left=element_text(angle=0))
p_allActualAbundancesXDAXDB
save_plot(paste0(plotdir, "/allActualAbundancesXDAXDB.pdf"), p_allActualAbundancesXDAXDB, 
          base_width=7, base_height=1.75)

p_actualAbundancesNoParentsXBAXBB <- theoretical_mixtures_df_categorized %>% 
  filter(mixtureType=="actual" & passage==5 & combo=="XBA-XBB") %>% 
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance_old, fill=Family), stat="identity", color="black", size=0.01) +
  scale_fill_manual(values=paletteVector) +
  xlab("Ratio") +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  facet_wrap(~inoculationReplicate, nrow=3) +
  DEFAULTS.THEME_PRINT +
  theme(strip.text.x=element_blank())
p_actualAbundancesNoParentsXBAXBB
save_plot(paste0(plotdir, "/actualAbundancesNoParentsXBAXBB.pdf"), 
          p_actualAbundancesNoParentsXBAXBB, base_width=1.25, base_height=2.45)

p_actualAbundancesParent1XBAXBB <- theoretical_mixtures_df_categorized %>% 
  filter(combo=="XBA-XBB" & ratio=="1:0" & passage==5) %>% 
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance_old, fill=Family), stat="identity", color="black", size=0.01) +
  scale_fill_manual(values=paletteVector) +
  xlab("Ratio") +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x=element_blank()) +
  facet_wrap(~inoculationReplicate, nrow=2) +
  DEFAULTS.THEME_PRINT +
  theme(strip.text.x=element_blank())
p_actualAbundancesParent1XBAXBB
save_plot(paste0(plotdir, "/actualAbundancesParent1XBAXBB.pdf"), 
          p_actualAbundancesParent1XBAXBB, base_width=0.75, base_height=2.2)

p_actualAbundancesParent2XBAXBB <- theoretical_mixtures_df_categorized %>% 
  filter(combo=="XBA-XBB" & ratio=="0:1" & passage==5) %>% 
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance_old, fill=Family), stat="identity", color="black", size=0.01) +
  scale_fill_manual(values=paletteVector) +
  xlab("Ratio") +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank()) +
  facet_wrap(~inoculationReplicate, nrow=2) +
  DEFAULTS.THEME_PRINT +
  theme(strip.text.x=element_blank())
p_actualAbundancesParent2XBAXBB
save_plot(paste0(plotdir, "/actualAbundancesParent2XBAXBB.pdf"), 
          p_actualAbundancesParent2XBAXBB, base_width=0.4, base_height=2.2)
  

# Genus abundance plots --------------------------------------------------

# Plot genus abundances for all mixtures, passage 3 with black lines.
p_allGenusAbundancesP3 <- theoretical_mixtures_df_bottomed %>%
  filter(passage==3) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(taxa)), color="black", size=0.001, stat="identity") +
  scale_fill_manual(values=taxaPaletteList) +
  xlab("Community") +
  ylab("Genus relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/allGenusAbundancesP3.png"), p_allGenusAbundancesP3, nrow=1.5, ncol=1.75)

# Plot genus abundances for all mixtures, passage 5 with black lines.
p_allGenusAbundancesP5 <- theoretical_mixtures_df_bottomed %>%
  filter(passage==5) %>%
  ggplot() +
  geom_bar(aes(x=ratio, y=rel_abundance, fill=factor(taxa)), color="black", size=0.001, stat="identity") +
  scale_fill_manual(values=taxaPaletteList) +
  xlab("Community") +
  ylab("Genus relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  facet_grid(~line~combo, scales="free") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10))
save_plot(paste0(plotdir, "/allGenusAbundancesP5.png"), p_allGenusAbundancesP5, nrow=1.5, ncol=1.75)
