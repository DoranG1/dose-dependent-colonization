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

# Plot family abundances, example plot with one theoretical mixture.
p_theoreticalFamilyAbundancesExample <- theoretical_mixtures_df_categorized %>%
  filter(mixtureType=="theoretical" & combo=="XFA-XBA" & inoculationReplicate==1 & passage==5) %>%
  #mutate(combo="Theoretical\nXFA-XBA") %>% 
  mutate(combo="Theoretical\nD1/A1") %>% 
  # Filter relative abundance for legend sizing.
  #filter(family_rel_abundance > 0.01) %>% 
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
save_plot(paste0(plotdir, "/theoreticalExampleLegend.pdf"), theoreticalExampleLegend, base_width=1.1, base_height=1.5)

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
  #theme(legend.position="none") +
  facet_grid(~combo~ratio, scales="free", switch="y") +
  DEFAULTS.THEME_PRINT +
  theme(strip.placement="outside",
        strip.text.y.left=element_text(angle=0))
p_allActualAbundances
save_plot(paste0(plotdir, "/allActualAbundances.pdf"), p_allActualAbundances, 
          base_width=7, base_height=6.5)

e0017Legend <- get_legend(p_allActualAbundances)
save_plot(paste0(plotdir, "/e0017Legend.pdf"), e0017Legend, base_width=5, base_height=0.75)
