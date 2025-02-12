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

# Reorder taxa by full taxonomy for correct color ordering.
theoretical_mixtures_df_categorized <- theoretical_mixtures_df_categorized %>% 
  left_join(taxaPalette %>% select(Family, taxa), by="Family")

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
                             "overcolonizing"="#CCA245",
                             "undercolonizing"="#8FB4DD",
                             "dose-dependent"="#762523")) +
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
#save_plot(paste0(plotdirParentDiff, "/ASVcategoriesNumberRBY.pdf"), p_ASVcategoriesNumberRBY, 
#          base_width=3.4, base_height=1.2)
save_plot(paste0(plotdirParentDiff, "/ASVcategoriesNumberRBY.png"), p_ASVcategoriesNumberRBY, 
          base_width=3.24, base_height=1.2)

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
