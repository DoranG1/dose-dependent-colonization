library(phyloseq)
library(tidyverse)
library(cowplot)
library(foreach)
library(RColorBrewer)
library(ape)
library(doParallel)
library(data.table)

# Allow parallel proccessing with foreach.
registerDoParallel(cores=6)

theme_set(theme_cowplot())
source("../../config/palettes/plotDefaults.R")

# Import data.
theoretical_mixtures_df_categorized <- fread("../mixtureDataframe.txt")

# Read in phylogenetic tree.
tree <- read.tree("../../workflow/out/e0017/DADA2_output/dsvs_msa.tree")

# Create phylogenetic distance matrix between ASVs from tree.
distMatrix <- cophenetic.phylo(tree)

# Turn distance matrix into dataframe.
distDf <- data.frame(distMatrix) %>% 
  mutate(OTU1 = row.names(distMatrix)) %>% 
  pivot_longer(cols=!OTU1, names_to="OTU2", values_to="distance") %>% 
  left_join(theoretical_mixtures_df_categorized %>% ungroup() %>%
              select(OTU, Family, OTUnum) %>% unique() %>%
              mutate(OTUID1 = paste0(Family, "-", OTUnum)) %>%
              rename(OTU1 = OTU) %>% select(OTU1, OTUID1)) %>%
  left_join(theoretical_mixtures_df_categorized %>% ungroup() %>%
              select(OTU, Family, OTUnum) %>% unique() %>%
              mutate(OTUID2 = paste0(Family, "-", OTUnum)) %>%
              rename(OTU2 = OTU) %>% select(OTU2, OTUID2))

combos <- c("XBA-XBB", "XCA-XCB", "XDA-XDB", "XFA-XFB", "XBA-XCA", "XCA-XDA", "XDA-XFA", "XFA-XBA")
passages <- c("3","5")
comboPassages <- expand.grid(combos,passages) %>% mutate(comboPassage=paste(Var1,Var2,sep="_"))
comboPassages <- comboPassages$comboPassage

# Analyze closest neighbors for all ASVs. ----------------------

# Create dataframe of phylogenetic distances between all DD ASVs in each passage/combo.
ASVDistances <- foreach(x = comboPassages, .combine="rbind") %do% {
  currCombo <- sub("_.*","",x)
  currPassage <- sub(".*_","",x)
  DDASVs <- theoretical_mixtures_df_categorized %>% 
    filter(combo==currCombo & passage==currPassage & category=="dose-dependent") %>%
    mutate(OTUID = paste0(Family, "-", OTUnum)) %>% 
    pull(OTUID) %>% 
    unique()
  # Calculate min distance to neighbor for each DD ASV.
  DDdistanceDf <- distDf %>%
    select(OTUID1, OTUID2, distance) %>% 
    filter(OTUID1 %in% DDASVs & OTUID2 %in% DDASVs & OTUID1!=OTUID2)
  minDDDistances <- DDdistanceDf %>% 
    group_by(OTUID1) %>% 
    filter(distance==min(distance)) %>% 
    ungroup() %>% 
    mutate(combo=currCombo, passage=currPassage, category="DD")
  OCASVs <- theoretical_mixtures_df_categorized %>% 
    filter(combo==currCombo & passage==currPassage & category=="overcolonizing") %>%
    mutate(OTUID = paste0(Family, "-", OTUnum)) %>% 
    pull(OTUID) %>% 
    unique()
  # Calculate min distance to neighbor for each OC ASV.
  OCdistanceDf <- distDf %>%
    select(OTUID1, OTUID2, distance) %>% 
    filter(OTUID1 %in% OCASVs & OTUID2 %in% OCASVs & OTUID1!=OTUID2)
  minOCDistances <- OCdistanceDf %>% 
    group_by(OTUID1) %>% 
    filter(distance==min(distance)) %>% 
    ungroup() %>% 
    mutate(combo=currCombo, passage=currPassage, category="OC")
  UCASVs <- theoretical_mixtures_df_categorized %>% 
    filter(combo==currCombo & passage==currPassage & category=="undercolonizing") %>%
    mutate(OTUID = paste0(Family, "-", OTUnum)) %>% 
    pull(OTUID) %>% 
    unique()
  # Calculate min distance to neighbor for each UC ASV.
  UCdistanceDf <- distDf %>%
    select(OTUID1, OTUID2, distance) %>% 
    filter(OTUID1 %in% UCASVs & OTUID2 %in% UCASVs & OTUID1!=OTUID2)
  minUCDistances <- UCdistanceDf %>% 
    group_by(OTUID1) %>% 
    filter(distance==min(distance)) %>% 
    ungroup() %>% 
    mutate(combo=currCombo, passage=currPassage, category="UC")
  # Calculate min distance to neighbor between OC and UC ASVs.
  OCUCdistanceDf <- distDf %>%
    select(OTUID1, OTUID2, distance) %>% 
    filter(OTUID1 %in% OCASVs & OTUID2 %in% UCASVs)
  minOCUCDistances <- OCUCdistanceDf %>% 
    group_by(OTUID1) %>% 
    filter(distance==min(distance)) %>% 
    ungroup() %>% 
    mutate(combo=currCombo, passage=currPassage, category="OC-UC")
  # Combine min distances from all categories.
  minDistances <- rbind(minDDDistances, minOCDistances, minUCDistances, minOCUCDistances)
}

# Permute ASVs in all categories and calculate nearest neighbor. ----------

# Permute ASVs in categories and calculate null distribution of closest phylogenetic distances.
# ASVDistancesPermuted <- foreach(x=comboPassages, .combine="rbind") %do% {
#   currCombo <- sub("_.*","",x)
#   currPassage <- sub(".*_","",x)
#   currAllASVsDf <- theoretical_mixtures_df_categorized %>% 
#     filter(passage==currPassage & combo==currCombo & category!="lowAbundance") %>% 
#     ungroup() %>% 
#     mutate(OTUID = paste0(Family, "-", OTUnum)) %>% 
#     select(OTU, OTUID, category) %>% 
#     unique() 
#   numDDASVs <- n_distinct(currAllASVsDf %>% filter(category=="dose-dependent") %>% pull(OTU))
#   numOCASVs <- n_distinct(currAllASVsDf %>% filter(category=="overcolonizing") %>% pull(OTU))
#   numUCASVs <- n_distinct(currAllASVsDf %>% filter(category=="undercolonizing") %>% pull(OTU))
#   categoriesVector <- c(rep("DD",numDDASVs), rep("OC", numOCASVs), rep("UC", numUCASVs))
#   numASVstotal <- length(categoriesVector)
#   minDistancesPermuted <- foreach(y=seq(1:1000), .combine="rbind") %do% {
#     # Permute DD, OC, and UC ASVs from full dataframe with no overlap.
#     currASVsDf <- currAllASVsDf %>% 
#       slice_sample(n = numASVstotal) %>% 
#       mutate(category = sample(categoriesVector))
#     # Calculate min distance to neighbor for each DD ASV.
#     DDASVsDf <- currASVsDf %>% 
#       filter(category=="DD")
#     currDDASVs <- DDASVsDf %>% pull(OTUID) %>% unique()
#     # Create a single dataframe with the distances between all simulated DD ASVs.
#     distanceDfCurrDDASVs <- distDf %>%
#       select(OTUID1, OTUID2, distance) %>% 
#       filter(OTUID1 %in% currDDASVs, OTUID2 %in% currDDASVs, OTUID1!=OTUID2)
#     minDDDistances <- distanceDfCurrDDASVs %>%
#       group_by(OTUID1) %>% 
#       filter(distance==min(distance)) %>% 
#       ungroup() %>% 
#       mutate(combo=currCombo, passage=currPassage, permutationReplicate=y, category="DD")
#     # Calculate min distance to neighbor for each OC ASV.
#     OCASVsDf <- currASVsDf %>% 
#       filter(category=="OC")
#     currOCASVs <- OCASVsDf %>% pull(OTUID) %>% unique()
#     # Create a single dataframe with the distances between all simulated OC ASVs.
#     distanceDfCurrOCASVs <- distDf %>%
#       select(OTUID1, OTUID2, distance) %>% 
#       filter(OTUID1 %in% currOCASVs, OTUID2 %in% currOCASVs, OTUID1!=OTUID2)
#     minOCDistances <- distanceDfCurrOCASVs %>%
#       group_by(OTUID1) %>% 
#       filter(distance==min(distance)) %>% 
#       ungroup() %>% 
#       mutate(combo=currCombo, passage=currPassage, permutationReplicate=y, category="OC")
#     # Calculate min distance to neighbor for each UC ASV.
#     UCASVsDf <- currASVsDf %>% 
#       filter(category=="UC")
#     currUCASVs <- UCASVsDf %>% pull(OTUID) %>% unique()
#     # Create a single dataframe with the distances between all simulated UC ASVs.
#     distanceDfCurrUCASVs <- distDf %>%
#       select(OTUID1, OTUID2, distance) %>% 
#       filter(OTUID1 %in% currUCASVs, OTUID2 %in% currUCASVs, OTUID1!=OTUID2)
#     minUCDistances <- distanceDfCurrUCASVs %>%
#       group_by(OTUID1) %>% 
#       filter(distance==min(distance)) %>% 
#       ungroup() %>% 
#       mutate(combo=currCombo, passage=currPassage, permutationReplicate=y, category="UC")
#     # Calculate min distance to neighbor between OC and UC ASVs.
#     # Create a single dataframe with the distances between all simulated UC ASVs.
#     distanceDfCurrOCUCASVs <- distDf %>%
#       select(OTUID1, OTUID2, distance) %>% 
#       filter(OTUID1 %in% currOCASVs, OTUID2 %in% currUCASVs)
#     minOCUCDistances <- distanceDfCurrOCUCASVs %>%
#       group_by(OTUID1) %>% 
#       filter(distance==min(distance)) %>% 
#       ungroup() %>% 
#       mutate(combo=currCombo, passage=currPassage, permutationReplicate=y, category="OC-UC")
#     # Combine min distances from all categories.
#     minDistances <- rbind(minDDDistances, minOCDistances, minUCDistances, minOCUCDistances)
#   }
# }
#write.table(ASVDistancesPermuted, "ASVDistancesPermuted.txt", quote=FALSE, row.names=FALSE, sep="\t")

ASVDistancesPermuted <- fread("ASVDistancesPermuted.txt")

# Plot categories separately. ---------------------------------------------

# Plot closest phylogenetic distances for DD ASVs, all points.
p_closestDDdistancesAll <- closestDDASVDistances %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_jitter(aes(x=combo, y=distance), width=0.15) +
  ylim(0,1.15) +
  DEFAULTS.THEME_PRINT
p_closestDDdistancesAll
save_plot("out/closestDDdistancesAll.png", p_closestDDdistancesAll)

# Plot closest phylogenetic distances for DD ASVs, points colored by single/pair.
p_closestDDdistancesPaired <- closestDDASVDistances %>% 
  filter(passage==5) %>% 
  ungroup() %>% 
  group_by(combo, distance) %>% 
  mutate(type = ifelse(n()==1, "single", "double")) %>% 
  slice(1) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(aes(x=combo, y=distance, color=type)) +
  ylim(0,1.15) +
  DEFAULTS.THEME_PRINT
p_closestDDdistancesPaired
save_plot("out/closestDDdistancesPaired.png", p_closestDDdistancesPaired)

# Plot closest phylogenetic distances for OC ASVs, all points.
p_closestOCdistancesAll <- closestOCASVDistances %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_jitter(aes(x=combo, y=distance), width=0.15) +
  ylim(0,1.15) +
  DEFAULTS.THEME_PRINT
p_closestOCdistancesAll
save_plot("out/closestOCdistancesAll.png", p_closestOCdistancesAll)

# Plot closest phylogenetic distances for OC ASVs, points colored by single/pair.
p_closestOCdistancesPaired <- closestOCASVDistances %>% 
  filter(passage==5) %>% 
  ungroup() %>% 
  group_by(combo, distance) %>% 
  mutate(type = ifelse(n()==1, "single", "double")) %>% 
  slice(1) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(aes(x=combo, y=distance, color=type)) +
  ylim(0,1.15) +
  DEFAULTS.THEME_PRINT
p_closestOCdistancesPaired
save_plot("out/closestOCdistancesPaired.png", p_closestOCdistancesPaired)

# Plot closest phylogenetic distances for UC ASVs, all points.
p_closestUCdistancesAll <- closestUCASVDistances %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_jitter(aes(x=combo, y=distance), width=0.15) +
  ylim(0,1.15) +
  DEFAULTS.THEME_PRINT
p_closestUCdistancesAll
save_plot("out/closestUCdistancesAll.png", p_closestUCdistancesAll)

# Plot closest phylogenetic distances for UC ASVs, points colored by single/pair.
p_closestUCdistancesPaired <- closestUCASVDistances %>% 
  filter(passage==5) %>% 
  ungroup() %>% 
  group_by(combo, distance) %>% 
  mutate(type = ifelse(n()==1, "single", "double")) %>% 
  slice(1) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(aes(x=combo, y=distance, color=type)) +
  ylim(0,1.15) +
  DEFAULTS.THEME_PRINT
p_closestUCdistancesPaired
save_plot("out/closestUCdistancesPaired.png", p_closestUCdistancesPaired)


# Plot closest phylogenetic distances for UC ASVs, all points.
p_closestOCUCdistancesAll <- closestOCUCASVDistances %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_jitter(aes(x=combo, y=distance), width=0.15) +
  ylim(0,1.2) +
  DEFAULTS.THEME_PRINT
p_closestOCUCdistancesAll
save_plot("out/closestOCUCdistancesAll.png", p_closestOCUCdistancesAll)

# Plot closest phylogenetic distances for UC ASVs, points colored by single/pair.
p_closestOCUCdistancesPaired <- closestOCUCASVDistances %>% 
  filter(passage==5) %>% 
  ungroup() %>% 
  group_by(combo, distance) %>% 
  mutate(type = ifelse(n()==1, "single", "double")) %>% 
  slice(1) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(aes(x=combo, y=distance, color=type)) +
  ylim(0,1.2) +
  DEFAULTS.THEME_PRINT
p_closestOCUCdistancesPaired
save_plot("out/closestOCUCdistancesPaired.png", p_closestOCUCdistancesPaired)

# Plot all categories together. -------------------------------------------

# Identify set of ASVs that are ever present above low abundance in data.
highAbundanceASVs <- theoretical_mixtures_df_categorized %>% 
  ungroup() %>% 
  #filter(category!="lowAbundance") %>% 
  filter(category %in% c("dose-dependent","overcolonizing","undercolonizing")) %>% 
  pull(OTU) %>% 
  unique()

# Create dataframe with combined phylogenetic distances from different comparisons.
closestDistancesCombined <- rbind(
  closestDDASVDistances %>% 
    mutate(category1="dose-dependent", category2="dose-dependent", type="DD"),
  closestOCASVDistances %>% 
    mutate(category1="overcolonizing", category2="overcolonizing", type="OC"),
  closestUCASVDistances %>% 
    mutate(category1="undercolonizing", category2="undercolonizing", type="UC"),
  closestOCUCASVDistances %>% 
    mutate(type="OC-UC"),
  # Add null distribution from all ASVs.
  distDf %>% 
    filter(OTU1 %in% highAbundanceASVs & OTU2 %in% highAbundanceASVs & OTU1!=OTU2) %>% 
    mutate(combo=NA, passage=NA, Family2=NA, OTUnum2=NA, Family1=NA, OTUnum1=NA,
           OTUID1=NA, category1=NA, category2=NA, type="All") %>% 
    group_by(OTU1) %>% 
    slice_min(distance, n=1)
  ) %>% 
  ungroup() %>% 
  mutate(type = fct_relevel(type, c("All", "DD","OC","UC","OC-UC")))

p_closestDistancesCombined <- closestDistancesCombined %>% 
  ggplot() +
  #geom_jitter(aes(x=type, y=distance), width=0.15) +
  geom_violin(aes(x=type, y=distance), size=0.25) +
  geom_boxplot(aes(x=type, y=distance), size=0.25, width=0.05, outlier.shape=NA) +
  #ylim(0,1.2) +
  DEFAULTS.THEME_PRINT
p_closestDistancesCombined
save_plot("out/closestDistancesCombined.png", p_closestDistancesCombined)

# Create dataframe with actual and permuted closest distances.
closestDistancesPermuted <- rbind(
  ASVDistances %>% mutate(permutationReplicate=NA, group="Actual"),
  ASVDistancesPermuted %>% mutate(group="Permuted")
  ) %>% 
  mutate(group = fct_relevel(group, c("Permuted","Actual")),
         category = case_when(category=="DD" ~ "DD",
                              category=="OC" ~ "S",
                              category=="UC" ~ "W",
                              category=="OC-UC" ~ "S-W"),
         category = fct_relevel(category, c("DD","S","W","S-W")))

p_closestDistancesPermutedDD <- closestDistancesPermuted %>% 
  #filter(category!="S") %>% 
  filter(category=="DD") %>% 
  ggplot() +
  geom_violin(aes(x=group, y=distance, fill=group), size=0.25) +
  geom_boxplot(aes(x=group, y=distance, group=group), 
               size=0.25, width=0.1, outlier.shape=NA) +
  scale_fill_manual("Type", values = c("#707070","#B8DBD9")) +
  xlab("Distribution") +
  ylab("Phylogenetic distance") +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        legend.position="top") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_closestDistancesPermutedDD
save_plot("out/closestDistancesPermutedDD.pdf", p_closestDistancesPermutedDD, base_width=1.75, base_height=1.5)

p_closestDistancesPermutedOther <- closestDistancesPermuted %>% 
  filter(category!="DD") %>% 
  #filter(category=="DD") %>% 
  ggplot() +
  geom_violin(aes(x=category, y=distance, fill=group), 
              size=0.25, position=position_dodge(width=0.75)) +
  geom_boxplot(aes(x=category, y=distance, group=interaction(category, group)), 
               size=0.25, width=0.1, outlier.shape=NA, position=position_dodge(width=0.75)) +
  scale_fill_manual("Type", values = c("#707070","#B8DBD9")) +
  xlab("Pattern") +
  ylab("Phylogenetic distance") +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        legend.position="top") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_closestDistancesPermutedOther
save_plot("out/closestDistancesPermutedOther.pdf", p_closestDistancesPermutedOther, base_width=2.5, base_height=1.5)

p_closestDistancesPermutedAll <- closestDistancesPermuted %>% 
  ggplot() +
  geom_violin(aes(x=category, y=distance, fill=group), 
              size=0.25, position=position_dodge(width=0.75)) +
  geom_boxplot(aes(x=category, y=distance, group=interaction(category, group)), 
               size=0.25, width=0.1, outlier.shape=NA, position=position_dodge(width=0.75)) +
  scale_fill_manual("Type", values = c("#707070","#B8DBD9")) +
  xlab("Pattern") +
  ylab("Phylogenetic distance") +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        legend.position="top") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_closestDistancesPermutedAll
save_plot("out/closestDistancesPermutedAll.pdf", p_closestDistancesPermutedAll, base_width=2.75, base_height=1.5)

closestDistancesLegend <- get_legend(p_closestDistancesPermuted)
save_plot("out/closestDistancesLegend.pdf", closestDistancesLegend, base_width=1, base_height=0.25)

# Check significance of distributions.
DDsig <- t.test(closestDistancesPermuted %>% filter(category=="DD" & group=="Actual") %>% pull(distance),
                closestDistancesPermuted %>% filter(category=="DD" & group=="Permuted") %>% pull(distance))
Wsig <- t.test(closestDistancesPermuted %>% filter(category=="W" & group=="Actual") %>% pull(distance),
                closestDistancesPermuted %>% filter(category=="W" & group=="Permuted") %>% pull(distance))
SWsig <- t.test(closestDistancesPermuted %>% filter(category=="S-W" & group=="Actual") %>% pull(distance),
                closestDistancesPermuted %>% filter(category=="S-W" & group=="Permuted") %>% pull(distance))
