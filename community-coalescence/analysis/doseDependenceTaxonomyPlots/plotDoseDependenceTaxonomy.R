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
library(reshape2)
library(patchwork)
library(plyr)
library(data.table)

source("analysis/generateMixturesDataframe.R")
plotdir <- "analysis/doseDependenceTaxonomyPlots"

# Get relative abundances of dose-dependent ASVs at key mixture ratios.
lowestTaxaRelAb <- theoretical_mixtures_df_categorized %>% 
  filter(category=="dose-dependent" & 
           (ratio %in% c("1:0","0:1") | (ratio=="1:1" & mixtureType=="actual")) & !is.na(rel_abundance_old)) %>% 
  group_by(passage, combo, ratio, Phylum, Class, Order, Family, Genus, OTUnum) %>% 
  dplyr::summarize(avgAbundance=ifelse(ratio=="1:1", sum(rel_abundance_old/3) , sum(rel_abundance_old/2))) %>% 
  unique()

# Quantify dose-dependence among taxa. ------------------------------------

# Multiple dose-dependent taxa, in opposite directions, at the same taxonomic level.

# Count number of cases in which opposite-direction dose-dependent ASVs within same taxa at various levels.
# Filter dataframe to just dose-dependent taxa, add marker for opposite-direction ASVs.
taxaDoseDependence_ASV <- theoretical_mixtures_df_categorized %>% 
  filter(category=="dose-dependent") %>% 
  mutate(direction=a_doseDifference > 0)
# Add marker for presence of opposite-direction ASVs at each taxonomic level.
taxaDoseDependence_ASV <- join_all(list(
  taxaDoseDependence_ASV %>% 
    group_by(passage, combo, Phylum) %>%
    dplyr::summarize(Class, Order, Family, Genus, OTUnum, oppositeDirectionsPhylum = TRUE %in% direction & FALSE %in% direction &
                       abs(sum(a_doseDifference)) <= 0.5) %>% 
    unique(),
  taxaDoseDependence_ASV %>% 
    group_by(passage, combo, Phylum, Class) %>%
    dplyr::summarize(Order, Family, Genus, OTUnum, oppositeDirectionsClass = TRUE %in% direction & FALSE %in% direction &
                       abs(sum(a_doseDifference)) <= 0.5) %>% 
    unique(),
  taxaDoseDependence_ASV %>%  
    group_by(passage, combo, Phylum, Class, Order) %>%
    dplyr::summarize(Family, Genus, OTUnum, oppositeDirectionsOrder = TRUE %in% direction & FALSE %in% direction & 
                       abs(sum(a_doseDifference)) <= 0.5) %>% 
    unique(),
  taxaDoseDependence_ASV %>% 
    group_by(passage, combo, Phylum, Class, Order, Family) %>%
    dplyr::summarize(Genus, OTUnum, oppositeDirectionsFamily = TRUE %in% direction & FALSE %in% direction &
                       abs(sum(a_doseDifference)) <= 0.5) %>% 
    unique(),
  taxaDoseDependence_ASV %>% 
    group_by(passage, combo, Phylum, Class, Order, Family, Genus) %>% 
    dplyr::summarize(OTUnum, oppositeDirectionsGenus = 
                       TRUE %in% direction & FALSE %in% direction & !is.na(Genus) &
                       abs(sum(a_doseDifference)) <= 0.5) %>% 
    unique()
  ), type='left') %>% 
  filter(!(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>% 
  # Set lowest rank at which opposite direction dose dependence occurs
  mutate(oppositeDirectionsPhylum = oppositeDirectionsPhylum & !oppositeDirectionsClass & 
           !oppositeDirectionsOrder & !oppositeDirectionsFamily & !oppositeDirectionsGenus,
         oppositeDirectionsClass = oppositeDirectionsClass & !oppositeDirectionsOrder & 
           !oppositeDirectionsFamily & !oppositeDirectionsGenus,
         oppositeDirectionsOrder = oppositeDirectionsOrder & !oppositeDirectionsFamily &
           !oppositeDirectionsGenus,
         oppositeDirectionsFamily = oppositeDirectionsFamily & !oppositeDirectionsGenus) %>% 
  group_by(passage, combo, Phylum, Class, Order, Family, Genus, OTUnum) %>% 
  dplyr::summarize(oppositeDirectionsPhylum, oppositeDirectionsClass, oppositeDirectionsOrder, 
                   oppositeDirectionsFamily, oppositeDirectionsGenus,
                   lowestTaxonomicLevel = case_when(
                     TRUE %in% oppositeDirectionsPhylum ~ "Phylum",
                     TRUE %in% oppositeDirectionsClass ~ "Class",
                     TRUE %in% oppositeDirectionsOrder ~ "Order",
                     TRUE %in% oppositeDirectionsFamily ~ "Family",
                     TRUE %in% oppositeDirectionsGenus ~ "Genus"),
                   lowestTaxonomicLevel = replace(lowestTaxonomicLevel, is.na(lowestTaxonomicLevel), "None"))
# Get relative abundances of dose-dependent ASVs at key mixture ratios.
lowestTaxaRelAb <- theoretical_mixtures_df_categorized %>% 
  filter(category=="dose-dependent" & 
           (ratio %in% c("1:0","0:1") | (ratio=="1:1" & mixtureType=="actual")) & !is.na(rel_abundance_old)) %>% 
  group_by(passage, combo, ratio, Phylum, Class, Order, Family, Genus, OTUnum) %>% 
  dplyr::summarize(avgAbundance=ifelse(ratio=="1:1", sum(rel_abundance_old/3) , sum(rel_abundance_old/2))) %>% 
  unique()
# Bind relative abundances to dose-dependent ASVs.
taxaDoseDependence_ASV <- taxaDoseDependence_ASV %>% right_join(lowestTaxaRelAb) %>% 
  ungroup() %>% 
  mutate(lowestTaxonomicLevel = fct_relevel(
    lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus"))) %>% 
  group_by(passage, combo, ratio, lowestTaxonomicLevel) %>% 
  arrange(-avgAbundance) %>% 
  ungroup()

# Count instances of each level of lowest shared taxa.
lowestTaxaOpposite <- taxaDoseDependence_ASV %>%
  select(-c(ratio, avgAbundance)) %>% 
  group_by(passage, combo, lowestTaxonomicLevel) %>% 
  unique() %>% 
  dplyr::summarize(count=n())

# Quantify dose-dependence among taxa, all ASV categories. --------

# Count number of cases in which opposite-direction dose-dependent ASVs within same taxa at various levels, 
# inlcuding all ASVs in each category.
# Add marker for opposite-direction dose-dependent ASVs.
# taxaDoseDependence_ASV_allCategories <- theoretical_mixtures_df_categorized %>%
#   filter(category!="lowAbundance")
#   #mutate(direction = case_when(category=="dose-dependent" & a_doseDifference > 0 ~ "positive",
#   #                             category=="dose-dependent" & a_doseDifference <= 0 ~ "negative",
#   #                             category!="dose-dependent" ~ "NA"))
# Add marker for presence of opposite-direction ASVs at each taxonomic level.
taxaDoseDependence_ASV_allCategories <- join_all(list(
  theoretical_mixtures_df_categorized %>%
    filter(category!="lowAbundance") %>% 
    group_by(passage, combo, Phylum) %>%
    dplyr::summarize(Class, Order, Family, Genus, OTUnum, category, 
                     oppositeDirectionsPhylum = "dose-dependent" %in% category & n_distinct(OTUnum) > 1 & 
                       abs(sum(unique(a_doseDifference))) <= 0.09) %>% 
    unique(),
  theoretical_mixtures_df_categorized %>%
    filter(category!="lowAbundance") %>% 
    group_by(passage, combo, Phylum, Class) %>%
    dplyr::summarize(Order, Family, Genus, OTUnum, category, 
                     oppositeDirectionsClass = "dose-dependent" %in% category & n_distinct(OTUnum) > 1 &
                       abs(sum(unique(a_doseDifference))) <= 0.09) %>% 
    unique(),
  theoretical_mixtures_df_categorized %>%
    filter(category!="lowAbundance") %>%
    group_by(passage, combo, Phylum, Class, Order) %>%
    dplyr::summarize(Family, Genus, OTUnum, category, 
                     oppositeDirectionsOrder = "dose-dependent" %in% category & n_distinct(OTUnum) > 1 &
                       abs(sum(unique(a_doseDifference))) <= 0.09) %>% 
    unique(),
  theoretical_mixtures_df_categorized %>%
    filter(category!="lowAbundance") %>% 
    group_by(passage, combo, Phylum, Class, Order, Family) %>%
    dplyr::summarize(Genus, OTUnum, category, 
                     oppositeDirectionsFamily = "dose-dependent" %in% category & n_distinct(OTUnum) > 1 &
                       abs(sum(unique(a_doseDifference))) <= 0.09) %>% 
    unique(),
  theoretical_mixtures_df_categorized %>%
    filter(category!="lowAbundance") %>% 
    group_by(passage, combo, Phylum, Class, Order, Family, Genus) %>% 
    dplyr::summarize(OTUnum, category, 
                     oppositeDirectionsGenus = "dose-dependent" %in% category & n_distinct(OTUnum) > 1 &
                       abs(sum(unique(a_doseDifference))) <= 0.09 & !is.na(Genus)) %>% 
    unique()
), type='left') %>% 
  filter(category=="dose-dependent") %>% 
  filter(!(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>% 
  # Set lowest rank at which opposite direction dose dependence occurs
  mutate(oppositeDirectionsPhylum = oppositeDirectionsPhylum & !oppositeDirectionsClass & 
           !oppositeDirectionsOrder & !oppositeDirectionsFamily & !oppositeDirectionsGenus,
         oppositeDirectionsClass = oppositeDirectionsClass & !oppositeDirectionsOrder & 
           !oppositeDirectionsFamily & !oppositeDirectionsGenus,
         oppositeDirectionsOrder = oppositeDirectionsOrder & !oppositeDirectionsFamily &
           !oppositeDirectionsGenus,
         oppositeDirectionsFamily = oppositeDirectionsFamily & !oppositeDirectionsGenus) %>% 
  group_by(passage, combo, Phylum, Class, Order, Family, Genus, OTUnum) %>% 
  dplyr::summarize(oppositeDirectionsPhylum, oppositeDirectionsClass, oppositeDirectionsOrder, 
                   oppositeDirectionsFamily, oppositeDirectionsGenus,
                   lowestTaxonomicLevel = case_when(
                     TRUE %in% oppositeDirectionsPhylum ~ "Phylum",
                     TRUE %in% oppositeDirectionsClass ~ "Class",
                     TRUE %in% oppositeDirectionsOrder ~ "Order",
                     TRUE %in% oppositeDirectionsFamily ~ "Family",
                     TRUE %in% oppositeDirectionsGenus ~ "Genus"),
                   lowestTaxonomicLevel = replace(lowestTaxonomicLevel, is.na(lowestTaxonomicLevel), "None"))
# Bind relative abundances to dose-dependent ASVs.
taxaDoseDependence_ASV_allCategories <- taxaDoseDependence_ASV_allCategories %>% right_join(lowestTaxaRelAb) %>% 
  ungroup() %>% 
  mutate(lowestTaxonomicLevel = fct_relevel(
    lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus"))) %>% 
  group_by(passage, combo, ratio, lowestTaxonomicLevel) %>% 
  arrange(-avgAbundance) %>% 
  ungroup()

# Count instances of each level of lowest shared taxa.
lowestTaxaOpposite_allCategories <- taxaDoseDependence_ASV_allCategories %>%
  select(-c(ratio, avgAbundance)) %>% 
  group_by(passage, combo, lowestTaxonomicLevel) %>% 
  unique() %>% 
  dplyr::summarize(count=n())

# Quantify dose-dependence complementation, all subsets of all categories. --------

DDASVs_comboPassages <- theoretical_mixtures_df_categorized %>% 
  filter(category=="dose-dependent") %>% 
  mutate(ID = paste0(Family, "-", OTUnum, ".", combo, ":", passage)) %>% 
  pull(ID) %>% 
  unique()

DDcomplementationDf <- foreach(x=DDASVs_comboPassages, .combine="rbind") %dopar% {
  currOTUID <- sub("\\..*","",x)
  currCombo <- sub(":.*","",sub(".*\\.","",x))
  currPassage <- sub(".*:","",x)
  currDf <- theoretical_mixtures_df_categorized %>% 
    filter(combo==currCombo & passage==currPassage & paste0(Family, "-", OTUnum)==currOTUID &
             !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus)))
  currDoseDiff <- unique(currDf$a_doseDifference)
  # Assess whether genus-level complementation is possible.
  currGenus <- unique(currDf$Genus)
  GenusDf <- theoretical_mixtures_df_categorized %>% 
    ungroup() %>% 
    filter(combo==currCombo & passage==currPassage & Genus==currGenus &
             category!="lowAbundance" & a_doseDifference!=0 &
             !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>% 
    mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
    filter(OTUID != currOTUID) %>%
    ungroup() %>%
    select(Phylum, Class, Order, Family, Genus, OTUnum, OTUID, a_doseDifference) %>%
    unique()
  complementationGenusDf <- foreach(x=seq(1:(nrow(GenusDf)*1000)), .combine="rbind") %dopar% {
    nOTUs <- sample(1:nrow(GenusDf), 1)
    tempGenusDf <- GenusDf %>%
      slice_sample(n=nOTUs) %>%
      mutate(rep=x,
             complemented=abs(sum(a_doseDifference) + currDoseDiff) <= 0.09)
    }
  if (TRUE %in% complementationGenusDf$complemented) {
    lowestTaxonomicLevel <- "Genus"
  }
  if (!(TRUE %in% complementationGenusDf$complemented)) {
    # Assess whether family-level complementation is possible.
    currFamily <- unique(currDf$Family)
    FamilyDf <- theoretical_mixtures_df_categorized %>%
      filter(combo==currCombo & passage==currPassage & Family==currFamily &
               category!="lowAbundance" & a_doseDifference!=0 &
               !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>%
      mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
      filter(OTUID != currOTUID) %>%
      ungroup() %>%
      select(Phylum, Class, Order, Family, Genus, OTUnum, OTUID, a_doseDifference) %>%
      unique()
    complementationFamilyDf <- foreach(x=seq(1:(nrow(FamilyDf)*1000)), .combine="rbind") %dopar% {
      nOTUs <- sample(1:nrow(FamilyDf), 1)
      tempFamilyDf <- FamilyDf %>%
        slice_sample(n=nOTUs) %>%
        mutate(rep=x,
               complemented=abs(sum(a_doseDifference) + currDoseDiff) <= 0.09)
    }
    if (TRUE %in% complementationFamilyDf$complemented) {
      lowestTaxonomicLevel <- "Family"
    }
    if (!(TRUE %in% complementationFamilyDf$complemented)) {
      # Assess whether order-level complementation is possible.
      currOrder <- unique(currDf$Order)
      OrderDf <- theoretical_mixtures_df_categorized %>%
        filter(combo==currCombo & passage==currPassage & Order==currOrder &
                 category!="lowAbundance" & a_doseDifference!=0 &
                 !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>%
        mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
        filter(OTUID != currOTUID) %>%
        ungroup() %>%
        select(Phylum, Class, Order, Family, Genus, OTUnum, OTUID, a_doseDifference) %>%
        unique()
      complementationOrderDf <- foreach(x=seq(1:(nrow(OrderDf)*1000)), .combine="rbind") %dopar% {
        nOTUs <- sample(1:nrow(OrderDf), 1)
        tempOrderDf <- OrderDf %>%
          slice_sample(n=nOTUs) %>%
          mutate(rep=x,
                 complemented=abs(sum(a_doseDifference) + currDoseDiff) <= 0.09)
      }
      if (TRUE %in% complementationOrderDf$complemented) {
        lowestTaxonomicLevel <- "Order"
      }
      if (!(TRUE %in% complementationOrderDf$complemented)) {
        # Assess whether class-level complementation is possible.
        currClass <- unique(currDf$Class)
        ClassDf <- theoretical_mixtures_df_categorized %>%
          filter(combo==currCombo & passage==currPassage & Class==currClass &
                   category!="lowAbundance" & a_doseDifference!=0 &
                   !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>%
          mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
          filter(OTUID != currOTUID) %>%
          ungroup() %>%
          select(Phylum, Class, Order, Family, Genus, OTUnum, OTUID, a_doseDifference) %>%
          unique()
        complementationClassDf <- foreach(x=seq(1:(nrow(ClassDf)*1000)), .combine="rbind") %dopar% {
          nOTUs <- sample(1:nrow(ClassDf), 1)
          tempClassDf <- ClassDf %>%
            slice_sample(n=nOTUs) %>%
            mutate(rep=x,
                   complemented=abs(sum(a_doseDifference) + currDoseDiff) <= 0.09)
        }
        if (TRUE %in% complementationClassDf$complemented) {
          lowestTaxonomicLevel <- "Class"
        }
        if (!(TRUE %in% complementationClassDf$complemented)) {
          # Assess whether phylum-level complementation is possible.
          currPhylum <- unique(currDf$Phylum)
          PhylumDf <- theoretical_mixtures_df_categorized %>%
            filter(combo==currCombo & passage==currPassage & Phylum==currPhylum &
                     category!="lowAbundance" & a_doseDifference!=0 &
                     !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>%
            mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
            filter(OTUID != currOTUID) %>%
            ungroup() %>%
            select(Phylum, Class, Order, Family, Genus, OTUnum, OTUID, a_doseDifference) %>%
            unique()
          complementationPhylumDf <- foreach(x=seq(1:(nrow(PhylumDf)*1000)), .combine="rbind") %dopar% {
            nOTUs <- sample(1:nrow(PhylumDf), 1)
            tempPhylumDf <- PhylumDf %>%
              slice_sample(n=nOTUs) %>%
              mutate(rep=x,
                     complemented=abs(sum(a_doseDifference) + currDoseDiff) <= 0.09)
          }
          if (TRUE %in% complementationPhylumDf$complemented) {
            lowestTaxonomicLevel <- "Phylum"
          }
          if (!(TRUE %in% complementationPhylumDf$complemented)) {
            lowestTaxonomicLevel <- "None"
          }
        }
      }
    }
  }
  lowestTaxaDf <- currDf %>%
    ungroup() %>%
    select(passage, combo, Phylum, Class, Order, Family, Genus, OTUnum) %>%
    unique() %>%
    mutate(lowestTaxonomicLevel = lowestTaxonomicLevel)
}

#write.table(DDcomplementationDf, paste0(plotdir, "/DDcomplementationDataframe.txt"), 
#            sep="\t", row.names=FALSE, quote=FALSE)

DDcomplementationDf <- fread(paste0(plotdir, "/DDcomplementationDataframe.txt"))

# Bind relative abundances to dose-dependent ASVs.
DDcomplementationDf <- DDcomplementationDf %>% right_join(lowestTaxaRelAb) %>% 
  ungroup() %>% 
  mutate(lowestTaxonomicLevel = fct_relevel(
    lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus"))) %>% 
  group_by(passage, combo, ratio, lowestTaxonomicLevel) %>% 
  arrange(-avgAbundance) %>% 
  ungroup()

# Count instances of each level of lowest shared taxa.
DDcomplementationLowestTaxaOpposite <- DDcomplementationDf %>%
  select(-c(ratio, avgAbundance)) %>% 
  group_by(passage, combo, lowestTaxonomicLevel) %>% 
  unique() %>% 
  dplyr::summarize(count=n())

testTaxa <- theoretical_mixtures_df_categorized %>% 
  ungroup() %>% 
  select(Phylum, Class, Order, Family, Genus, OTUnum) %>% 
  unique() %>% 
  select(!OTUnum) %>% 
  slice_sample(1)

# Permute taxa of each DD ASV and recalculate complementation.
DDcomplementationDf_DDtaxaPermuted <- foreach(x=DDASVs_comboPassages, .combine="rbind") %dopar% {
  currOTUID <- sub("\\..*","",x)
  currCombo <- sub(":.*","",sub(".*\\.","",x))
  currPassage <- sub(".*:","",x)
  foreach(x=seq(1:1000, .combine="rbind")) %dopar% {
    currDf <- theoretical_mixtures_df_categorized %>% 
      filter(combo==currCombo & passage==currPassage & paste0(Family, "-", OTUnum)==currOTUID &
               !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>% 
      select(-c(Phylum, Class, Order, Genus)) %>% 
      left_join(theoretical_mixtures_df_categorized %>% 
                  ungroup() %>% 
                  select(Phylum, Class, Order, Family, Genus, OTUnum) %>% 
                  unique() %>% 
                  select(!OTUnum) %>% 
                  slice_sample(1)) %>% 
      mutate(permutationReplicate=x)
  }
  currDoseDiff <- unique(currDf$a_doseDifference)
  # Assess whether genus-level complementation is possible.
  currGenus <- unique(currDf$Genus)
  GenusDf <- theoretical_mixtures_df_categorized %>% 
    ungroup() %>% 
    filter(combo==currCombo & passage==currPassage & Genus==currGenus &
             category!="lowAbundance" & a_doseDifference!=0 &
             !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>% 
    mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
    filter(OTUID != currOTUID) %>%
    ungroup() %>%
    select(Phylum, Class, Order, Family, Genus, OTUnum, OTUID, a_doseDifference) %>%
    unique()
  complementationGenusDf <- foreach(x=seq(1:(nrow(GenusDf)*1000)), .combine="rbind") %dopar% {
    nOTUs <- sample(1:nrow(GenusDf), 1)
    tempGenusDf <- GenusDf %>%
      slice_sample(n=nOTUs) %>%
      mutate(rep=x,
             complemented=abs(sum(a_doseDifference) + currDoseDiff) <= 0.09)
  }
  if (TRUE %in% complementationGenusDf$complemented) {
    lowestTaxonomicLevel <- "Genus"
  }
  if (!(TRUE %in% complementationGenusDf$complemented)) {
    # Assess whether family-level complementation is possible.
    currFamily <- unique(currDf$Family)
    FamilyDf <- theoretical_mixtures_df_categorized %>%
      filter(combo==currCombo & passage==currPassage & Family==currFamily &
               category!="lowAbundance" & a_doseDifference!=0 &
               !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>%
      mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
      filter(OTUID != currOTUID) %>%
      ungroup() %>%
      select(Phylum, Class, Order, Family, Genus, OTUnum, OTUID, a_doseDifference) %>%
      unique()
    complementationFamilyDf <- foreach(x=seq(1:(nrow(FamilyDf)*1000)), .combine="rbind") %dopar% {
      nOTUs <- sample(1:nrow(FamilyDf), 1)
      tempFamilyDf <- FamilyDf %>%
        slice_sample(n=nOTUs) %>%
        mutate(rep=x,
               complemented=abs(sum(a_doseDifference) + currDoseDiff) <= 0.09)
    }
    if (TRUE %in% complementationFamilyDf$complemented) {
      lowestTaxonomicLevel <- "Family"
    }
    if (!(TRUE %in% complementationFamilyDf$complemented)) {
      # Assess whether order-level complementation is possible.
      currOrder <- unique(currDf$Order)
      OrderDf <- theoretical_mixtures_df_categorized %>%
        filter(combo==currCombo & passage==currPassage & Order==currOrder &
                 category!="lowAbundance" & a_doseDifference!=0 &
                 !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>%
        mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
        filter(OTUID != currOTUID) %>%
        ungroup() %>%
        select(Phylum, Class, Order, Family, Genus, OTUnum, OTUID, a_doseDifference) %>%
        unique()
      complementationOrderDf <- foreach(x=seq(1:(nrow(OrderDf)*1000)), .combine="rbind") %dopar% {
        nOTUs <- sample(1:nrow(OrderDf), 1)
        tempOrderDf <- OrderDf %>%
          slice_sample(n=nOTUs) %>%
          mutate(rep=x,
                 complemented=abs(sum(a_doseDifference) + currDoseDiff) <= 0.09)
      }
      if (TRUE %in% complementationOrderDf$complemented) {
        lowestTaxonomicLevel <- "Order"
      }
      if (!(TRUE %in% complementationOrderDf$complemented)) {
        # Assess whether class-level complementation is possible.
        currClass <- unique(currDf$Class)
        ClassDf <- theoretical_mixtures_df_categorized %>%
          filter(combo==currCombo & passage==currPassage & Class==currClass &
                   category!="lowAbundance" & a_doseDifference!=0 &
                   !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>%
          mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
          filter(OTUID != currOTUID) %>%
          ungroup() %>%
          select(Phylum, Class, Order, Family, Genus, OTUnum, OTUID, a_doseDifference) %>%
          unique()
        complementationClassDf <- foreach(x=seq(1:(nrow(ClassDf)*1000)), .combine="rbind") %dopar% {
          nOTUs <- sample(1:nrow(ClassDf), 1)
          tempClassDf <- ClassDf %>%
            slice_sample(n=nOTUs) %>%
            mutate(rep=x,
                   complemented=abs(sum(a_doseDifference) + currDoseDiff) <= 0.09)
        }
        if (TRUE %in% complementationClassDf$complemented) {
          lowestTaxonomicLevel <- "Class"
        }
        if (!(TRUE %in% complementationClassDf$complemented)) {
          # Assess whether phylum-level complementation is possible.
          currPhylum <- unique(currDf$Phylum)
          PhylumDf <- theoretical_mixtures_df_categorized %>%
            filter(combo==currCombo & passage==currPassage & Phylum==currPhylum &
                     category!="lowAbundance" & a_doseDifference!=0 &
                     !(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>%
            mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
            filter(OTUID != currOTUID) %>%
            ungroup() %>%
            select(Phylum, Class, Order, Family, Genus, OTUnum, OTUID, a_doseDifference) %>%
            unique()
          complementationPhylumDf <- foreach(x=seq(1:(nrow(PhylumDf)*1000)), .combine="rbind") %dopar% {
            nOTUs <- sample(1:nrow(PhylumDf), 1)
            tempPhylumDf <- PhylumDf %>%
              slice_sample(n=nOTUs) %>%
              mutate(rep=x,
                     complemented=abs(sum(a_doseDifference) + currDoseDiff) <= 0.09)
          }
          if (TRUE %in% complementationPhylumDf$complemented) {
            lowestTaxonomicLevel <- "Phylum"
          }
          if (!(TRUE %in% complementationPhylumDf$complemented)) {
            lowestTaxonomicLevel <- "None"
          }
        }
      }
    }
  }
  lowestTaxaDf <- currDf %>%
    ungroup() %>%
    select(passage, combo, Phylum, Class, Order, Family, Genus, OTUnum) %>%
    unique() %>%
    mutate(lowestTaxonomicLevel = lowestTaxonomicLevel)
}

# Quantify complementation for overcolonizing and undercolonizing ASVs. --------


# Quantify dose-dependence among taxa compared to category permutation. --------

# Permute random combinations of dose-dependent ASVs and calculate lowest shared taxonomic levels.
permuteDDtaxa <- function(comboPassage) {
  currCombo <- sub("_.*","",comboPassage)
  currPassage <- sub(".*_","",comboPassage)
  set.seed(0)
  subset_df <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
    filter(passage==currPassage & combo==currCombo & category!="lowAbundance") %>% 
    filter(!(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>% 
    ungroup() %>% 
    dplyr::select(passage, combo, Phylum, Class, Order, Family, Genus, OTU, a_doseDifference) %>% 
    unique()
  numASVs <- n_distinct((theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
                           filter(passage==currPassage & combo==currCombo & category=="dose-dependent"))$OTU)
  DD_df <- foreach(x=seq(1:1000), .combine="rbind") %do% {
      subset_df %>% 
      slice_sample(n=numASVs) %>% 
      mutate(permutationReplicate = x, direction = a_doseDifference > 0)
  }
  DD_df <- join_all(list(
    DD_df %>% 
      group_by(permutationReplicate, Phylum) %>%
      dplyr::summarize(passage, combo, Class, Order, Family, Genus, OTU, oppositeDirectionsPhylum = TRUE %in% direction & FALSE %in% direction) %>% 
      unique(),
    DD_df %>% 
      group_by(permutationReplicate, Phylum, Class) %>%
      dplyr::summarize(passage, combo, Order, Family, Genus, OTU, oppositeDirectionsClass = TRUE %in% direction & FALSE %in% direction) %>% 
      unique(),
    DD_df %>%  
      group_by(permutationReplicate, Phylum, Class, Order) %>%
      dplyr::summarize(passage, combo, Family, Genus, OTU, oppositeDirectionsOrder = TRUE %in% direction & FALSE %in% direction) %>% 
      unique(),
    DD_df %>% 
      group_by(permutationReplicate, Phylum, Class, Order, Family) %>%
      dplyr::summarize(passage, combo, Genus, OTU, oppositeDirectionsFamily = TRUE %in% direction & FALSE %in% direction) %>% 
      unique(),
    DD_df %>% 
      group_by(permutationReplicate, Phylum, Class, Order, Family, Genus) %>% 
      dplyr::summarize(passage, combo, OTU, oppositeDirectionsGenus = TRUE %in% direction & FALSE %in% direction & !is.na(Genus)) %>% 
      unique()
    ), type='left') %>% 
    mutate(oppositeDirectionsPhylum = oppositeDirectionsPhylum & !oppositeDirectionsClass & 
             !oppositeDirectionsOrder & !oppositeDirectionsFamily & !oppositeDirectionsGenus,
           oppositeDirectionsClass = oppositeDirectionsClass & !oppositeDirectionsOrder & 
             !oppositeDirectionsFamily & !oppositeDirectionsGenus,
           oppositeDirectionsOrder = oppositeDirectionsOrder & !oppositeDirectionsFamily &
             !oppositeDirectionsGenus,
           oppositeDirectionsFamily = oppositeDirectionsFamily & !oppositeDirectionsGenus) %>% 
    group_by(passage, combo, permutationReplicate, Phylum, Class, Order, Family, Genus, OTU) %>% 
    dplyr::summarize(oppositeDirectionsPhylum, oppositeDirectionsClass, oppositeDirectionsOrder, 
                     oppositeDirectionsFamily, oppositeDirectionsGenus,
                     lowestTaxonomicLevel = case_when(
                       TRUE %in% oppositeDirectionsPhylum ~ "Phylum",
                       TRUE %in% oppositeDirectionsClass ~ "Class",
                       TRUE %in% oppositeDirectionsOrder ~ "Order",
                       TRUE %in% oppositeDirectionsFamily ~ "Family",
                       TRUE %in% oppositeDirectionsGenus ~ "Genus"),
                     lowestTaxonomicLevel = replace(lowestTaxonomicLevel, is.na(lowestTaxonomicLevel), "None"))
  return(DD_df)
}
# Generate all permuted dataframes in a list.
permutedTaxa <- list.rbind(lapply(comboPassages, permuteDDtaxa))

# Get relative abundances of all ASVs at key mixture ratios.
lowestTaxaRelAbAll <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
  filter(category!="lowAbundance" & 
           (ratio %in% c("1:0","0:1") | (ratio=="1:1" & mixtureType=="actual")) & !is.na(rel_abundance_old)) %>% 
  group_by(passage, combo, ratio, Phylum, Class, Order, Family, Genus, OTUnum) %>% 
  dplyr::summarize(avgAbundance=ifelse(ratio=="1:1", sum(rel_abundance_old/3) , sum(rel_abundance_old/2))) %>% 
  unique()
# Bind relative abundances to permuted ASVs.
# Sum relative abundances among shared taxonomic levels.
permutedTaxaRelAb <- permutedTaxa %>% right_join(lowestTaxaRelAbAll) %>%  
  ungroup() %>% 
  group_by(passage, combo, ratio, permutationReplicate, lowestTaxonomicLevel) %>% 
  dplyr::summarize(sumAbundance = sum(avgAbundance)) %>%
  ungroup() %>% 
  complete(passage, combo, ratio, permutationReplicate, lowestTaxonomicLevel, fill=list(sumAbundance=0)) %>% 
  dplyr::rename(permutedAbundance=sumAbundance) %>% 
  left_join(taxaDoseDependence_ASV %>% 
              group_by(passage, combo, ratio, lowestTaxonomicLevel) %>% 
              dplyr::summarize(realAbundance=sum(avgAbundance))) %>% 
  ungroup() %>% 
  mutate(realAbundance=replace(realAbundance, is.na(realAbundance), 0),
         lowestTaxonomicLevel = fct_relevel(
           lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))

# Count permuted instances of each level of lowest shared taxa and bind to actual instances.
permutedLowestTaxaOpposite <- permutedTaxa %>%
  group_by(passage, combo, permutationReplicate, lowestTaxonomicLevel) %>% 
  unique() %>% 
  dplyr::summarize(count=n()) %>%
  ungroup() %>% 
  complete(passage, combo, permutationReplicate, lowestTaxonomicLevel, fill=list(count=0)) %>% 
  dplyr::rename(permutedCount=count) %>% 
  left_join(lowestTaxaOpposite %>% dplyr::rename(realCount=count)) %>% 
  mutate(realCount=replace(realCount, is.na(realCount), 0)) %>% 
  ungroup() %>% 
  mutate(lowestTaxonomicLevel = fct_relevel(
    lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))

# Quantify dose-dependence among taxa compared to category permutation, all categories. --------

# Permute random combinations of dose-dependent ASVs and calculate lowest shared taxonomic levels, including
# ASVs from all categories.
permuteDDtaxa_allCategories <- function(comboPassage) {
  currCombo <- sub("_.*","",comboPassage)
  currPassage <- sub(".*_","",comboPassage)
  set.seed(0)
  subset_df <- theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
    filter(passage==currPassage & combo==currCombo & category!="lowAbundance") %>% 
    filter(!(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>% 
    ungroup() %>% 
    dplyr::select(passage, combo, Phylum, Class, Order, Family, Genus, OTUnum, a_doseDifference) %>% 
    unique()
  numASVs <- n_distinct((theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>% 
                           filter(passage==currPassage & combo==currCombo & category=="dose-dependent"))$OTU)
  DD_df <- foreach(x=seq(1:1000), .combine="rbind") %do% {
    subset_df %>% 
      dplyr::mutate(permutationReplicate = x, 
                    category = "other",
                    category = replace(category, sample(n(), numASVs), "dose-dependent"),
                    direction = case_when(category=="dose-dependent" & a_doseDifference > 0 ~ "positive",
                              category=="dose-dependent" & a_doseDifference <= 0 ~ "negative",
                              category!="dose-dependent" ~ "NA")) 
  }
  DD_df <- join_all(list(
    DD_df %>% 
      group_by(permutationReplicate, Phylum) %>%
      dplyr::summarize(passage, combo, Class, Order, Family, Genus, OTUnum, oppositeDirectionsPhylum = "positive" %in% direction & 
                         "negative" %in% direction & abs(sum(unique(a_doseDifference))) <= 0.31) %>% 
      unique(),
    DD_df %>% 
      group_by(permutationReplicate, Phylum, Class) %>%
      dplyr::summarize(passage, combo, Order, Family, Genus, OTUnum, oppositeDirectionsClass = "positive" %in% direction & 
                         "negative" %in% direction & abs(sum(unique(a_doseDifference))) <= 0.31) %>% 
      unique(),
    DD_df %>%  
      group_by(permutationReplicate, Phylum, Class, Order) %>%
      dplyr::summarize(passage, combo, Family, Genus, OTUnum, oppositeDirectionsOrder = "positive" %in% direction & 
                         "negative" %in% direction & abs(sum(unique(a_doseDifference))) <= 0.31) %>% 
      unique(),
    DD_df %>% 
      group_by(permutationReplicate, Phylum, Class, Order, Family) %>%
      dplyr::summarize(passage, combo, Genus, OTUnum, oppositeDirectionsFamily = "positive" %in% direction & 
                         "negative" %in% direction & abs(sum(unique(a_doseDifference))) <= 0.31) %>% 
      unique(),
    DD_df %>% 
      group_by(permutationReplicate, Phylum, Class, Order, Family, Genus) %>% 
      dplyr::summarize(passage, combo, OTUnum, oppositeDirectionsGenus = "positive" %in% direction & 
                         "negative" %in% direction & abs(sum(unique(a_doseDifference))) <= 0.31 & !is.na(Genus)) %>% 
      unique()
  ), type='left') %>% 
    mutate(oppositeDirectionsPhylum = oppositeDirectionsPhylum & !oppositeDirectionsClass & 
             !oppositeDirectionsOrder & !oppositeDirectionsFamily & !oppositeDirectionsGenus,
           oppositeDirectionsClass = oppositeDirectionsClass & !oppositeDirectionsOrder & 
             !oppositeDirectionsFamily & !oppositeDirectionsGenus,
           oppositeDirectionsOrder = oppositeDirectionsOrder & !oppositeDirectionsFamily &
             !oppositeDirectionsGenus,
           oppositeDirectionsFamily = oppositeDirectionsFamily & !oppositeDirectionsGenus) %>% 
    group_by(passage, combo, permutationReplicate, Phylum, Class, Order, Family, Genus, OTUnum) %>% 
    dplyr::summarize(oppositeDirectionsPhylum, oppositeDirectionsClass, oppositeDirectionsOrder, 
                     oppositeDirectionsFamily, oppositeDirectionsGenus,
                     lowestTaxonomicLevel = case_when(
                       TRUE %in% oppositeDirectionsPhylum ~ "Phylum",
                       TRUE %in% oppositeDirectionsClass ~ "Class",
                       TRUE %in% oppositeDirectionsOrder ~ "Order",
                       TRUE %in% oppositeDirectionsFamily ~ "Family",
                       TRUE %in% oppositeDirectionsGenus ~ "Genus"),
                     lowestTaxonomicLevel = replace(lowestTaxonomicLevel, is.na(lowestTaxonomicLevel), "None"))
  return(DD_df)
}
# Generate all permuted dataframes in a list.
permutedTaxa_allCategories <- list.rbind(lapply(comboPassages, permuteDDtaxa_allCategories))

# Bind relative abundances to permuted ASVs.
# Sum relative abundances among shared taxonomic levels.
permutedTaxaRelAb_allCategories <- permutedTaxa_allCategories %>% right_join(lowestTaxaRelAbAll) %>%  
  ungroup() %>% 
  group_by(passage, combo, ratio, permutationReplicate, lowestTaxonomicLevel) %>% 
  dplyr::summarize(sumAbundance = sum(avgAbundance)) %>%
  ungroup() %>% 
  complete(passage, combo, ratio, permutationReplicate, lowestTaxonomicLevel, fill=list(sumAbundance=0)) %>% 
  dplyr::rename(permutedAbundance=sumAbundance) %>% 
  left_join(taxaDoseDependence_ASV %>% 
              group_by(passage, combo, ratio, lowestTaxonomicLevel) %>% 
              dplyr::summarize(realAbundance=sum(avgAbundance))) %>% 
  ungroup() %>% 
  mutate(realAbundance=replace(realAbundance, is.na(realAbundance), 0),
         lowestTaxonomicLevel = fct_relevel(
           lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))

# Count permuted instances of each level of lowest shared taxa and bind to actual instances.
permutedLowestTaxaOpposite_allCategories <- permutedTaxa_allCategories %>%
  group_by(passage, combo, permutationReplicate, lowestTaxonomicLevel) %>% 
  unique() %>% 
  dplyr::summarize(count=n()) %>%
  ungroup() %>% 
  complete(passage, combo, permutationReplicate, lowestTaxonomicLevel, fill=list(count=0)) %>% 
  dplyr::rename(permutedCount=count) %>% 
  left_join(lowestTaxaOpposite %>% dplyr::rename(realCount=count)) %>% 
  mutate(realCount=replace(realCount, is.na(realCount), 0)) %>% 
  ungroup() %>% 
  mutate(lowestTaxonomicLevel = fct_relevel(
    lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))

# Quantify dose-dependence among taxa compared to direction permutation. --------

# Permute direction of dose-dependent ASVs and calculate lowest shared taxonomic levels.
set.seed(0)
permutedDirectionsTaxa <- foreach(x=seq(1:1000), .combine="rbind") %do% {
  theoretical_mixtures_HighAbundanceOneParent_summaryFiltered %>%
    filter(category=="dose-dependent") %>% 
    filter(!(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>%
    ungroup() %>% 
    select(passage, combo, Phylum, Class, Order, Family, Genus, OTU) %>% 
    unique() %>% 
    dplyr::mutate(permutationReplicate = x, direction = sample(c(TRUE, FALSE), n(), replace=TRUE))
}
# Add marker for presence of permuted opposite-direction ASVs at each taxonomic level.
permutedDirectionsTaxa <- join_all(list(
  permutedDirectionsTaxa %>% 
    group_by(passage, combo, permutationReplicate, Phylum) %>%
    dplyr::summarize(Class, Order, Family, Genus, OTU, oppositeDirectionsPhylum = TRUE %in% direction & FALSE %in% direction) %>% 
    unique(),
  permutedDirectionsTaxa %>% 
    group_by(passage, combo, permutationReplicate, Phylum, Class) %>%
    dplyr::summarize(Order, Family, Genus, OTU, oppositeDirectionsClass = TRUE %in% direction & FALSE %in% direction) %>% 
    unique(),
  permutedDirectionsTaxa %>%  
    group_by(passage, combo, permutationReplicate, Phylum, Class, Order) %>%
    dplyr::summarize(Family, Genus, OTU, oppositeDirectionsOrder = TRUE %in% direction & FALSE %in% direction) %>% 
    unique(),
  permutedDirectionsTaxa %>% 
    group_by(passage, combo, permutationReplicate, Phylum, Class, Order, Family) %>%
    dplyr::summarize(Genus, OTU, oppositeDirectionsFamily = TRUE %in% direction & FALSE %in% direction) %>% 
    unique(),
  permutedDirectionsTaxa %>% 
    group_by(passage, combo, permutationReplicate, Phylum, Class, Order, Family, Genus) %>% 
    dplyr::summarize(OTU, oppositeDirectionsGenus = TRUE %in% direction & FALSE %in% direction & !is.na(Genus)) %>% 
    unique()
  ), type='left') %>% 
  # Set lowest rank at which opposite direction dose dependence occurs
  mutate(oppositeDirectionsPhylum = oppositeDirectionsPhylum & !oppositeDirectionsClass & 
           !oppositeDirectionsOrder & !oppositeDirectionsFamily & !oppositeDirectionsGenus,
         oppositeDirectionsClass = oppositeDirectionsClass & !oppositeDirectionsOrder & 
           !oppositeDirectionsFamily & !oppositeDirectionsGenus,
         oppositeDirectionsOrder = oppositeDirectionsOrder & !oppositeDirectionsFamily &
           !oppositeDirectionsGenus,
         oppositeDirectionsFamily = oppositeDirectionsFamily & !oppositeDirectionsGenus) %>% 
  group_by(passage, combo, permutationReplicate, Phylum, Class, Order, Family, Genus, OTU) %>% 
  dplyr::summarize(oppositeDirectionsPhylum, oppositeDirectionsClass, oppositeDirectionsOrder, 
                   oppositeDirectionsFamily, oppositeDirectionsGenus,
                   lowestTaxonomicLevel = case_when(
                     TRUE %in% oppositeDirectionsPhylum ~ "Phylum",
                     TRUE %in% oppositeDirectionsClass ~ "Class",
                     TRUE %in% oppositeDirectionsOrder ~ "Order",
                     TRUE %in% oppositeDirectionsFamily ~ "Family",
                     TRUE %in% oppositeDirectionsGenus ~ "Genus"),
                   lowestTaxonomicLevel = replace(lowestTaxonomicLevel, is.na(lowestTaxonomicLevel), "None"))
# Bind relative abundances to dose-dependent ASVs with permuted directions.
permutedDirectionsTaxaRelAb <- permutedDirectionsTaxa %>% right_join(lowestTaxaRelAb) %>% 
  ungroup() %>% 
  group_by(passage, ratio, permutationReplicate, Phylum, lowestTaxonomicLevel) %>% 
  dplyr::summarize(sumAbundance = sum(avgAbundance)) %>% 
  ungroup() %>% 
  complete(passage, ratio, permutationReplicate, Phylum, lowestTaxonomicLevel, fill=list(sumAbundance=0)) %>% 
  dplyr::rename(permutedAbundance=sumAbundance) %>% 
  left_join(taxaDoseDependence_ASV %>% 
              group_by(passage, ratio, Phylum, lowestTaxonomicLevel) %>% 
              dplyr::summarize(realAbundance=sum(avgAbundance))) %>% 
  ungroup() %>% 
  mutate(realAbundance=replace(realAbundance, is.na(realAbundance), 0),
         lowestTaxonomicLevel = fct_relevel(
           lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))

# Count instances of each level of lowest shared taxa with permuted directions and bind to actual instances.
permutedDirectionsLowestTaxaOpposite <- permutedDirectionsTaxa %>%
  group_by(passage, permutationReplicate, Phylum, lowestTaxonomicLevel) %>% 
  unique() %>% 
  dplyr::summarize(count=n()) %>%
  ungroup() %>% 
  complete(passage, permutationReplicate, Phylum, lowestTaxonomicLevel, fill=list(count=0)) %>% 
  dplyr::rename(permutedCount=count) %>% 
  left_join(taxaDoseDependence_ASV %>%
              select(-c(ratio, avgAbundance)) %>% 
              group_by(passage, Phylum, lowestTaxonomicLevel) %>% 
              unique() %>% 
              dplyr::summarize(realCount=n())) %>% 
  mutate(realCount=replace(realCount, is.na(realCount), 0)) %>% 
  ungroup() %>% 
  mutate(lowestTaxonomicLevel = fct_relevel(
    lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))

# Investigate direction permutation for specific taxa below phyla. --------

# Bacteroidota phylum with genus grouping.
permutedDirectionsTaxaRelAb_BacteroidotaGenus <- permutedDirectionsTaxa %>% right_join(lowestTaxaRelAb) %>% 
  ungroup() %>% 
  filter(Phylum=="Bacteroidota") %>% 
  group_by(passage, ratio, permutationReplicate, Genus, lowestTaxonomicLevel) %>% 
  dplyr::summarize(permutedAbundance = sum(avgAbundance)) %>% 
  ungroup() %>% 
  complete(passage, ratio, permutationReplicate, Genus, 
           lowestTaxonomicLevel=c("None","Phylum","Class","Order","Family","Genus"), 
           fill=list(permutedAbundance=0)) %>% 
  left_join(taxaDoseDependence_ASV %>% 
              group_by(passage, ratio, Genus, lowestTaxonomicLevel) %>% 
              dplyr::summarize(realAbundance=sum(avgAbundance))) %>% 
  ungroup() %>% 
  mutate(realAbundance=replace(realAbundance, is.na(realAbundance), 0),
         lowestTaxonomicLevel = fct_relevel(
           lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))

# Bacteroidota phylum with genus grouping.
permutedDirectionsLowestTaxaOpposite_BacteroidotaGenus <- permutedDirectionsTaxa %>% 
  filter(Phylum=="Bacteroidota") %>% 
  group_by(passage, permutationReplicate, Genus, lowestTaxonomicLevel) %>% 
  unique() %>% 
  dplyr::summarize(permutedCount=n()) %>% 
  ungroup() %>% 
  complete(passage, permutationReplicate, Genus, 
           lowestTaxonomicLevel=c("None","Phylum","Class","Order","Family","Genus"),
           fill=list(permutedCount=0)) %>% 
  left_join(taxaDoseDependence_ASV %>%
              select(-c(ratio, avgAbundance)) %>% 
              group_by(passage, Genus, lowestTaxonomicLevel) %>% 
              unique() %>% 
              dplyr::summarize(realCount=n())) %>% 
  mutate(realCount=replace(realCount, is.na(realCount), 0)) %>% 
  ungroup() %>% 
  mutate(lowestTaxonomicLevel = fct_relevel(
    lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))

# Firmicutes phylum with class grouping.
permutedDirectionsTaxaRelAb_FirmicutesClass <- permutedDirectionsTaxa %>% right_join(lowestTaxaRelAb) %>% 
  ungroup() %>% 
  filter(Phylum=="Firmicutes") %>% 
  group_by(passage, ratio, permutationReplicate, Class, lowestTaxonomicLevel) %>% 
  dplyr::summarize(permutedAbundance = sum(avgAbundance)) %>% 
  ungroup() %>% 
  complete(passage, ratio, permutationReplicate, Class, 
           lowestTaxonomicLevel=c("None","Phylum","Class","Order","Family","Genus"), 
           fill=list(permutedAbundance=0)) %>% 
  left_join(taxaDoseDependence_ASV %>% 
              group_by(passage, ratio, Class, lowestTaxonomicLevel) %>% 
              dplyr::summarize(realAbundance=sum(avgAbundance))) %>% 
  ungroup() %>% 
  mutate(realAbundance=replace(realAbundance, is.na(realAbundance), 0),
         lowestTaxonomicLevel = fct_relevel(
           lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))

# Firmicutes phylum with class grouping.
permutedDirectionsLowestTaxaOpposite_FirmicutesClass <- permutedDirectionsTaxa %>% 
  filter(Phylum=="Firmicutes") %>% 
  group_by(passage, permutationReplicate, Class, lowestTaxonomicLevel) %>% 
  unique() %>% 
  dplyr::summarize(permutedCount=n()) %>% 
  ungroup() %>% 
  complete(passage, permutationReplicate, Class, 
           lowestTaxonomicLevel=c("None","Phylum","Class","Order","Family","Genus"),
           fill=list(permutedCount=0)) %>% 
  left_join(taxaDoseDependence_ASV %>%
              select(-c(ratio, avgAbundance)) %>% 
              group_by(passage, Class, lowestTaxonomicLevel) %>% 
              unique() %>% 
              dplyr::summarize(realCount=n())) %>% 
  mutate(realCount=replace(realCount, is.na(realCount), 0)) %>% 
  ungroup() %>% 
  mutate(lowestTaxonomicLevel = fct_relevel(
    lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))

# Firmicutes phylum with order grouping.
permutedDirectionsTaxaRelAb_FirmicutesOrder <- permutedDirectionsTaxa %>% right_join(lowestTaxaRelAb) %>% 
  ungroup() %>% 
  filter(Phylum=="Firmicutes") %>% 
  group_by(passage, ratio, permutationReplicate, Order, lowestTaxonomicLevel) %>% 
  dplyr::summarize(permutedAbundance = sum(avgAbundance)) %>% 
  ungroup() %>% 
  complete(passage, ratio, permutationReplicate, Order, 
           lowestTaxonomicLevel=c("None","Phylum","Class","Order","Family","Genus"), 
           fill=list(permutedAbundance=0)) %>% 
  left_join(taxaDoseDependence_ASV %>% 
              group_by(passage, ratio, Order, lowestTaxonomicLevel) %>% 
              dplyr::summarize(realAbundance=sum(avgAbundance))) %>% 
  ungroup() %>% 
  mutate(realAbundance=replace(realAbundance, is.na(realAbundance), 0),
         lowestTaxonomicLevel = fct_relevel(
           lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))

# Firmicutes phylum with order grouping.
permutedDirectionsLowestTaxaOpposite_FirmicutesOrder <- permutedDirectionsTaxa %>% 
  filter(Phylum=="Firmicutes") %>% 
  group_by(passage, permutationReplicate, Order, lowestTaxonomicLevel) %>% 
  unique() %>% 
  dplyr::summarize(permutedCount=n()) %>% 
  ungroup() %>% 
  complete(passage, permutationReplicate, Order, 
           lowestTaxonomicLevel=c("None","Phylum","Class","Order","Family","Genus"),
           fill=list(permutedCount=0)) %>% 
  left_join(taxaDoseDependence_ASV %>%
              select(-c(ratio, avgAbundance)) %>% 
              group_by(passage, Order, lowestTaxonomicLevel) %>% 
              unique() %>% 
              dplyr::summarize(realCount=n())) %>% 
  mutate(realCount=replace(realCount, is.na(realCount), 0)) %>% 
  ungroup() %>% 
  mutate(lowestTaxonomicLevel = fct_relevel(
    lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")))
  

# Investigate complementation of Strep and Entero in XFA-XFB. --------------

df_Strep7 <- theoretical_mixtures_df_categorized %>% 
  filter(passage==5 & combo=="XFA-XFB" & Family=="Streptococcaceae" & OTUnum %in% c(7, 17))
sum(unique(df_Strep7$a_doseDifference))
# 0.01926611

df_Entero3 <- theoretical_mixtures_df_categorized %>% 
  filter(passage==5 & combo=="XFA-XFB" & 
           ((Family=="Streptococcaceae" & OTUnum==17) | (Family=="Enterococcaceae" & OTUnum==3)))
sum(unique(df_Entero3$a_doseDifference))
# 0.1153041

df_Entero23 <- theoretical_mixtures_df_categorized %>% 
  filter(passage==5 & combo=="XFA-XFB" & 
           ((Family=="Streptococcaceae" & OTUnum==17) | (Family=="Enterococcaceae" & OTUnum %in% c(2,3))))
sum(unique(df_Entero23$a_doseDifference))
# -0.02079666

df_EnteroAll <- theoretical_mixtures_df_categorized %>% 
  filter(passage==5 & combo=="XFA-XFB" & 
           ((Family=="Streptococcaceae" & OTUnum==17) | (Family=="Enterococcaceae" & category!="lowAbundance")))
sum(unique(df_EnteroAll$a_doseDifference))
# -0.0437508

# Plot dose-dependent taxa with opposite partners. ------------------------

# Plot number of each level of lowest shared taxa.
p_lowestTaxaOpposite <- lowestTaxaOpposite %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_bar(aes(x=combo, y=count, fill=lowestTaxonomicLevel), stat="identity", color="black") +
  #scale_fill_brewer(name="Lowest shared\ntaxonomic level", palette="BuPu") +
  scale_fill_manual(name="Lowest shared\ntaxonomic level", values=c("#EDF8FB","#810F7C")) +
  ylim(c(0,10)) +
  xlab("Mixture") +
  ylab("Count") +
  DEFAULTS.THEME_PRES
p_lowestTaxaOpposite
save_plot(paste0(plotdir, "/lowestTaxaOpposite.png"), p_lowestTaxaOpposite, nrow=2, ncol=1.5)

# Plot relative abundances of each level of lowest shared taxa.
p_lowestTaxaOppositeRelAb <- taxaDoseDependence_ASV %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_bar(aes(x=ratio, y=avgAbundance, fill=lowestTaxonomicLevel), stat="identity", color="black") +
  #scale_fill_brewer(name="Lowest shared\ntaxonomic level", palette="BuPu") +
  scale_fill_manual(name="Lowest shared\ntaxonomic level", values=c("#EDF8FB","#810F7C")) +
  xlab("Ratio") +
  ylab("Relative abundance") +
  facet_grid(~combo) +
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeRelAb
save_plot(paste0(plotdir, "/lowestTaxaOppositeRelAb.png"), p_lowestTaxaOppositeRelAb, nrow=2, ncol=1.5)

# Plot number of each level of lowest shared taxa, all categories.
p_lowestTaxaOpposite_allCategories <- lowestTaxaOpposite_allCategories %>% 
  filter(passage==5) %>% 
  dplyr::mutate(combo=case_when(
    combo=="XBA-XBB" ~ "XBA\nXBB",
    combo=="XCA-XCB" ~ "XCA\nXCB",
    combo=="XDA-XDB" ~ "XDA\nXDB",
    combo=="XFA-XFB" ~ "XFA\nXFB",
    combo=="XBA-XCA" ~ "XBA\nXCA",
    combo=="XCA-XDA" ~ "XCA\nXDA",
    combo=="XDA-XFA" ~ "XDA\nXFA",
    combo=="XFA-XBA" ~ "XFA\nXBA"
  ),
  combo = fct_relevel(combo, c("XBA\nXBB","XCA\nXCB","XDA\nXDB","XFA\nXFB","XBA\nXCA","XCA\nXDA","XDA\nXFA","XFA\nXBA"))) %>% 
  ungroup() %>% 
  complete(combo, lowestTaxonomicLevel=c("None","Phylum","Class","Order","Family","Genus"), fill=list(count=0)) %>%
  mutate(lowestTaxonomicLevel = fct_relevel(lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus"))) %>% 
  ggplot() +
  geom_bar(aes(x=combo, y=count, fill=lowestTaxonomicLevel), stat="identity", color="black", size=0.25) +
  scale_fill_brewer(name="Lowest shared\ntaxonomic level", palette="BuPu") +
  scale_y_continuous("Count", breaks=c(0,2,4,6,8,10), limits=c(0,10)) +
  xlab("Mixture") +
  ylab("Count") +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_lowestTaxaOpposite_allCategories
#save_plot(paste0(plotdir, "/lowestTaxaOpposite_allCategories.png"), 
#          p_lowestTaxaOpposite_allCategories, nrow=2, ncol=1.5)
save_plot(paste0(plotdir, "/lowestTaxaOpposite_allCategories.pdf"), 
          p_lowestTaxaOpposite_allCategories, base_width=3, base_height=1.75)

lowestTaxonomicLevelLegend <- get_legend(p_lowestTaxaOpposite_allCategories)
save_plot(paste0(plotdir, "/lowestTaxonomicLevelLegend.pdf"), lowestTaxonomicLevelLegend, base_width=0.5, base_height=0.75)

# Plot relative abundances of each level of lowest shared taxa, all categories.
p_lowestTaxaOppositeRelAb_allCategories <- taxaDoseDependence_ASV_allCategories %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_bar(aes(x=ratio, y=avgAbundance, fill=lowestTaxonomicLevel), stat="identity", color="black") +
  scale_fill_brewer(name="Lowest shared\ntaxonomic level", palette="BuPu") +
  xlab("Ratio") +
  ylab("Relative abundance") +
  facet_grid(~combo) +
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeRelAb_allCategories
save_plot(paste0(plotdir, "/lowestTaxaOppositeRelAb_allCategories.png"), 
          p_lowestTaxaOppositeRelAb_allCategories, nrow=2, ncol=1.5)

# Plot number of each level of lowest shared taxa, all subsets of all categories.
p_lowestTaxaOpposite_allSubsetsAllCategories <- DDcomplementationLowestTaxaOpposite %>%
  ungroup() %>% 
  filter(passage==5) %>% 
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
  ungroup() %>% 
  complete(combo, lowestTaxonomicLevel=c("None","Phylum","Class","Order","Family","Genus"), fill=list(count=0)) %>%
  mutate(lowestTaxonomicLevel = fct_relevel(lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus"))) %>% 
  ggplot() +
  geom_bar(aes(x=combo, y=count, fill=lowestTaxonomicLevel), stat="identity", color="black", size=0.25) +
  scale_fill_brewer(name="Lowest shared\ntaxonomic level", palette="BuPu") +
  # scale_fill_manual(values=c("None"="#808080",
  #                            "Phylum"="#FFC2FF",
  #                            "Class"="#D15CFF",
  #                            "Order"="#8A59EE",
  #                            "Family"="#4F0CD4",
  #                            "Genus"="#04214E")) +
  scale_y_continuous("# DD ASVs", breaks=c(0,2,4,6,8,10), limits=c(0,10)) +
  xlab("Mixture") +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_lowestTaxaOpposite_allSubsetsAllCategories
#save_plot(paste0(plotdir, "/lowestTaxaOpposite_allCategories.png"), 
#          p_lowestTaxaOpposite_allCategories, nrow=2, ncol=1.5)
save_plot(paste0(plotdir, "/lowestTaxaOpposite_allSubsetsAllCategories.pdf"), 
          p_lowestTaxaOpposite_allSubsetsAllCategories, base_width=2.7, base_height=1.6)

lowestTaxonomicLevelLegend <- get_legend(p_lowestTaxaOpposite_allCategories)
save_plot(paste0(plotdir, "/lowestTaxonomicLevelLegend.pdf"), lowestTaxonomicLevelLegend, base_width=0.5, base_height=0.75)

# Plot relative abundance of each level of lowest shared taxa, all subsets of all categories.
p_lowestTaxaOppositeRelAb_allSubsetsAllCategories <- DDcomplementationDf %>%
  ungroup() %>% 
  filter(passage==5) %>% 
  dplyr::mutate(combo=case_when(
    combo=="XBA-XBB" ~ "XBA\nXBB",
    combo=="XCA-XCB" ~ "XCA\nXCB",
    combo=="XDA-XDB" ~ "XDA\nXDB",
    combo=="XFA-XFB" ~ "XFA\nXFB",
    combo=="XBA-XCA" ~ "XBA\nXCA",
    combo=="XCA-XDA" ~ "XCA\nXDA",
    combo=="XDA-XFA" ~ "XDA\nXFA",
    combo=="XFA-XBA" ~ "XFA\nXBA"
  ),
  combo = fct_relevel(combo, c("XBA\nXBB","XCA\nXCB","XDA\nXDB","XFA\nXFB","XBA\nXCA","XCA\nXDA","XDA\nXFA","XFA\nXBA"))) %>% 
  ungroup() %>% 
  complete(combo, ratio=c("1:0","1:1","0:1"), lowestTaxonomicLevel=c("None","Phylum","Class","Order","Family","Genus"), fill=list(avgAbundance=0)) %>%
  mutate(lowestTaxonomicLevel = fct_relevel(lowestTaxonomicLevel, c("None","Phylum","Class","Order","Family","Genus")),
         ratio = fct_relevel(ratio, c("1:0","1:1","0:1"))) %>% 
  ggplot() +
  geom_bar(aes(x=ratio, y=avgAbundance, fill=lowestTaxonomicLevel), stat="identity", color="black", size=0.25) +
  scale_fill_brewer(name="Lowest shared\ntaxonomic level", palette="BuPu") +
  # scale_fill_manual(values=c("None"="#808080",
  #                            "Phylum"="#FFC2FF",
  #                            "Class"="#D15CFF",
  #                            "Order"="#8A59EE",
  #                            "Family"="#4F0CD4",
  #                            "Genus"="#04214E")) +
  xlab("Mixture") +
  ylab("Relative abundance") +
  ylim(0,1) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~combo, nrow=1, strip.position="bottom") +
  theme(strip.placement="outside") +
  DEFAULTS.THEME_PRINT
p_lowestTaxaOppositeRelAb_allSubsetsAllCategories
save_plot(paste0(plotdir, "/lowestTaxaOppositeRelAb_allSubsetsAllCategories.png"), 
          p_lowestTaxaOppositeRelAb_allSubsetsAllCategories, base_width=4, base_height=1.75)

# Plot number of each level of lowest shared taxa compared to ASV permutation.
p_lowestTaxaOppositePermuted <- permutedLowestTaxaOpposite %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_point(aes(x=lowestTaxonomicLevel, y=realCount), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedCount), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Count") +
  facet_wrap(~combo, nrow=2) +
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositePermuted
save_plot(paste0(plotdir, "/lowestTaxaOppositePermuted.png"), p_lowestTaxaOppositePermuted, nrow=2, ncol=2)

# Plot relative abundances of each level of lowest shared taxa compared to ASV permutation.
p_lowestTaxaOppositeRelAbPermuted <- permutedTaxaRelAb %>% 
  filter(passage==5 & ratio=="1:1") %>% 
  ggplot() + 
  geom_point(aes(x=lowestTaxonomicLevel, y=realAbundance), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedAbundance), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Relative abundance") +
  facet_wrap(~combo, nrow=2) + 
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeRelAbPermuted
save_plot(paste0(plotdir, "/lowestTaxaOppositeRelAbPermuted.png"), p_lowestTaxaOppositeRelAbPermuted, 
          nrow=2, ncol=2)

# Plot number of each level of lowest shared taxa compared to ASV permutation, all categories.
p_lowestTaxaOppositePermuted_allCategories <- permutedLowestTaxaOpposite_allCategories %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_point(aes(x=lowestTaxonomicLevel, y=realCount), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedCount), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Count") +
  facet_wrap(~combo, nrow=2) +
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositePermuted_allCategories
save_plot(paste0(plotdir, "/lowestTaxaOppositePermuted_allCategories.png"), 
          p_lowestTaxaOppositePermuted_allCategories, nrow=2, ncol=2)

# Plot relative abundances of each level of lowest shared taxa compared to ASV permutation.
p_lowestTaxaOppositeRelAbPermuted_allCategories <- permutedTaxaRelAb_allCategories %>% 
  filter(passage==5 & ratio=="1:1") %>% 
  ggplot() + 
  geom_point(aes(x=lowestTaxonomicLevel, y=realAbundance), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedAbundance), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Relative abundance") +
  facet_wrap(~combo, nrow=2) + 
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeRelAbPermuted_allCategories
save_plot(paste0(plotdir, "/lowestTaxaOppositeRelAbPermuted_allCategories.png"),
          p_lowestTaxaOppositeRelAbPermuted_allCategories, nrow=2, ncol=2)

# Plot number of each level of lowest shared taxa compared to direction permutation.
p_lowestTaxaOppositeDirectionPermuted <- permutedDirectionsLowestTaxaOpposite %>% 
  filter(passage==5 & Phylum %in% c("Bacteroidota","Firmicutes","Proteobacteria")) %>% 
  ggplot() +
  geom_point(aes(x=lowestTaxonomicLevel, y=realCount), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedCount), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Count") +
  facet_grid(~Phylum) + 
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeDirectionPermuted
save_plot(paste0(plotdir, "/lowestTaxaOppositeDirectionPermuted.png"), p_lowestTaxaOppositeDirectionPermuted, 
          nrow=2, ncol=2)

# Plot relative abundance of each level of lowest shared taxa compared to direction permutation.
p_lowestTaxaOppositeRelAbDirectionPermuted <- permutedDirectionsTaxaRelAb %>% 
  filter(passage==5 & ratio=="1:1" & Phylum %in% c("Bacteroidota", "Firmicutes", "Proteobacteria")) %>% 
  ggplot() +
  geom_point(aes(x=lowestTaxonomicLevel, y=realAbundance), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedAbundance), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Relative abundance") +
  facet_grid(~Phylum) +
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeRelAbDirectionPermuted
save_plot(paste0(plotdir, "/lowestTaxaOppositeRelAbDirectionPermuted.png"), p_lowestTaxaOppositeRelAbDirectionPermuted, 
          nrow=2, ncol=2)

# Plot number of each level of lowest shared taxa compared to direction permutation.
# Phylum Bacteroidota, grouped by genus.
p_lowestTaxaOppositeDirectionPermuted_BacteroidotaGenus <- permutedDirectionsLowestTaxaOpposite_BacteroidotaGenus %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_point(aes(x=lowestTaxonomicLevel, y=realCount), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedCount), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Count") +
  facet_grid(~Genus) + 
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeDirectionPermuted_BacteroidotaGenus
save_plot(paste0(plotdir, "/lowestTaxaOppositeDirectionPermuted_BacteroidotaGenus.png"), 
          p_lowestTaxaOppositeDirectionPermuted_BacteroidotaGenus, nrow=2, ncol=2)

# Plot relative abundance of each level of lowest shared taxa compared to direction permutation.
# Phylum Bacteroidota, grouped by genus.
p_lowestTaxaOppositeRelAbDirectionPermuted_BacteroidotaGenus <- permutedDirectionsTaxaRelAb_BacteroidotaGenus %>% 
  filter(passage==5 & ratio=="1:1") %>% 
  ggplot() +
  geom_point(aes(x=lowestTaxonomicLevel, y=realAbundance), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedAbundance), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Relative abundance") +
  facet_grid(~Genus) +
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeRelAbDirectionPermuted_BacteroidotaGenus
save_plot(paste0(plotdir, "/lowestTaxaOppositeRelAbDirectionPermuted_BacteroidotaGenus.png"), 
          p_lowestTaxaOppositeRelAbDirectionPermuted_BacteroidotaGenus, nrow=2, ncol=2)  

# Plot number of each level of lowest shared taxa compared to direction permutation.
# Phylum Firmicutes, grouped by class.
p_lowestTaxaOppositeDirectionPermuted_FirmicutesClass <- permutedDirectionsLowestTaxaOpposite_FirmicutesClass %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_point(aes(x=lowestTaxonomicLevel, y=realCount), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedCount), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Count") +
  facet_grid(~Class) + 
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeDirectionPermuted_FirmicutesClass
save_plot(paste0(plotdir, "/lowestTaxaOppositeDirectionPermuted_FirmicutesClass.png"), 
          p_lowestTaxaOppositeDirectionPermuted_FirmicutesClass, nrow=2, ncol=2)

# Plot relative abundance of each level of lowest shared taxa compared to direction permutation.
# Phylum Firmicutes, grouped by class.
p_lowestTaxaOppositeRelAbDirectionPermuted_FirmicutesClass <- permutedDirectionsTaxaRelAb_FirmicutesClass %>% 
  filter(passage==5 & ratio=="1:1") %>% 
  ggplot() +
  geom_point(aes(x=lowestTaxonomicLevel, y=realAbundance), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedAbundance), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Relative abundance") +
  facet_grid(~Class) +
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeRelAbDirectionPermuted_FirmicutesClass
save_plot(paste0(plotdir, "/lowestTaxaOppositeRelAbDirectionPermuted_FirmicutesClass.png"), 
          p_lowestTaxaOppositeRelAbDirectionPermuted_FirmicutesClass, nrow=2, ncol=2)  

# Plot number of each level of lowest shared taxa compared to direction permutation.
# Phylum Firmicutes, grouped by order.
p_lowestTaxaOppositeDirectionPermuted_FirmicutesOrder <- permutedDirectionsLowestTaxaOpposite_FirmicutesOrder %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_point(aes(x=lowestTaxonomicLevel, y=realCount), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedCount), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Count") +
  facet_wrap(~Order, nrow=2) + 
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeDirectionPermuted_FirmicutesOrder
save_plot(paste0(plotdir, "/lowestTaxaOppositeDirectionPermuted_FirmicutesOrder.png"), 
          p_lowestTaxaOppositeDirectionPermuted_FirmicutesOrder, nrow=2, ncol=2)

# Plot relative abundance of each level of lowest shared taxa compared to direction permutation.
# Phylum Bacteroidota, grouped by genus.
p_lowestTaxaOppositeRelAbDirectionPermuted_FirmicutesOrder <- permutedDirectionsTaxaRelAb_FirmicutesOrder %>% 
  filter(passage==5 & ratio=="1:1") %>% 
  ggplot() +
  geom_point(aes(x=lowestTaxonomicLevel, y=realAbundance), color="orange") +
  geom_boxplot(aes(x=lowestTaxonomicLevel, y=permutedAbundance), color="black", alpha=0) +
  xlab("Lowest shared taxonomic level") +
  ylab("Relative abundance") +
  facet_wrap(~Order, nrow=2) +
  DEFAULTS.THEME_PRES
p_lowestTaxaOppositeRelAbDirectionPermuted_FirmicutesOrder
save_plot(paste0(plotdir, "/lowestTaxaOppositeRelAbDirectionPermuted_FirmicutesOrder.png"), 
          p_lowestTaxaOppositeRelAbDirectionPermuted_FirmicutesOrder, nrow=2, ncol=2)  
