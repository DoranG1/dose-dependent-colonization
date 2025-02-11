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
library(data.table)

# Import data and sample sheet, set up paths. ---------------------------------------------

# Set working directory for running locally.
setwd("../")
# Set the plot theme. 
theme_set(theme_cowplot())
source("config/palettes/plotDefaults.R")
# Set plot output directories.
outdir <- "workflow/out/e0017"
output_path <- file.path(outdir,"DADA2_output")

# Import read data.
ps <- readRDS(file.path(output_path,"ps_all.rds"))
# Import coverage data.
coverage <- read.table(file.path(output_path, "read_summary.txt"), header=TRUE, 
                       col.names=c("filename", "readsIn", "readsOut")) %>%
  mutate(filename=gsub("_.*", "", filename))
# Import samplesheet.
ss <- read.table("config/210819-e0015-e0020-16S-samplesheet.txt", header=TRUE, sep="\t") %>% 
  filter(group=="e0017") %>%
  select(-c(sample, group, plate, media, 
            cultureID, stockID, bioTyperDate, bioTyperClassification))
# Import taxa color palette.
palette <- read.table("config/palettes/KCHcolors-Silva-partial.txt",
                      header=TRUE, stringsAsFactors = FALSE) %>% 
  mutate(taxaShort=gsub("*.*\\.","",taxa))

# e0017 global variables
combos <- c("XBA-XBB", "XCA-XCB", "XDA-XDB", "XFA-XFB", "XBA-XCA", "XCA-XDA", "XDA-XFA", "XFA-XBA")
passages <- c(3,5)
ratios <- c("1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000")
ratiosFull <- c("1:0", "1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000", "0:1")

# Create dataframe. ------------------------------

# Extract OTU table and taxonomy tables from phyloseq object.
df_otu_raw <- as.data.frame(otu_table(ps))
df_tax_raw <- as.data.frame(tax_table(ps))

# Create a column containing the sample names in each table.
df_otu <- df_otu_raw %>%
  mutate(sample=rownames(df_otu_raw))
df_tax <- df_tax_raw %>%
  mutate(OTU=rownames(df_tax_raw))

# Tidy the OTU data frame.
df_otu <- df_otu %>%
  gather(OTU, count, -sample)

# Combine ASV and taxonomy dataframes. Split full sample string into sample and filename strings. 
df <- left_join(df_otu, df_tax) %>%
  mutate(filename=gsub("_.*", "", sample)) %>%
  select(-sample) %>%
  left_join(ss, by="filename")

# Generate a list of unique IDs for each sample for SRA upload.
dfSRA <- left_join(df_otu, df_tax) %>% 
  mutate(filename=gsub("_.*", "", sample)) %>%
  left_join(ss, by="filename") %>% 
  mutate(sample_ID = sample) %>% 
  filter(!is.na(biosample1) | !is.na(biosample2)) %>%
  filter(extractionReplicate==1) %>% 
  mutate(sample_name = ifelse(!is.na(ratio), 
                              paste0(sub("-.*","",biosample1),"_",sub("-.*","",biosample2),"_",sub(":","-",ratio),"_p",passage,"_rep",inoculationReplicate),
                              paste0(sub("-.*","",biosample1),"_unmixed_p",passage,"_rep",inoculationReplicate))) %>% 
  select(sample_ID, sample_name, biosample1, biosample2, ratio, passage, inoculationReplicate) %>% 
  unique()
write.table(dfSRA, "config/SRAsamples.txt", row.names=FALSE, quote=FALSE, sep="\t")

# Remove timepoint from sample name, group samples into control or mixture.
df <- df %>%
  mutate(biosample1=gsub("-.*", "", biosample1), biosample2=gsub("-.*", "", biosample2),
         group=case_when(is.na(biosample1) & is.na(biosample2) ~"blank", 
                         !is.na(biosample1) & is.na(biosample2) ~ "control",
                         !is.na(biosample1) & !is.na(biosample2) ~ "mixture")) %>%
  filter(group!="blank")

# Remove low-coverage samples.
df <- df %>%
  left_join(coverage %>% select(c(filename, readsOut)), by="filename")
  filter(readsOut > 10000)

# Generate labels
df_mixtures <- foreach(x=combos, .combine=rbind) %do% {
  sample1 <- gsub("-.*", "", x)
  sample2 <- gsub(".*-", "", x)
  dfTemp <- df %>%
    filter(biosample1 %in% c(sample1, sample2) & biosample2 %in% c(sample1, sample2)) %>%
    mutate(ratio1=gsub(":.*", "", ratio), ratio2=gsub(".*:", "", ratio)) %>%
    mutate(label = case_when(biosample1==sample1 & !is.na(biosample2) ~ paste(sample1, ratio, sample2, sep="-"),
                             biosample1!=sample1 & !is.na(biosample2) ~ paste(sample1, paste0(ratio2, ":", ratio1), 
                                                                              sample2, sep="-"))) %>%
    mutate(ratio=sub(".*?-", "", label), ratio=sub("-.*", "", ratio)) %>%
    select(-c(ratio1, ratio2))
  return(dfTemp)
}
df <- rbind(df %>%
              filter(is.na(biosample2)) %>%
              mutate(label=biosample1),
            df_mixtures) %>%
  mutate(mixtureType="actual")

# Calculate relative abundance and join to dataframe.
df <- df %>% 
  group_by(filename) %>%
  mutate(rel_abundance=count/sum(count)) %>% 
  ungroup()

# Create stripped-down color palette containing only families present in dataset.
palette <- palette %>% 
  # Filter to family present at high cumulative relative abundance.
  filter(taxaShort %in% sort(unique(df %>% 
                                      group_by(passage, label, inoculationReplicate, extractionReplicate, mixtureType, Family) %>% 
                                      mutate(familyRelAbundance = sum(rel_abundance, na.rm=TRUE)) %>% 
                                      filter(familyRelAbundance > 0.01) %>% 
                                      pull(Family))))
  # Filter to ASV in family present at high relative abundance.
  #filter(taxaShort %in% sort(unique(df %>% filter(rel_abundance > 0.01) %>% pull(Family))))
# Make a named list.
paletteVector <- palette$hex
names(paletteVector) <- palette$taxaShort

# Create dataframe with no zero abundances.
dfAbundant <- df %>%
  filter(count > 0)

# Create dataframe with >0.1% relative abundance.
dfRelAbundant <- df %>%
  filter(rel_abundance>0.001)

# Generate theoretical mixtures. ---------------------------------------------------

comboPassages <- expand.grid(combos,passages) %>% mutate(comboPassage=paste(Var1,Var2,sep="_"))
comboPassages <- comboPassages$comboPassage

# Function for generating a mini-dataframe of relative abundances for a specific theoretical mixture at one ratio.
addMixtureRatio <- function(ratio, otu_df, sample1, sample2) {
  ratio1 <- as.integer(gsub(":.*", "", ratio))
  ratio2 <- as.integer(gsub(".*:", "", ratio))
  otu_df <- otu_df %>%
    group_by(inoculationReplicate) %>%
    mutate(rel_abundance = (A_rel_abundance/ratio2) + (B_rel_abundance/ratio1),
           ratio=ratio,
           label=paste(sample1, ratio, sample2, sep="-"))
  return(otu_df)
}

# Function for generating the full theoretical mixture dataframe with all controls, ratios, replicates for a specific mixture.
generateTheoreticalMixture <- function(comboPassage) {
  # Parse combo, samples, and passage from given comboPassage parameter.
  combo <- gsub("_.*", "", comboPassage)
  sample1 <- gsub("-.*", "", combo)
  sample2 <- gsub(".*-", "", combo)
  labelOrder <- paste(sample1, ratiosFull, sample2, sep="-")
  passageNum <- gsub(".*_", "", comboPassage)
  A_otu <- paste0(sample1, "_otu")
  assign(A_otu, df %>%
           filter(extractionReplicate==1 & biosample1==sample1 & is.na(biosample2) & passage==passageNum))
  B_otu <- paste0(sample2, "_otu")
  assign(B_otu, df %>%
           filter(extractionReplicate==1 & biosample1==sample2 & is.na(biosample2) & passage==passageNum))
  # Generate dataframe of all OTUs to be combined, annotated by presence in each subject.
  mix_otus <- full_join(get(A_otu) %>%
                          select(c(OTU, rel_abundance, inoculationReplicate, passage)),
                        get(B_otu) %>% 
                          select(c(OTU, rel_abundance, inoculationReplicate, passage)),
                        by=c("OTU", "inoculationReplicate", "passage")) %>%
    left_join(df_tax %>% select(!Kingdom), by="OTU") %>%
    rename(A_rel_abundance=rel_abundance.x, B_rel_abundance=rel_abundance.y) %>%
    mutate(A_rel_abundance=replace(A_rel_abundance, is.na(A_rel_abundance), 0),
           B_rel_abundance=replace(B_rel_abundance, is.na(B_rel_abundance), 0))
  # Generate all theoretical mixtures.
  theoretical_mixture <- lapply(ratios, FUN=addMixtureRatio,
                                otu_df=mix_otus, sample1=sample1, sample2=sample2) %>%
    list.rbind() %>% 
    # Add relative abundances back to mixtures.
    group_by(inoculationReplicate, ratio) %>%
    mutate(rel_abundance=rel_abundance/sum(rel_abundance), group="mixture") %>%
    ungroup() %>%
    select(-c(A_rel_abundance, B_rel_abundance)) %>%
    # Add controls.
    rbind(get(A_otu) %>%
            select(c(OTU, Phylum, Class, Order, Family, Genus, rel_abundance, inoculationReplicate,
                     passage, ratio, label, group)) %>%
            mutate(ratio="1:0", label=paste(sample1, "1:0", sample2, sep="-")),
          get(B_otu) %>%
            select(c(OTU, Phylum, Class, Order, Family, Genus, rel_abundance, inoculationReplicate,
                     passage, ratio, label, group)) %>%
            mutate(ratio="0:1", label=paste(sample1, "0:1", sample2, sep="-"))
    ) %>%
    mutate(mixtureType="theoretical") %>%
    # Bind actual to theoretical mixture.
    rbind(df %>%
            filter(extractionReplicate==1 & passage==passageNum &
                     biosample1 %in% c(sample1, sample2) & biosample2 %in% c(sample1, sample2)) %>%
            select(c(OTU, Phylum, Class, Order, Family, Genus, rel_abundance, inoculationReplicate,
                     passage, ratio, label, group, mixtureType))
    ) %>%
    mutate(combo=combo,
           line=paste(mixtureType, inoculationReplicate, sep="-")) %>%
    # Reorder samples in correct label order.
    ungroup() %>%
    mutate(label=fct_relevel(label, labelOrder)) %>%
    # Add grouping for origin of OTUs.
    mutate(parent=NA,
           parent=replace(parent,
                          is.na(parent) &
                            OTU %in% get(A_otu)$OTU &
                            OTU %in% get(B_otu)$OTU,
                          "both"),
           parent=replace(parent,
                          is.na(parent) &
                            OTU %in% get(A_otu)$OTU,
                          "A"),
           parent=replace(parent,
                          is.na(parent) &
                            OTU %in% get(B_otu)$OTU,
                          "B"),
           parent=replace(parent, is.na(parent), "neither"))
  return(theoretical_mixture)
}

# Generate all 8 theoretical mixtures + actual mixtures in a list.
theoretical_mixtures_unfiltered <- lapply(comboPassages, generateTheoreticalMixture)

# Add all 8 theoretical mixtures + actual mixtures to a dataframe.
theoretical_mixtures_df_unfiltered <- list.rbind(theoretical_mixtures_unfiltered) %>% 
  # Remove poorly sequenced XBA replicate.
  filter(!(combo %in% c("XBA-XBB", "XBA-XCA", "XFA-XBA") & 
             inoculationReplicate==1 & passage==3 & mixtureType=="theoretical") &
           rel_abundance>0) %>%
  # Reorder some columns for future sorting.
  ungroup() %>%
  mutate(line=fct_relevel(line, c("theoretical-1", "theoretical-2", "actual-1", "actual-2", "actual-3")),
         combo=fct_relevel(combo, combos),
         ratio=fct_relevel(ratio, c("1:0","1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000", "0:1")),
         logAbundance=log10(rel_abundance)) %>% 
  # Add family level relative abundances.
  group_by(combo, passage, Family, line, ratio) %>% 
  mutate(family_rel_abundance=sum(rel_abundance),
         logFamilyAbundance=log10(family_rel_abundance))

# Annotate dataframe with theoretical communities. ----------------------------------

theoretical_mixtures_df <- theoretical_mixtures_df_unfiltered %>% 
  # Remove replicates that sequenced poorly.
  filter(line!="theoretical-3" & 
           !(line=="theoretical-1" & combo %in% c("XBA-XBB","XBA-XCA", "XFA-XBA") & passage==3) &
           !(line=="actual-1" & combo=="XFA-XBA" & passage==5 & ratio=="100:1")) %>%
  # Annotate abundances.
  mutate(rel_abundance_old = rel_abundance) %>% 
  # Add pseudopoints at 10^-3 rel abundance.
  group_by(combo, passage) %>%
  complete(nesting(OTU, Family),
           line=c("theoretical-1","theoretical-2","actual-1","actual-2","actual-3"),
           ratio=ratiosFull,
           fill=list(logAbundance=-3, rel_abundance=10^-3)) %>% 
  # Set abundances of poorly-sequenced XFA-XBA 100:1 replicate to NA to keep them in plots.
  mutate(logAbundance = replace(
    logAbundance, passage==5 & combo=="XFA-XBA" & ratio=="100:1" & line=="actual-1", NA),
    rel_abundance = replace(
      rel_abundance, passage==5 & combo=="XFA-XBA" & ratio=="100:1" & line=="actual-1", NA)) %>% 
  ungroup() %>% 
  # Remove pseudopoints at parent communities.
  filter(!(ratio %in% c("1:0","0:1")) | line %in% c("theoretical-1","theoretical-2")) %>% 
  # Add family-level relative abundances from before filtering.
  select(-c(family_rel_abundance, logFamilyAbundance)) %>% 
  left_join(theoretical_mixtures_df_unfiltered %>% 
              group_by(Family, combo, passage, line, ratio) %>% 
              summarize(logFamilyAbundance, family_rel_abundance) %>%
              unique()) %>% 
  # Redo certain columns that were removed by complete().
  # Set relative abundances lower than 10^-3 to 10^-3.
  mutate(ratio=fct_relevel(ratio, ratiosFull),
         inoculationReplicate=sub(".*-","",line),
         mixtureType=sub("-.*","",line),
         logAbundance=replace(logAbundance, logAbundance<=-3, -3),
         rel_abundance=replace(rel_abundance, rel_abundance<=10^-3, 10^-3),
         logFamilyAbundance=replace(
           logFamilyAbundance, logFamilyAbundance<=-3 | is.na(logFamilyAbundance), -3),
         family_rel_abundance=replace(
           family_rel_abundance, family_rel_abundance<=10^-3 | is.na(family_rel_abundance), 10^-3),
         pointType=if_else(logAbundance!=-3, "circle", "openCircle"),
         familyPointType=if_else(logFamilyAbundance!=-3, "circle", "openCircle"),
         group=if_else(ratio %in% c("1:0","0:1"), "control", "mixture"),
         label=paste(sub("-.*","",combo), ratio, sub(".*-","",combo), sep="-"),
         filterRelAb3 = rel_abundance >= 10^-3) %>% 
  # Determine parent origin for all ASVs above 10^-3 relative abundance.
  group_by(combo, passage, OTU) %>% 
  mutate(parent1 = sum(rel_abundance[ratio=="1:0" & inoculationReplicate==1]) > 10^-3 | 
           sum(rel_abundance[ratio=="1:0" & inoculationReplicate==2]) > 10^-3,
         parent2 = sum(rel_abundance[ratio=="0:1" & inoculationReplicate==1]) > 10^-3 |
           sum(rel_abundance[ratio=="0:1" & inoculationReplicate==2]) > 10^-3,
         origin = case_when(parent1 & parent2 ~ "both",
                            parent1 & !parent2 ~ "A",
                            !parent1 & parent2 ~ "B",
                            !parent1 & !parent2 ~ "neither")) %>% 
  select(-c(parent1, parent2)) %>% 
  ungroup() %>% 
  # Add ASV counts above 10^-3 relative abundance.
  left_join(theoretical_mixtures_df_unfiltered %>% 
              filter(rel_abundance >= 10^-3) %>% 
              group_by(passage, label, line) %>% 
              summarize(ASVcountRelAb3 = n_distinct(OTU))) %>% 
  # Calculate difference between theoretical parent communities.
  group_by(combo, passage, OTU, line) %>%
  mutate(leftParent = ifelse(ratio=="1:0", logAbundance, ifelse(mixtureType=="theoretical", 0, NA)),
         rightParent = ifelse(ratio=="0:1", logAbundance, ifelse(mixtureType=="theoretical", 0, NA)),
         leftParent = ifelse(sum(leftParent)!=0, sum(leftParent), -3),
         rightParent = ifelse(sum(rightParent)!=0, sum(rightParent), -3),
         highAbundanceParentLeft = leftParent >= -2.6,
         highAbundanceParentRight = rightParent >= -2.6,
         parentDiff = rightParent - leftParent) %>% 
  ungroup() %>% 
  group_by(combo, passage, OTU) %>% 
  mutate(parentDiffAvg = mean(unique(parentDiff), na.rm=TRUE),
         # Replace NA families with character "NA" to add OTU numbers.
         Family=replace(Family, is.na(Family), "NA")) %>% 
  ungroup()

# Add OTU numbers to OTUs in each family.
theoretical_mixtures_df <- foreach(x=unique(theoretical_mixtures_df$Family), .combine=rbind) %do% {
  dfTemp <- theoretical_mixtures_df %>% 
    filter(Family==x) %>% 
    group_by(OTU) %>%
    mutate(OTUnum=cur_group_id())
  return(dfTemp)
} 

# Summary statistics. -------------------------------------------------

# Calculate average theoretical sd and range of ASV/family log abundances across all ratios.
t_sum <- theoretical_mixtures_df %>%
  filter(mixtureType=="theoretical") %>%
  # Remove poorly sequenced XBA replicate from passage 3.
  filter(!(passage==3 & combo %in% c("XBA-XBB","XBA-XCA","XFA-XBA") & line=="theoretical-1")) %>% 
  group_by(Family, label, passage) %>%
  mutate(familyRangeAbundance=max(logFamilyAbundance, na.rm=TRUE)-min(logFamilyAbundance, na.rm=TRUE),
         familyAvgAbundance=mean(logFamilyAbundance, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(OTU, label, passage) %>% 
  mutate(rangeAbundance=max(logAbundance, na.rm=TRUE)-min(logAbundance, na.rm=TRUE),
         avgAbundance=mean(logAbundance, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(OTU, combo, line, passage) %>% 
  arrange(line, ratio) %>% 
  mutate(doseDist=abs(logAbundance-lag(logAbundance)),
         doseDistTotal=mean(doseDist, na.rm=TRUE),
         doseDiff=logAbundance-lag(logAbundance),
         doseDiffTotal=mean(doseDiff, na.rm=TRUE),
         leftAbundance = ifelse(ratio=="1000:1", logAbundance, 0),
         rightAbundance = ifelse(ratio=="1:1000", logAbundance, 0),
         leftAbundance = sum(leftAbundance),
         rightAbundance = sum(rightAbundance),
         doseSumDist = abs(rightAbundance - leftAbundance),
         doseSumDiff = rightAbundance - leftAbundance,
         familyDoseDist=abs(logFamilyAbundance-lag(logFamilyAbundance)),
         familyDoseDistTotal=mean(familyDoseDist, na.rm=TRUE),
         familyDoseDiff=logFamilyAbundance-lag(logFamilyAbundance),
         familyDoseDiffTotal=mean(familyDoseDiff, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(OTU, combo, passage) %>%
  summarize(Family, OTUnum, ratio, group, label, avgAbundance, rangeAbundance,
            t_rangeAvg=mean(rangeAbundance),
            t_doseDistance=mean(doseDistTotal), 
            t_doseDifference=mean(doseDiffTotal),
            t_doseSumDistance=mean(doseSumDist),
            t_doseSumDifference=mean(doseSumDiff),
            familyAvgAbundance, familyRangeAbundance,
            t_familyRangeAvg=mean(familyRangeAbundance),
            t_familyDoseDistance=mean(familyDoseDistTotal), 
            t_familyDoseDifference=mean(familyDoseDiffTotal)) %>%
  distinct()

# Calculate average actual sd and range of log abundances across all ratios.
a_sum <- theoretical_mixtures_df %>%
  filter(mixtureType=="actual" & group=="mixture") %>%
  group_by(Family, label, passage) %>%
  mutate(familyRangeAbundance=max(logFamilyAbundance, na.rm=TRUE)-min(logFamilyAbundance, na.rm=TRUE),
         familyAvgAbundance=mean(logFamilyAbundance, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(OTU, label, passage) %>% 
  mutate(rangeAbundance=max(logAbundance, na.rm=TRUE)-min(logAbundance, na.rm=TRUE),
         avgAbundance=mean(logAbundance, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(OTU, combo, line, passage) %>% 
  arrange(line, ratio) %>% 
  # Correct for missing replicate in XFA-XBA mixture.
  mutate(logAbundance=ifelse(
    passage==5 & combo=="XFA-XBA" & line=="actual-1" & ratio=="100:1", lag(logAbundance), logAbundance)) %>% 
  mutate(doseDist=ifelse(
    !(passage==5 & combo=="XFA-XBA" & line=="actual-1" & ratio=="100:1"), abs(logAbundance-lag(logAbundance)), NA),
    doseDistTotal=mean(doseDist, na.rm=TRUE),
    doseDiff=ifelse(
      !(passage==5 & combo=="XFA-XBA" & line=="actual-1" & ratio=="100:1"), logAbundance-lag(logAbundance), NA),
    doseDiffTotal=mean(doseDiff, na.rm=TRUE),
    leftAbundance = ifelse(ratio=="1000:1", logAbundance, 0),
    rightAbundance = ifelse(ratio=="1:1000", logAbundance, 0),
    leftAbundance = sum(leftAbundance),
    rightAbundance = sum(rightAbundance),
    doseSumDist = abs(rightAbundance - leftAbundance),
    doseSumDiff = rightAbundance - leftAbundance,
    familyDoseDist=abs(logFamilyAbundance-lag(logFamilyAbundance)),
    familyDoseDistTotal=mean(familyDoseDist, na.rm=TRUE),
    familyDoseDiff=logFamilyAbundance-lag(logFamilyAbundance),
    familyDoseDiffTotal=mean(familyDoseDiff, na.rm=TRUE)) %>% 
  ungroup() %>%
  group_by(OTU, combo, passage) %>%
  summarize(Family, OTUnum, ratio, group, label, avgAbundance, rangeAbundance,
            a_rangeAvg=mean(rangeAbundance),
            a_doseDistance=mean(doseDistTotal), 
            a_doseDifference=mean(doseDiffTotal),
            a_doseSumDistance=mean(doseSumDist),
            a_doseSumDifference=mean(doseSumDiff),
            familyAvgAbundance, familyRangeAbundance,
            a_familyRangeAvg=mean(familyRangeAbundance),
            a_familyDoseDistance=mean(familyDoseDistTotal), 
            a_familyDoseDifference=mean(familyDoseDiffTotal)) %>%
  distinct()

# Calculate differences between theoretical and actual abundances for each ratio and averages.
summary_df <- full_join(t_sum, a_sum, by=c("OTU", "Family", "OTUnum", "passage", "combo", "ratio", "group", "label")) %>% 
  rename(t_abundance=avgAbundance.x, a_abundance=avgAbundance.y,
         t_rangeAbundance=rangeAbundance.x, a_rangeAbundance=rangeAbundance.y,
         t_familyAbundance=familyAvgAbundance.x, a_familyAbundance=familyAvgAbundance.y,
         t_familyRangeAbundance=familyRangeAbundance.x, a_familyRangeAbundance=familyRangeAbundance.y) %>%
  mutate(difference=a_abundance-t_abundance, 
         distance=abs(difference),
         familyDifference=a_familyAbundance-t_familyAbundance,
         familyDistance=abs(familyDifference)) %>%
  group_by(OTU, combo, passage) %>%
  mutate(avgdifference=mean(difference, na.rm=TRUE),
         avgdistance=mean(distance, na.rm=TRUE),
         familyAvgDifference=mean(familyDifference, na.rm=TRUE), 
         familyAvgDistance=mean(familyDistance, na.rm=TRUE))

# Pull out specific summary statistics for each ASV/combo.
ASV_sum <- summary_df %>%
  group_by(OTU, Family, OTUnum, combo, passage) %>%
  summarize(t_rangeAvg, a_rangeAvg, 
            t_doseDistance, t_doseDifference, t_doseSumDistance, t_doseSumDifference,
            a_doseDistance, a_doseDifference, a_doseSumDistance, a_doseSumDifference,
            avgdistance, avgdifference, 
            t_familyRangeAvg, a_familyRangeAvg, 
            t_familyDoseDistance, t_familyDoseDifference, a_familyDoseDistance, a_familyDoseDifference, 
            familyAvgDifference, familyAvgDistance) %>% 
  distinct() %>% 
  # Remove NAs caused by inclusion of parent communities.
  filter(!is.na(a_rangeAvg))

# Group ASVs into behavior categories.
theoretical_mixtures_df_categorized <- theoretical_mixtures_df %>% 
  left_join(ASV_sum, by=c("OTU","Family","OTUnum","combo","passage")) %>% 
  mutate(actualLogAbundance = ifelse(mixtureType=="actual", logAbundance, NA),
         absDoseDifference = ifelse(abs(a_doseDifference) < 0.1, 0.1, abs(a_doseDifference))) %>%
  group_by(combo, passage, OTU) %>%
  mutate(category = ifelse(
    # High abundance?
    TRUE %in% highAbundanceParentLeft | TRUE %in% highAbundanceParentRight, ifelse(
    #mean(actualLogAbundance, na.rm = TRUE) > -2.6, ifelse(
      # Low noise?
      a_doseDistance / absDoseDifference <= 2.25, ifelse(
        # High dose-difference?
        abs(a_doseDifference) > 0.09, "dose-dependent",
        ifelse(
          # High dose distance?
          #a_doseDistance > 0.3, ifelse(
            # Theoretical difference positive or negative?
            #avgdifference > 0, "DD-over", "DD-under"
          #), ifelse(
            # High theoretical distance?
            avgdistance > 0.3 & abs(avgdifference) >= 0.2, ifelse(
              # Theoretical difference positive or negative?
              avgdifference > 0, "overcolonizing", "undercolonizing"
            ), "flat"
          #)
        )
      ), "noisy"
    ), "lowAbundance"))

# Export dataframe.
#write.table(theoretical_mixtures_df_categorized, "analysis/mixtureDataframe.txt",
#            quote=FALSE, row.names=FALSE, sep="\t")

# Generate new list of all unique families in dataframe to plot over.
familyList <- na.exclude(unique(theoretical_mixtures_HighAbundanceOneParent$Family))

# Create theoretical phyloseq. --------------------------------------------

# Generate OTU table with OTUs as rownames.
# Convert relative abundance back to counts: scale relative abundance to median counts across all communities.
df_otu_raw_theoretical <- theoretical_mixtures_df %>% 
  group_by(passage, label, line) %>% 
  mutate(rel_abundance_old = replace(rel_abundance_old, is.na(rel_abundance_old), 0),
         count=ceiling(rel_abundance_old * 131515)) %>% 
  ungroup() %>% 
  group_by() %>% 
  mutate(sample=paste(line, label, passage, sep="_")) %>% 
  select(OTU, count, sample) %>% 
  pivot_wider(names_from=sample, values_from=count)
df_otu_raw_theoretical[is.na(df_otu_raw_theoretical)] <- 0
otu_rownames <- df_otu_raw_theoretical$OTU
df_otu_raw_theoretical <- df_otu_raw_theoretical %>% select(!OTU)
rownames(df_otu_raw_theoretical) <- otu_rownames

# Generate sample metadata using new sample names.
df_sam_theoretical <- theoretical_mixtures_df %>% 
  ungroup() %>% 
  select(inoculationReplicate, passage, ratio, label, group, mixtureType, combo, line) %>% 
  unique() %>% 
  mutate(sample=paste(line, label, passage, sep="_"))
sample_rownames <- df_sam_theoretical$sample
sam_theoretical <- df_sam_theoretical %>% select(!sample)
sam_theoretical <- sample_data(df_sam_theoretical)
sample_names(sam_theoretical) <- sample_rownames

# Create theoretical phyloseq.
ps_theoretical <- phyloseq(otu_table(as.matrix(df_otu_raw_theoretical), taxa_are_rows=TRUE),
                           tax_table(as.matrix(df_tax_raw)),
                           sam_theoretical)

# Export theoretical phyloseq as individual text files.
#write.table(df_otu_raw_theoretical, "analysis/OTUtableTheoretical.txt", quote=FALSE, row.names=TRUE, sep="\t")
#write.table(sam_theoretical, "analysis/samTableTheoretical.txt", quote=FALSE, row.names=FALSE, sep="\t")
