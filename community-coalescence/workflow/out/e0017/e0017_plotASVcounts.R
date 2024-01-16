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

# Import data and sample sheet, set up paths. ---------------------------------------------

# Set working directory for running locally.
setwd("/oak/stanford/groups/relman/users/dorang/210827-e0017-16S")
# Set the plot theme.
theme_set(theme_cowplot())
source("config/palettes/plotDefaults.R")
# Set plot output directories.
outdir <- "workflow/out/e0017"
output_path <- file.path(outdir,"DADA2_output")

# Import read data.
ps <- readRDS(file.path(output_path,"ps_taxa.rds"))
# Import coverage data.
coverage <- read.table(file.path(output_path, "read_summary.txt"), header=TRUE, 
                       col.names=c("filename", "readsIn", "readsOut")) %>%
  mutate(filename=gsub("_.*", "", filename))
# Import samplesheet.
ss <- read_tsv("config/210819-e0015-e0020-16S-samplesheet.txt") %>%
  filter(group=="e0017") %>%
  select(-c(sample, group, plate, media, 
            cultureID, stockID, bioTyperDate, bioTyperClassification))
# Import taxa color palette.
palette <- read.table("config/palettes/KCHcolors-Silva-partial.txt",
                      header=TRUE, stringsAsFactors = FALSE)

# e0017 global variables
combos <- c("XBA-XBB", "XCA-XCB", "XDA-XDB", "XFA-XFB", "XBA-XCA", "XCA-XDA", "XDA-XFA", "XFA-XBA")
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

# Remove timepoint from sample name, group samples into control or mixture.
df <- df %>%
  mutate(biosample1=gsub("-.*", "", biosample1), biosample2=gsub("-.*", "", biosample2),
         group=case_when(is.na(biosample1) & is.na(biosample2) ~"blank", 
                         !is.na(biosample1) & is.na(biosample2) ~ "control",
                         !is.na(biosample1) & !is.na(biosample2) ~ "mixture")) %>%
  filter(group!="blank")

# Remove low-coverage samples.
df <- df %>%
  left_join(coverage %>% select(c(filename, readsOut)), by="filename") %>%
  filter(readsOut > 10000)

# Generate labels
df_mixtures <- foreach(x=combos, .combine=rbind) %dopar% {
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
  left_join(df %>%
              group_by(filename) %>%
              summarize(total_counts=sum(count))) %>%
  mutate(rel_abundance=count/total_counts)

# Bind taxa colors by family.
df <- left_join(df, 
                palette %>% dplyr::select(taxashort, hex) %>%
                  filter(taxashort %in% df$Family) %>% 
                  rename(Family=taxashort), 
                by="Family")

# Extract color palette. Convert unknown taxa to dark gray.
taxaPalette <- df %>% ungroup() %>%
  dplyr::select(Family, hex) %>%
  unique() %>% mutate(hex=ifelse(is.na(hex),"#615c5c", hex)) %>%
  arrange(Family)
taxaPaletteList <- list(taxaPalette$hex)[[1]]
names(taxaPaletteList) <- taxaPalette$Family

# Create dataframe with no zero abundances.
dfAbundant <- df %>%
  filter(count > 0)

# Generate theoretical mixtures. ---------------------------------------------------

# Function for generating a mini-dataframe of relative abundances for a specific theoretical mixture at one ratio.
addMixtureRatio <- function(ratio, otu_df, sample1, sample2) {
  ratio1 <- as.integer(gsub(":.*", "", ratio))
  ratio2 <- as.integer(gsub(".*:", "", ratio))
  otu_df <- otu_df %>%
    group_by(inoculationReplicate) %>%
    mutate(rel_abundance = (A_abundance/ratio2) + (B_abundance/ratio1),
           ratio=ratio,
           label=paste(sample1, ratio, sample2, sep="-"))
  return(otu_df)
}

# Function for generating the full theoretical mixture dataframe with all controls, ratios, replicates for a specific mixture.
generateTheoreticalMixture <- function(combo) {
  sample1 <- gsub("-.*", "", combo)
  sample2 <- gsub(".*-", "", combo)
  labelOrder <- c(paste(sample1, "1:0", sample2, sep="-"), 
                  paste(sample1, "1000:1", sample2, sep="-"), 
                  paste(sample1, "100:1", sample2, sep="-"),
                  paste(sample1, "10:1", sample2, sep="-"), 
                  paste(sample1, "1:1", sample2, sep="-"),
                  paste(sample1, "1:10", sample2, sep="-"), 
                  paste(sample1, "1:100", sample2, sep="-"),
                  paste(sample1, "1:1000", sample2, sep="-"), 
                  paste(sample1, "0:1", sample2, sep="-"))
  A_otu <- paste0(sample1, "_otu")
  assign(A_otu, dfAbundant %>%
           filter(extractionReplicate==1 & biosample1==sample1 & is.na(biosample2) & passage==3))
  B_otu <- paste0(sample2, "_otu")
  assign(B_otu, dfAbundant %>%
           filter(extractionReplicate==1 & biosample1==sample2 & is.na(biosample2) & passage==3))
  # Generate dataframe of all OTUs to be combined, annotated by presence in each subject.
  mix_otus <- full_join(get(A_otu) %>%
                          select(c(OTU, rel_abundance, inoculationReplicate)),
                        get(B_otu) %>% 
                          select(c(OTU, rel_abundance, inoculationReplicate)),
                        by=c("OTU", "inoculationReplicate")) %>%
    left_join(df_tax, by="OTU") %>%
    rename(A_abundance=rel_abundance.x, B_abundance=rel_abundance.y) %>%
    mutate(A_abundance=replace(A_abundance, is.na(A_abundance), 0),
           B_abundance=replace(B_abundance, is.na(B_abundance), 0))
  # Generate all theoretical mixtures.
  theoretical_mixture <- lapply(ratios, FUN=addMixtureRatio,
                                otu_df=mix_otus, sample1=sample1, sample2=sample2) %>%
    list.rbind() %>%
    mutate(group="mixture") %>%
    select(-c(A_abundance, B_abundance)) %>%
    # Add controls.
    rbind(get(A_otu) %>%
            select(c(OTU, Phylum, Class, Order, Family, Genus, rel_abundance, inoculationReplicate, ratio, label, group)) %>%
            mutate(ratio="1:0", label=paste(sample1, "1:0", sample2, sep="-")),
          get(B_otu) %>%
            select(c(OTU, Phylum, Class, Order, Family, Genus, rel_abundance, inoculationReplicate, ratio, label, group)) %>%
            mutate(ratio="0:1", label=paste(sample1, "0:1", sample2, sep="-"))
    ) %>%
    mutate(mixtureType="theoretical") %>%
    # Bind actual to theoretical mixture.
    rbind(dfAbundant %>%
            filter(extractionReplicate==1 & passage==3 &
                     biosample1 %in% c(sample1, sample2) & biosample2 %in% c(sample1, sample2)) %>%
            select(c(OTU, Phylum, Class, Order, Family, Genus, rel_abundance, inoculationReplicate, ratio, label, group, mixtureType)),
          get(A_otu) %>%
            select(c(OTU, Phylum, Class, Order, Family, Genus, rel_abundance, inoculationReplicate, ratio, label, group)) %>%
            mutate(mixtureType="actual", ratio="1:0", label=paste(sample1, "1:0", sample2, sep="-")),
          get(B_otu) %>%
            select(c(OTU, Phylum, Class, Order, Family, Genus, rel_abundance, inoculationReplicate, ratio, label, group)) %>%
            mutate(mixtureType="actual", ratio="0:1", label=paste(sample1, "0:1", sample2, sep="-"))
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
           parent=replace(parent, is.na(parent), "neither")) %>%
    # Add ASV counts.
    group_by(label, line) %>%
    mutate(total_rel=sum(rel_abundance), rel_abundance=rel_abundance/total_rel, ASVcount=length(unique(OTU))) %>%
    select(!total_rel)
  return(theoretical_mixture)
}

# Generate all 8 theoretical mixtures + actual mixtures in a list.
theoretical_mixtures_unfiltered <- lapply(combos, generateTheoreticalMixture)

# Add all 8 theoretical mixtures + actual mixtures to a dataframe.
theoretical_mixtures_df_unfiltered <- list.rbind(theoretical_mixtures_unfiltered) %>%
  # Remove poorly sequenced XBA replicate.
  filter(!(combo %in% c("XBA-XBB", "XBA-XCA", "XFA-XBA") & inoculationReplicate==1)) %>%
  # Reorder some columns for future sorting.
  ungroup() %>%
  mutate(line=fct_relevel(line, c("theoretical-1", "theoretical-2", "actual-1", "actual-2", "actual-3")),
         combo=fct_relevel(combo, combos),
         ratio=fct_relevel(ratio, c("1:0","1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000", "0:1")))

# Create filtered versions of dataframe. ----------------------------------

# Filter to 10^-4 relative abundance minimum.
theoretical_mixtures_df_bottomed <- theoretical_mixtures_df_unfiltered %>%
  ungroup() %>%
  group_by(label, line) %>%
  filter(rel_abundance >= 10^-4) %>%
  mutate(total_rel=sum(rel_abundance), rel_abundance=rel_abundance/total_rel,
         ASVcount=n_distinct(OTU)) %>%
  select(!total_rel)

# Filter out low quality ASVs.
theoretical_mixtures_df_filtered <- theoretical_mixtures_df_bottomed %>%
  # Remove ASV replicates with missing ratios in the middle.
  ungroup() %>%
  mutate(keep=TRUE) %>%
  arrange(OTU, factor(combo), factor(line), factor(ratio)) %>%
  group_by(OTU, combo, line) %>%
  rowwise() %>%
  mutate(ratioOrder=which(ratio==ratiosFull)) %>%
  ungroup() %>%
  group_by(OTU, combo, line) %>%
  mutate(ratioOrderBefore=ratioOrder-1,
         ratioOrderAfter=ratioOrder+1,
         before=if_else(ratioOrderBefore==lag(ratioOrder), TRUE, FALSE),
         after=if_else(ratioOrderAfter==lead(ratioOrder), TRUE, FALSE),
         keep=replace(keep, (FALSE %in% after & head(which(after==FALSE))!=length(after)) |
                        (FALSE %in% before & tail(which(before==FALSE))!=1), FALSE)) %>%
  filter(keep==TRUE) %>%
  select(-c(keep, ratioOrder, ratioOrderBefore, ratioOrderAfter, before, after)) %>%
  # Filter out ASVs not detected in at least three adjacent ratios in at least one actual replicate.
  ungroup() %>%
  mutate(keep=FALSE) %>%
  arrange(OTU, factor(combo), factor(line), factor(ratio)) %>%
  group_by(OTU, combo, line) %>%
  mutate(before=ifelse(!is.na(lag(ratio)), TRUE, FALSE),
         after=ifelse(!is.na(lead(ratio)), TRUE, FALSE),
         keep=replace(keep, before==TRUE & after==TRUE, TRUE)) %>%
  ungroup() %>%
  group_by(OTU, combo) %>%
  mutate(ASVcomboKeep=ifelse(keep==TRUE & mixtureType=="actual", TRUE, FALSE)) %>%
  filter(TRUE %in% ASVcomboKeep) %>%
  # Remove ASVs that do not appear in the parent communities.
  mutate(ASVcomboKeep=ifelse(ratio %in% c("1:0", "0:1") & mixtureType=="actual", TRUE, FALSE)) %>%
  filter(TRUE %in% ASVcomboKeep) %>%
  select(-c(keep, before, after, ASVcomboKeep)) %>%
  # Add log abundances, averages and sd after filtering for each combo/ratio.
  ungroup() %>%
  group_by(OTU, label, mixtureType) %>%
  mutate(logAbundance=log10(rel_abundance),
         sdAbundance=sd(logAbundance),
         avgAbundance=mean(logAbundance)) %>%
  # Add ASV counts after filtering.
  ungroup() %>%
  group_by(label, line) %>%
  mutate(ASVcount=length(unique(OTU)))

# Create a filtered dataframe with 10^-4 for the theoretical bottom, without removing values.
# Filter low quality ASVs.
theoretical_mixtures_df_stats <- theoretical_mixtures_df_unfiltered %>%
  # Filter to 10^-4 relative abundance minimum.
  ungroup() %>%
  group_by(label, line) %>%
  mutate(rel_abundance=replace(rel_abundance, rel_abundance <= 10^-4 & mixtureType=="theoretical", 0),
         rel_abundance=replace(rel_abundance, rel_abundance <= 10^-4 & mixtureType=="actual", NA)) %>%
  filter(!is.na(rel_abundance)) %>%
  mutate(total_rel=sum(rel_abundance), 
         rel_abundance=rel_abundance/total_rel,
         rel_abundance=replace(rel_abundance, rel_abundance==0, 10^-4)) %>%
  select(!total_rel) %>%
  # Remove ASV replicates with missing ratios in the middle.
  ungroup() %>%
  mutate(keep=TRUE) %>%
  arrange(OTU, factor(combo), factor(line), factor(ratio)) %>%
  group_by(OTU, combo, line) %>%
  rowwise() %>%
  mutate(ratioOrder=which(ratio==ratiosFull)) %>%
  ungroup() %>%
  group_by(OTU, combo, line) %>%
  mutate(ratioOrderBefore=ratioOrder-1,
         ratioOrderAfter=ratioOrder+1,
         before=if_else(ratioOrderBefore==lag(ratioOrder), TRUE, FALSE),
         after=if_else(ratioOrderAfter==lead(ratioOrder), TRUE, FALSE),
         keep=replace(keep, (FALSE %in% after & head(which(after==FALSE))!=length(after)) |
                        (FALSE %in% before & tail(which(before==FALSE))!=1), FALSE)) %>%
  filter(keep==TRUE) %>%
  select(-c(keep, ratioOrder, ratioOrderBefore, ratioOrderAfter, before, after)) %>%
  # Filter out ASVs not detected in at least three adjacent ratios in at least one actual replicate.
  ungroup() %>%
  mutate(keep=FALSE) %>%
  arrange(OTU, factor(combo), factor(line), factor(ratio)) %>%
  group_by(OTU, combo, line) %>%
  mutate(before=ifelse(!is.na(lag(ratio)), TRUE, FALSE),
         after=ifelse(!is.na(lead(ratio)), TRUE, FALSE),
         keep=replace(keep, before==TRUE & after==TRUE, TRUE)) %>%
  ungroup() %>%
  group_by(OTU, combo) %>%
  mutate(ASVcomboKeep=ifelse(keep==TRUE & mixtureType=="actual", TRUE, FALSE)) %>%
  filter(TRUE %in% ASVcomboKeep) %>%
  # Remove ASVs that do not appear in the parent communities.
  mutate(ASVcomboKeep=ifelse(ratio %in% c("1:0", "0:1") & mixtureType=="actual", TRUE, FALSE)) %>%
  filter(TRUE %in% ASVcomboKeep) %>%
  select(-c(keep, before, after, ASVcomboKeep)) %>%
  # Add log abundances, averages and sd after filtering for each combo/ratio.
  ungroup() %>%
  group_by(OTU, label, mixtureType) %>%
  mutate(logAbundance=log10(rel_abundance),
         sdAbundance=sd(logAbundance),
         avgAbundance=mean(logAbundance)) %>%
  # Add ASV counts after filtering.
  ungroup() %>%
  group_by(label, line) %>%
  mutate(ASVcount=length(unique(OTU)))

# Create dataframe with >0.1% relative abundance.
dfRelAbundant <- df %>%
  filter(rel_abundance>0.001)

# Calculate summary statistics. -----------------------------------------------

# Calculate average theoretical sd across all ratios.
t_sum <- theoretical_mixtures_df_stats %>%
  filter(mixtureType=="theoretical") %>%
  group_by(OTU, combo) %>%
  summarize(Family, ratio, group, label, avgAbundance, sdAbundance, t_sdAvg=mean(sdAbundance, na.rm=TRUE)) %>%
  mutate(t_sdAvg=ifelse(is.nan(t_sdAvg), NA, t_sdAvg)) %>%
  distinct()

# Calculate average actual sd across all ratios.
a_sum <- theoretical_mixtures_df_stats %>%
  filter(mixtureType=="actual") %>%
  group_by(OTU, combo) %>%
  summarize(Family, ratio, group, label, avgAbundance, sdAbundance, a_sdAvg=mean(sdAbundance, na.rm=TRUE)) %>%
  mutate(a_sdAvg=ifelse(is.nan(a_sdAvg), NA, a_sdAvg)) %>%
  distinct()

# Calculate differences between theoretical and actual abundances for each ratio and averages.
summary_df <- full_join(t_sum, a_sum, by=c("OTU", "Family", "combo", "ratio", "group", "label")) %>%
  rename(theoretical_abundance=avgAbundance.x, actual_abundance=avgAbundance.y,
         theoretical_sdAbundance=sdAbundance.x, actual_sdAbundance=sdAbundance.y) %>%
  # Remove parent communities which will always be identical to theoretical.
  filter(group=="mixture") %>%
  mutate(difference=actual_abundance-theoretical_abundance, 
         distance=abs(difference)) %>%
  group_by(OTU, combo) %>%
  mutate(avgdifference=mean(difference, na.rm=TRUE), avgdistance=mean(distance, na.rm=TRUE))

# Pull out specific summary statistics for each ASV/combo.
ASV_sum <- summary_df %>%
  group_by(OTU, Family, combo) %>%
  summarize(t_sdAvg, a_sdAvg, avgdistance, avgdifference) %>%
  # Keeps only non-NA values if there are both NA and non-NA values in sd calculations.
  slice(ifelse(all(is.na(t_sdAvg) | is.na(a_sdAvg)),
               1,
               which(!is.na(t_sdAvg) & !is.na(a_sdAvg)))) %>%
  distinct()

# Remove ASVs with standard deviation >0.25 in either theoretical or actual or only one actual replicate.
ASV_sumFiltered <- ASV_sum %>%
  filter((t_sdAvg<=0.4 | is.na(t_sdAvg)) & (a_sdAvg<=0.25))

theoretical_mixtures_df_statsFiltered <- theoretical_mixtures_df_stats %>%
  left_join(ASV_sum, by=c("OTU", "Family", "combo")) %>%
  filter((t_sdAvg<=0.4 | is.na(t_sdAvg)) & (a_sdAvg<=0.25))