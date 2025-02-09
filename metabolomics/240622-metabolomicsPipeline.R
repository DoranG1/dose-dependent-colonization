library(data.table)
library(tidyverse)
library(cowplot)
library(foreach)
library(viridis)
library(scales)

theme_set(theme_cowplot())
source("../../config/palettes/plotDefaults.R")

# Import data and filter out extra empty columns.
dataNeg <- fread("../../data/KCHU015_CombineSubmit_April2024-neg.txt") %>% 
  filter(INCHIKEY!="Internal Standard") %>% 
  select(-c(MSI, `Alignment ID`, `Average Rt(min)`, `Average Mz`, INCHIKEY, `Adduct type`)) %>% 
  mutate(metaboliteNumber = row_number())
dataPos <- fread("../../data/KCHU015_CombineSubmit_April2024-pos.txt") %>% 
  filter(INCHIKEY!="Internal Standard") %>% 
  select(-c(MSI, `Alignment ID`, `Average Rt(min)`, `Average Mz`, INCHIKEY, `Adduct type`)) %>% 
  mutate(metaboliteNumber = row_number())

# Turn data into long format.
dataNeg <- dataNeg %>% pivot_longer(cols=2:67, names_to="sample", values_to="count") %>% 
  rename(metabolite=`Metabolite name`)
dataPos <- dataPos %>% pivot_longer(cols=2:67, names_to="sample", values_to="count") %>% 
  rename(metabolite=`Metabolite name`)

# Combine positive and negative counts for analysis.
dataCombined <- rbind(
  dataNeg %>% mutate(polarity="neg"),
  dataPos %>% mutate(polarity="pos")
)

# Extract sample information.
dataCombined <- dataCombined %>% 
  mutate(strain = case_when(grepl("BK",sample) ~ "BK",
                            grepl("QC",sample) ~ "QC",
                            !(grepl("BK",sample) | grepl("QC",sample)) ~ sub("_.*","",sample)),
         rep = case_when(grepl("BK",sample) & polarity=="pos" ~ sub("BK_Pos_","",sample),
                         grepl("BK",sample) & polarity=="neg" ~ sub("BK_Neg_","",sample),
                         grepl("QC",sample) & polarity=="pos" ~ sub("QC_Pos_","",sample),
                         grepl("QC",sample) & polarity=="neg" ~ sub("QC_Neg_","",sample),
                         grepl("F1",sample) ~ sub(".*F1_R","",sample)),
         rep = sub("_.*","",rep),
         timepoint = ifelse(grepl("F1",sample), sub(".*R*_T","",sample), NA),
         timepoint = sub("_.*","",timepoint),
         uniqueMetabolite = paste0(metaboliteNumber, "_", polarity)) %>% 
  select(-c(metaboliteNumber, polarity))

# Subtract the average blank values for each metabolite.
dataBlankSubtracted <- left_join(
  dataCombined %>% 
    filter(!(strain %in% c("BK","QC","Lachn00058","Bacte0003"))),
    #filter((timepoint==48 & !(strain %in% c("BK","QC","Lachn0058","Bacte0003"))) |
    #         strain=="mBHI"),
  dataCombined %>% 
    filter(strain=="BK") %>% 
    group_by(uniqueMetabolite) %>% 
    summarize(blankAvg = mean(count))
) %>%
  mutate(blankSubtractedCount = count - blankAvg,
         blankSubtractedCount = replace(blankSubtractedCount, count-blankAvg < 1, 1))

# Plot distribution of blank-subtracted counts.
dataBlankSubtracted %>% 
  ggplot() +
  geom_histogram(aes(x=log10(blankSubtractedCount)), bins=50)

# Annotate metabolites with >10-fold difference between replicate blank-subtracted counts.
dataBlankSubtracted <- dataBlankSubtracted %>% 
  group_by(strain, uniqueMetabolite, timepoint) %>%
  mutate(#replicateRatio = (max(blankSubtractedCount)+1000)/(min(blankSubtractedCount)+1000),
    replicateRatio = (max(blankSubtractedCount)+100)/(min(blankSubtractedCount)+100),
    divergent = replicateRatio > 10)

minDivergedReplicates <- dataBlankSubtracted %>% 
  group_by(strain, timepoint, uniqueMetabolite) %>% 
  filter("1" %in% as.character(blankSubtractedCount)) %>% 
  filter(blankSubtractedCount!=1)

minDivergedReplicates %>% 
  ggplot() +
  geom_violin(aes(x=rep, y=log10(blankSubtractedCount)))

dataBlankSubtracted %>%
  filter(timepoint %in% c(2,4,48)) %>% 
  ggplot() +
  geom_histogram(aes(x=log10(replicateRatio)), bins=50) +
  facet_wrap(~timepoint)

dataBlankSubtracted %>% 
  filter(timepoint %in% c(2,4,48)) %>% 
  ungroup() %>% 
  group_by(strain, uniqueMetabolite, timepoint) %>% 
  # Remove one replicate.
  slice(1) %>% 
  ungroup() %>% 
  group_by(uniqueMetabolite, timepoint) %>% 
  mutate(divergentCount = sum(divergent)) %>% 
  ggplot() +
  geom_histogram(aes(x=divergentCount)) +
  facet_wrap(~timepoint)

dataBlankSubtracted <- dataBlankSubtracted %>% 
  ungroup() %>% 
  # Annotate metabolites with divergent counts across at least 4 strains/communities.
  group_by(uniqueMetabolite, timepoint) %>% 
  mutate(divergentAll = sum(divergent[strain!="mBHI"]) > 8) %>% 
  ungroup() %>% 
  group_by(uniqueMetabolite) %>% 
  mutate(divergentAll = TRUE %in% divergentAll)

# Count number of divergent metabolite-strain pairs and divergentAll metabolites.
sum(dataBlankSubtracted$divergent)
n_distinct(dataBlankSubtracted %>% filter(divergentAll & timepoint==4) %>% pull(uniqueMetabolite))
sum(dataBlankSubtracted %>% filter(strain=="mBHI") %>% pull(divergent))

datamBHI <- dataBlankSubtracted %>%
  filter(strain=="mBHI") %>%
  filter(!divergent & !divergentAll) %>% 
  filter(blankSubtractedCount >= 50000) %>% 
  #filter(blankSubtractedCount >= 20000) %>% 
  group_by(uniqueMetabolite) %>% 
  mutate(mBHICount = mean(blankSubtractedCount))

# Count number of unique metabolites in mBHI.
n_distinct(datamBHI$uniqueMetabolite)

dataConsumed <- left_join(
  datamBHI %>% select(uniqueMetabolite, mBHICount) %>% unique(), 
  dataBlankSubtracted %>% select(!sample) %>% filter(strain!="mBHI")) %>% 
  mutate(foldChange = blankSubtractedCount/mBHICount)

# Write data to table for easier i/o later.
# write.table(dataConsumed, "dataConsumed.txt", quote=FALSE, row.names=FALSE, sep="\t")
dataConsumed <- fread("dataConsumed.txt")

sum(dataConsumed$divergent)
n_distinct(dataConsumed %>% filter(divergentAll) %>% pull(uniqueMetabolite))


samplesRaw <- c("Bacte0126","Tanne0007","EnteC0002","EnteC0003","Strep0007","Strep0017","XFA","XFB")
  
samplesClean <- c("B. fragilis","P. goldsteinii","E. faecalis",
             "E. casseliflavus","L. garvieae","L. lactis",
             "Community D1","Community D2")

# Plot distibution of metabolite significance for each strain.
p_significanceDistributionT48 <- dataConsumed %>%
  filter(timepoint==48 & strain %in% samplesRaw) %>% 
  filter(!divergent) %>% 
  group_by(strain, uniqueMetabolite) %>% 
  summarize(significance = ifelse(max(log10(foldChange))<=-4, "depleted", "none"),
            avgFoldChange = mean(foldChange)) %>%
  ungroup() %>% 
  mutate(strain = case_when(
    strain=="Bacte0126" ~ "B. fragilis",
    strain=="Tanne0007" ~ "P. goldsteinii",
    strain=="EnteC0002" ~ "E. faecalis",
    strain=="EnteC0003" ~ "E. casseliflavus",
    strain=="Strep0007" ~ "L. garvieae",
    strain=="Strep0017" ~ "L. lactis",
    strain=="XFA" ~ "Community D1",
    strain=="XFB" ~ "Community D2"),
    strain = fct_relevel(strain, c("E. faecalis", "E. casseliflavus", "L. garvieae", "L. lactis",
                                   "B. fragilis", "P. goldsteinii", "Community D1", "Community D2"))) %>% 
  ggplot() +
  geom_jitter(aes(x=strain, y=log10(avgFoldChange), color=significance), width=0.25, size=0.5, alpha=0.25) +
  geom_violin(aes(x=strain, y=log10(avgFoldChange)), size=0.5) +
  xlab("Strain") +
  ylab("Feature fold change from media") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  theme(legend.key.size=unit(0.5,"lines")) +
  scale_y_continuous(limits=c(-10,6), breaks=c(-10,-8,-6,-4,-2,0,2,4,6)) +
  theme(axis.ticks.x = element_line(size=0.25),
        axis.line.x = element_line(size=0.25),
        axis.ticks.y = element_line(size=0.25),
        axis.line.y = element_line(size=0.25)) +
  theme(legend.position="none") +
  scale_color_manual(name = "Depletion threshold",
                     values = c("#8EA604","#C2C2C2")) +
  DEFAULTS.THEME_PRINT
p_significanceDistributionT48
save_plot("out/significanceDistributionT48.tiff", p_significanceDistributionT48, 
          base_width=2.2, base_height=2.3, dpi=400)

# Plot scatterplot of fold-changes across both replicates for each strain.
replicateFoldChangeComparisonCor <- foreach(x = samplesRaw, .combine="rbind") %do% {
  currDf <- dataConsumed %>% 
    filter(timepoint==48 & strain==x) %>% 
    select(uniqueMetabolite, strain, rep, divergent, foldChange) %>% 
    pivot_wider(names_from=rep, values_from=foldChange) %>%
    rename(rep1 = `1`, rep2 = `2`) %>% 
    mutate(rep1 = log10(rep1), rep2 = log10(rep2)) %>% 
    mutate(coloring = case_when(divergent ~ "Divergent",
                                !divergent & rep1 < -4 & rep2 < -4 ~ ">10000-fold depleted",
                                .default = "None"),
           strain = case_when(strain=="Bacte0126" ~ "B. fragilis",
                              strain=="Tanne0007" ~ "P. goldsteinii",
                              strain=="EnteC0002" ~ "E. faecalis",
                              strain=="EnteC0003" ~ "E. casseliflavus",
                              strain=="Strep0007" ~ "L. garvieae",
                              strain=="Strep0017" ~ "L. lactis",
                              strain=="XFA" ~ "Community D1",
                              strain=="XFB" ~ "Community D2",))
  cor <- cor.test(currDf$rep1, currDf$rep2, method="pearson")
  p_replicateFoldChangeScatter <- currDf %>% 
    ggplot() + 
    geom_point(aes(x=rep1, y=rep2, color=coloring), size=0.5, alpha=0.75) +
    geom_abline(aes(slope=1, intercept=0), linetype="dashed", size=0.5) +
    scale_x_continuous(limits=c(-10,6), breaks=c(-10,-8,-6,-4,-2,0,2,4,6)) +
    scale_y_continuous(limits=c(-10,6), breaks=c(-10,-8,-6,-4,-2,0,2,4,6)) +
    xlab("Feature fold-change\nfrom media, replicate 1") +
    ylab("Feature fold-change\nfrom media, replicate 2") +
    scale_color_manual(name = "Depletion threshold",
                       values = c("#8EA604","#8F6593","#C2C2C2")) +
    theme(legend.key.size=unit(0.5,"lines")) +
    theme(legend.position="none") +
    facet_wrap(~strain) +
    DEFAULTS.THEME_PRINT
  p_replicateFoldChangeScatter
  save_plot(paste0("out/replicateFoldChangeCorrelations/",x,".png"), p_replicateFoldChangeScatter)
  save_plot(paste0("out/replicateFoldChangeCorrelations/",x,".pdf"), p_replicateFoldChangeScatter,
            base_width=1.75, base_height=1.75)
  corDf <- data.frame(strain=x, r = cor$estimate, p = cor$p.value)
}

foldChangeScatterLegend <- get_legend(p_replicateFoldChangeScatter)
save_plot("out/foldChangeScatterLegend.pdf", foldChangeScatterLegend, base_width=1, base_height=0.5)

# Extract correlation for just E. faecalis, with and without divergent points.
faecalisCorDf <- dataConsumed %>% 
  filter(timepoint==48 & strain=="EnteC0002") %>% 
  select(uniqueMetabolite, rep, divergent, foldChange) %>% 
  pivot_wider(names_from=rep, values_from=foldChange) %>%
  rename(rep1 = `1`, rep2 = `2`) %>% 
  mutate(rep1 = log10(rep1), rep2 = log10(rep2)) %>% 
  mutate(coloring = case_when(divergent ~ "Divergent",
                              !divergent & rep1 < -4 & rep2 < -4 ~ ">10000-fold depleted",
                              .default = "None"))
faecalisCorDivergent <- cor.test(faecalisCorDf$rep1, faecalisCorDf$rep2, method="pearson")
faecalisCorNonDivergent <- cor.test((faecalisCorDf %>% filter(coloring!="Divergent"))$rep1, 
                                 (faecalisCorDf %>% filter(coloring!="Divergent"))$rep2, method="pearson")

# Calculate overlaps between different strains with a -10000-fold cutoff.

# For each pair of strains, get metabolites in each category of shared or unique depletion.
strainOverlaps <- foreach(x=samplesRaw, .combine="rbind") %do% {
  foreach(y=samplesRaw, .combine="rbind") %do% {
    if (x!=y) {
      # Get list of metabolite-polarity pairs that are divergent in either strain to remove from analysis.
      divergentX <- dataConsumed %>%
        filter(timepoint==48) %>% 
        filter(strain==x & divergent) %>% 
        pull(uniqueMetabolite) %>% 
        unique()
      divergentY <- dataConsumed %>%
        filter(timepoint==48) %>% 
        filter(strain==y & divergent) %>% 
        pull(uniqueMetabolite) %>% 
        unique()
      # Get consumed metabolites.
      overlapData <- dataConsumed %>%
        filter(timepoint==48) %>% 
        filter(strain %in% c(x,y)) %>%
        filter(!(uniqueMetabolite %in% divergentX) & !(uniqueMetabolite %in% divergentY)) %>% 
        group_by(strain, uniqueMetabolite) %>%
        mutate(unifiedSignificance = case_when(strain==x & max(log10(foldChange))<=-4 ~ "depleted",
                                               strain==y & max(log10(foldChange))<=-4 ~ "depleted",
                                               .default = "none")) %>%
        ungroup() %>% 
        filter(unifiedSignificance=="depleted") %>%
        select(strain, uniqueMetabolite, unifiedSignificance) %>%
        unique() %>% 
        mutate(comparison = paste0(x, "-", y),
               strain = ifelse(strain==x, "strain1", "strain2")) %>% 
        pivot_wider(names_from=strain, values_from=unifiedSignificance)
      print(paste0(x,"-",y))
      sharedConsumed <- overlapData %>% filter(strain1=="depleted" & strain2=="depleted")
      consumedX <- overlapData %>% filter(strain1=="depleted" & is.na(strain2))
      consumedY <- overlapData %>% filter(is.na(strain1) & strain2=="depleted")
      intersection <- nrow(sharedConsumed)
      union <- nrow(consumedX) + nrow(consumedY) + intersection
      consumptionData <- rbind(sharedConsumed, consumedX, consumedY) %>% 
        mutate(numConsumed1 = (nrow(consumedX) + intersection), numConsumed2 = (nrow(consumedY) + intersection), 
               overlap=intersection, overlapJaccard=round(intersection/union, 6),
               overlapAsymmetric = round(intersection / (nrow(consumedX) + intersection), 6))
    }
  }
}
# write.table(strainOverlaps, "strainOverlaps.txt", row.names=FALSE, quote=FALSE, sep="\t")
strainOverlaps <- fread("strainOverlaps.txt")

# Calculate overlaps between two replicates of each strain with different thresholds.
# Set lower significance to -2 for same cutoff between replicates.
replicateOverlaps <- foreach(x=samplesRaw, .combine="rbind") %do% {
  overlapData <- dataConsumed %>%
    filter(timepoint==48) %>% 
    filter(!divergent) %>% 
    filter(strain==x) %>%
    group_by(uniqueMetabolite) %>%
    mutate(unifiedSignificance1 = ifelse(rep==1 & log10(foldChange) <= -4, 
                                         "depleted", "none"),
           lowerSignificance1 = ifelse(rep==1 & log10(foldChange) <= -4, 
                                       "depleted", "none"),
           unifiedSignificance1 = ifelse("depleted" %in% unifiedSignificance1, "depleted", "none"),
           lowerSignificance1 = ifelse("depleted" %in% lowerSignificance1, "depleted", "none"),
           unifiedSignificance2 = ifelse(rep==2 & log10(foldChange) <= -4, 
                                         "depleted", "none"),
           lowerSignificance2 = ifelse(rep==2 & log10(foldChange) <= -4, "depleted", "none"),
           unifiedSignificance2 = ifelse("depleted" %in% unifiedSignificance2, "depleted", "none"),
           lowerSignificance2 = ifelse("depleted" %in% lowerSignificance2, "depleted", "none")) %>% 
    ungroup() %>%
    filter(lowerSignificance1=="depleted" | lowerSignificance2=="depleted") %>% 
    select(strain, uniqueMetabolite, 
           unifiedSignificance1, unifiedSignificance2, lowerSignificance1, lowerSignificance2) %>% 
    unique() 
  significantConsumed1 <- which((overlapData %>% pull(unifiedSignificance1))=="depleted")
  significantConsumed2 <- which((overlapData %>% pull(unifiedSignificance2))=="depleted")
  lowerConsumed1 <- which((overlapData %>% pull(lowerSignificance1))=="depleted")
  lowerConsumed2 <- which((overlapData %>% pull(lowerSignificance2))=="depleted")
  intersection1 <- length(intersect(significantConsumed1, lowerConsumed2))
  intersection2 <- length(intersect(lowerConsumed1, significantConsumed2))
  union1 <- length(significantConsumed1) + length(lowerConsumed2)
  union2 <- length(lowerConsumed1) + length(significantConsumed2)
  df <- data.frame(strain=x, 
                   depleted1 = length(significantConsumed1), depleted2 = length(significantConsumed2),
                   overlap = (intersection1+intersection2)/2, 
                   overlapJaccard = ((intersection1/union1)+(intersection2/union2))/2, 
                   overlapAsymmetric = (
                     ((intersection1/length(significantConsumed1))+(intersection2/length(significantConsumed2)))/2))
}
# write.table(replicateOverlaps, "replicateOverlaps.txt", row.names=FALSE, quote=FALSE, sep="\t")
replicateOverlaps <- fread("replicateOverlaps.txt")

# Plot heatmap of niche overlaps between strains.
p_strainOverlaps <- strainOverlaps %>%
  mutate(strain1 = sub("-.*","",comparison), strain2 = sub(".*-","",comparison)) %>% 
  select(strain1, strain2, overlapAsymmetric) %>% 
  unique() %>% 
  filter(strain1!=strain2) %>%
  #filter(!(strain1 %in% c("XFA","XFB")) & !(strain2 %in% c("XFA","XFB"))) %>% 
  # Add avg overlap between replicates on diagonal.
  rbind(replicateOverlaps %>% 
          #filter(!(strain %in% c("XFA","XFB"))) %>% 
          select(strain, overlapAsymmetric) %>% 
          rename(strain1 = strain) %>% 
          mutate(strain2 = strain1)) %>% 
  ungroup() %>% 
  mutate(strain1 = fct_relevel(strain1, samplesClean),
         strain2 = fct_relevel(strain2, rev(samplesClean)),
         overlapAsymmetricRounded = round(overlapAsymmetric,2)) %>% 
  mutate(strain1 = case_when(strain1=="Bacte0126" ~ "B. fragilis",
                             strain1=="Tanne0007" ~ "P. goldsteinii",
                             strain1=="EnteC0002" ~ "E. faecalis",
                             strain1=="EnteC0003" ~ "E. casseliflavus",
                             strain1=="Strep0007" ~ "L. garvieae",
                             strain1=="Strep0017" ~ "L. lactis",
                             strain1=="XFA" ~ "Community D1",
                             strain1=="XFB" ~ "Community D2"),
         strain2 = case_when(strain2=="Bacte0126" ~ "B. fragilis",
                             strain2=="Tanne0007" ~ "P. goldsteinii",
                             strain2=="EnteC0002" ~ "E. faecalis",
                             strain2=="EnteC0003" ~ "E. casseliflavus",
                             strain2=="Strep0007" ~ "L. garvieae",
                             strain2=="Strep0017" ~ "L. lactis",
                             strain2=="XFA" ~ "Community D1",
                             strain2=="XFB" ~ "Community D2"),
         strain1 = fct_relevel(strain1, rev(c("E. faecalis","E. casseliflavus","L. garvieae","L. lactis",
                                              "B. fragilis","P. goldsteinii","Community D1","Community D2"))),
         strain2 = fct_relevel(strain2, c("E. faecalis","E. casseliflavus","L. garvieae","L. lactis",
                                          "B. fragilis","P. goldsteinii","Community D1","Community D2")),
         textColor = ifelse(overlapAsymmetricRounded>0.5, TRUE, FALSE)) %>% 
  ggplot() +
  geom_tile(aes(x=strain2, y=strain1, fill=overlapAsymmetric), color="grey30", size=0.01) +
  geom_text(aes(x=strain2, y=strain1, label=overlapAsymmetricRounded, color=textColor), size=2) +
  scale_color_manual(values=c("white","black"), guide="none") +
  # xlab("\n\nFocal strain") +
  # ylab("Comparison strain\n\n") +
  xlab("\n\nComparison strain") +
  ylab("Focal strain\n\n") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  theme(legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  scale_fill_viridis(name = "Niche overlap", option="viridis",
                     limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1)) +
  DEFAULTS.THEME_PRINT
p_strainOverlaps
save_plot("out/strainOverlap.png", p_strainOverlaps)
save_plot("out/strainOverlap.pdf", p_strainOverlaps, base_width=2.9, base_height=2.9)


# Analyze resource consumption rates for each group of resources in each comparison.
avgConsumptionRatesT4 <- foreach(x=samplesRaw, .combine="rbind") %do% {
  foreach(y=samplesRaw, .combine="rbind") %do% {
    if (x!=y) {
      sharedConsumed <- strainOverlaps %>% 
        filter(comparison==paste0(x, "-", y) & strain1=="depleted" & strain2=="depleted") %>% 
        pull(uniqueMetabolite)
      consumedX <- strainOverlaps %>% 
        filter(comparison==paste0(x, "-", y) & strain1=="depleted" & is.na(strain2)) %>% 
        pull(uniqueMetabolite)
      consumedY <- strainOverlaps %>% 
        filter(comparison==paste0(x, "-", y) & is.na(strain1) & strain2=="depleted") %>% 
        pull(uniqueMetabolite)
      avgFoldChangeSharedXrep1 <- median(log10(dataConsumed %>%
                                                 filter(strain==x & !divergent & uniqueMetabolite %in% sharedConsumed & 
                                                          timepoint==4 & rep==1) %>% 
                                           pull(foldChange)))
      avgFoldChangeSharedXrep2 <- median(log10(dataConsumed %>%
                                                 filter(strain==x & !divergent & uniqueMetabolite %in% sharedConsumed & 
                                                          timepoint==4 & rep==2) %>% 
                                                 pull(foldChange)))
      avgFoldChangeSharedYrep1 <- median(log10(dataConsumed %>%
                                                 filter(strain==y & !divergent & uniqueMetabolite %in% sharedConsumed & 
                                                          timepoint==4 & rep==1) %>% 
                                                 pull(foldChange)))
      avgFoldChangeSharedYrep2 <- median(log10(dataConsumed %>%
                                                 filter(strain==y & !divergent & uniqueMetabolite %in% sharedConsumed & 
                                                          timepoint==4 & rep==2) %>% 
                                                 pull(foldChange)))
      avgFoldChangeUniqueXrep1 <- median(log10(dataConsumed %>%
                                                 filter(strain==x & !divergent & uniqueMetabolite %in% consumedX & 
                                                          timepoint==4 & rep==1) %>% 
                                                 pull(foldChange)))
      avgFoldChangeUniqueXrep2 <- median(log10(dataConsumed %>%
                                                 filter(strain==x & !divergent & uniqueMetabolite %in% consumedX & 
                                                          timepoint==4 & rep==2) %>% 
                                                 pull(foldChange)))
      avgFoldChangeUniqueYrep1 <- median(log10(dataConsumed %>%
                                                 filter(strain==y & !divergent & uniqueMetabolite %in% consumedY & 
                                                          timepoint==4 & rep==1) %>% 
                                                 pull(foldChange)))
      avgFoldChangeUniqueYrep2 <- median(log10(dataConsumed %>%
                                                 filter(strain==y & !divergent & uniqueMetabolite %in% consumedY & 
                                                          timepoint==4 & rep==2) %>% 
                                                 pull(foldChange)))
      print(paste0(x,"-",y))
      df <- data.frame(comparison = paste0(x, "-", y), 
                       rateUnique1rep1 = avgFoldChangeUniqueXrep1,
                       rateUnique1rep2 = avgFoldChangeUniqueXrep2,
                       rateUnique2rep1 = avgFoldChangeUniqueYrep1,
                       rateUnique2rep2 = avgFoldChangeUniqueYrep2,
                       rateShared1rep1 = avgFoldChangeSharedXrep1,
                       rateShared1rep2 = avgFoldChangeSharedXrep2,
                       rateShared2rep1 = avgFoldChangeSharedYrep1,
                       rateShared2rep2 = avgFoldChangeSharedYrep2)
    }
  }
}

overlapConsumptionRatesDfT4 <- left_join(strainOverlaps %>% select(comparison, overlapAsymmetric),
                                       avgConsumptionRatesT4) %>% 
  unique()
# write.table(overlapConsumptionRatesDfT4, "overlapConsumptionRatesT4.txt", quote=FALSE, row.names=FALSE, sep="\t")
overlapConsumptionRatesDf <- fread("overlapConsumptionRatesT4.txt")

overlapConsumptionRatesDfCleaned <- overlapConsumptionRatesDf %>% 
  mutate(strain1 = sub("-.*","",comparison), 
         strain2 = sub(".*-","",comparison)) %>% 
  filter(!(strain1 %in% c("XFA","XFB")) & !(strain2 %in% c("XFA","XFB"))) %>% 
  mutate(strain1 = case_when(strain1=="Bacte0126" ~ "B. fragilis",
                             strain1=="Tanne0007" ~ "P. goldsteinii",
                             strain1=="EnteC0002" ~ "E. faecalis",
                             strain1=="EnteC0003" ~ "E. casseliflavus",
                             strain1=="Strep0007" ~ "L. garvieae",
                             strain1=="Strep0017" ~ "L. lactis"),
         strain2 = case_when(strain2=="Bacte0126" ~ "B. fragilis",
                             strain2=="Tanne0007" ~ "P. goldsteinii",
                             strain2=="EnteC0002" ~ "E. faecalis",
                             strain2=="EnteC0003" ~ "E. casseliflavus",
                             strain2=="Strep0007" ~ "L. garvieae",
                             strain2=="Strep0017" ~ "L. lactis"),
         strain1 = fct_relevel(strain1, c("E. faecalis", "E. casseliflavus", "L. garvieae", "L. lactis",
                                          "B. fragilis", "P. goldsteinii")),
         strain2 = fct_relevel(strain2, c("E. faecalis", "E. casseliflavus", "L. garvieae", "L. lactis",
                                          "B. fragilis", "P. goldsteinii")),
         comparison = paste0(strain1, "/", strain2),
         R1a = round(abs((rateUnique1rep1 + rateUnique1rep2)/2), 2),
         R1ab = round(abs((rateShared1rep1 + rateShared1rep2)/2), 2),
         R2ab = round(abs((rateShared2rep1 + rateShared2rep2)/2), 2),
         R2b = round(abs((rateUnique2rep1 + rateUnique2rep2)/2), 2)) %>%
  filter(strain1!=strain2 & strain1!="P. goldsteinii" &
           (strain1=="E. faecalis" |
              strain1=="E. casseliflavus" & strain2!="E. faecalis" |
              strain1=="L. garvieae" & !(strain2 %in% c("E. faecalis", "E. casseliflavus")) |
              strain1=="L. lactis" & !(strain2 %in% c("E. faecalis", "E. casseliflavus", "L. garvieae")) |
              strain1=="B. fragilis" & strain2=="P. goldsteinii")) %>%
  select(comparison, R1a, R1ab, R2ab, R2b) %>% 
  unique()
write.table(overlapConsumptionRatesDfCleaned, "overlapConsumptionRatesCleaned.txt", quote=FALSE, row.names=FALSE, sep="\t")

# Plot heatmap of consumption rates.
p_consumptionRateHeatmapT4 <- overlapConsumptionRatesDf %>% 
  pivot_longer(cols = c("rateUnique1rep1","rateUnique1rep2","rateUnique2rep1","rateUnique2rep2",
                        "rateShared1rep1","rateShared1rep2","rateShared2rep1","rateShared2rep2"), 
               names_to = "resourceRate", values_to = "rate") %>% 
  filter(comparison %in% c("EnteC0002-Strep0017","EnteC0003-Strep0017","Strep0007-Strep0017",
                           "Bacte0126-Strep0017","Bacte0126-Tanne0007")) %>% 
  mutate(comparison = case_when(comparison=="Bacte0126-Tanne0007" ~ "B. fragilis / P. goldsteinii",
                                comparison=="Bacte0126-Strep0017" ~ "B. fragilis / L. lactis",
                                comparison=="Strep0007-Strep0017" ~ "L. garvieae / L. lactis",
                                comparison=="EnteC0003-Strep0017" ~ "E. casseliflavus / L. lactis",
                                comparison=="EnteC0002-Strep0017" ~ "E. faecalis / L. lactis"),
         comparison = fct_relevel(comparison, c("E. faecalis / L. lactis", "L. garvieae / L. lactis",
                                                "E. casseliflavus / L. lactis","B. fragilis / P. goldsteinii",
                                                "B. fragilis / L. lactis")),
         rep = ifelse(grepl("rep1", resourceRate), "Replicate 1", "Replicate 2"),
         resourceRate = sub("rep.*","",resourceRate)) %>% 
  group_by(comparison, resourceRate) %>% 
  summarize(avgRate = -mean(rate)) %>% 
  mutate(textColor = ifelse(avgRate>0.5, TRUE, FALSE),
         resourceRate = case_when(resourceRate=="rateUnique1" ~ "R1a",
                                  resourceRate=="rateUnique2" ~ "R2b",
                                  resourceRate=="rateShared1" ~ "R1ab",
                                  resourceRate=="rateShared2" ~ "R2ab"),
         resourceRate = fct_relevel(resourceRate, c("R1a","R1ab","R2ab","R2b"))) %>% 
  ggplot() +
  geom_tile(aes(x=resourceRate, y=comparison, fill=avgRate), color="gray30", size=0.01) +
  geom_text(aes(x=resourceRate, y=comparison, label=round(avgRate, 2), color=textColor), size=2) +
  scale_fill_viridis(name = "Consumption rate", option="inferno",
                     limits = c(0,0.75), breaks=c(0,0.25,0.5,0.75),
                     labels = c(0,0.25,0.5,0.75)) +
  scale_color_manual(values=c("white","black"), guide="none") +
  xlab("Consumption rate") +
  ylab("Pairwise mixture") +
  #theme(legend.position="bottom") +
  theme(legend.position="none") +
  theme(legend.key.size=unit(0.5,"lines")) +
  #facet_grid(~resourceGrouping, scales="free") +
  DEFAULTS.THEME_PRINT +
  theme(strip.text.x=element_blank(),
        strip.background=element_blank())
p_consumptionRateHeatmapT4
save_plot("out/consumptionRateHeatmapT4.png", p_consumptionRateHeatmapT4)
save_plot("out/consumptionRateHeatmapT4.pdf", p_consumptionRateHeatmapT4, base_width=2.3, base_height=1.7)

consumptionRateLegend <- get_legend(p_consumptionRateHeatmapT4)
save_plot("out/consumptionRateLegend.pdf",consumptionRateLegend, base_width=2, base_height=0.5)

# Plot distribution of T4 fold changes as violin plots.
EnteC2Strep17RatesAll <- rbind(
  dataConsumed %>% 
    filter(strain=="EnteC0002" & timepoint==4 & !divergent & 
             uniqueMetabolite %in% (strainOverlaps %>% 
                                      filter(comparison=="EnteC0002-Strep0017" & 
                                               strain1=="depleted" & is.na(strain2)) %>% 
                                      pull(uniqueMetabolite))) %>% 
    mutate(resource = "a"),
  dataConsumed %>% 
    filter(strain=="EnteC0002" & timepoint==4 & !divergent &
             uniqueMetabolite %in% (strainOverlaps %>% 
                                      filter(comparison=="EnteC0002-Strep0017" & 
                                               strain1=="depleted" & strain2=="depleted") %>% 
                                      pull(uniqueMetabolite))) %>% 
    mutate(resource = "1ab1"),
  dataConsumed %>% 
    filter(strain=="Strep0017" & timepoint==4 & !divergent &
             uniqueMetabolite %in% (strainOverlaps %>% 
                                      filter(comparison=="EnteC0002-Strep0017" & 
                                               strain1=="depleted" & strain2=="depleted") %>% 
                                      pull(uniqueMetabolite))) %>% 
    mutate(resource = "2ab2"),
  dataConsumed %>% 
    filter(strain=="Strep0017" & timepoint==4 & !divergent &
             uniqueMetabolite %in% (strainOverlaps %>% 
                                      filter(comparison=="EnteC0002-Strep0017" & 
                                               is.na(strain1) & strain2=="depleted") %>% 
                                      pull(uniqueMetabolite))) %>% 
    mutate(resource = "b")) %>% 
  group_by(uniqueMetabolite, strain, resource) %>% 
  summarize(avgFoldChange = mean(foldChange))
  
p_EnteC2Strep17Rates <- EnteC2Strep17RatesAll %>% 
  ungroup() %>% 
  mutate(resource = fct_relevel(resource, c("a","1ab1","2ab2","b"))) %>% 
  ggplot() +
  geom_violin(aes(x=resource, y=log10(avgFoldChange), fill=strain), size=0.25) +
  xlab("Resource") +
  scale_y_continuous(name="Fold change from media at 4 hr",
                     limits=c(-6.8,1.25), breaks=c(-6,-5,-4,-3,-2,-1,0,1),
                     labels=label_math(10^.x)) +
  scale_fill_manual(values = c("#C1A92F","#343795")) +
  theme(legend.title = element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_EnteC2Strep17Rates
save_plot("out/EnteC2Strep17Rates.pdf", p_EnteC2Strep17Rates, base_width=2, base_height=1.67)

# Check statistical significance between different resource consumption rates for EnteC2 and Strep17 resources.
EnteC2Strep17ratesTestShared <- wilcox.test(EnteC2Strep17RatesAll %>% filter(resource=="1ab1") %>% pull(avgFoldChange),
                                              EnteC2Strep17RatesAll %>% filter(resource=="2ab2") %>% pull(avgFoldChange),
                                            alternative="less")
EnteC2Strep17ratesTestEnteC <- wilcox.test(EnteC2Strep17RatesAll %>% filter(resource=="a") %>% pull(avgFoldChange),
                                            EnteC2Strep17RatesAll %>% filter(resource=="1ab1") %>% pull(avgFoldChange),
                                           alternative="greater")
EnteC2Strep17ratesTestStrep17 <- wilcox.test(EnteC2Strep17RatesAll %>% filter(resource=="2ab2") %>% pull(avgFoldChange),
                                           EnteC2Strep17RatesAll %>% filter(resource=="b") %>% pull(avgFoldChange),
                                           alternative="less")
EnteC2Strep17ratesTestUnique <- wilcox.test(EnteC2Strep17RatesAll %>% filter(resource=="a") %>% pull(avgFoldChange),
                                             EnteC2Strep17RatesAll %>% filter(resource=="b") %>% pull(avgFoldChange),
                                             alternative="less")
