library(tidyverse)
library(phyloseq)
library(cowplot)
library(data.table)
library(foreach)
library(ape)
#library(ggtree)
library(treeio)
#library(tidytree)
library(ggdendro)
library(phylogram)

theme_set(theme_cowplot())
source("../../config/palettes/plotDefaults.R")

data <- fread("../mixtureDataframe.txt")

# Import family 5-letter codes.
familyCoding <- fread("../../config/familyCodingFull.txt")

# Import taxa color palette.
palette <- read.table("../../config/palettes/KCHcolors-Silva-partial.txt",
                      header=TRUE, stringsAsFactors = FALSE) %>% 
  mutate(taxashort=gsub("*.*\\.","",taxa))

# Bind taxa colors by family.
data <- data %>% left_join(palette %>% select(taxashort, taxa, hex) %>%
                         filter(taxashort %in% data$Family) %>% 
                         rename(Family=taxashort) %>% 
                         group_by(Family) %>% 
                         mutate(hexnum=row_number()) %>% 
                         filter(hexnum==1) %>% 
                         select(!hexnum),
                       by="Family")

# Create stripped-down color palette containing only families present in dataset.
palette <- palette %>% 
  filter(taxashort %in% sort(unique(data %>% filter(rel_abundance > 0.01) %>% pull(Family))))
# Make a named list.
paletteVector <- palette$hex
names(paletteVector) <- palette$taxashort

# Annotate number of distinct patterns per mixture.
dataPatterns <- data %>% 
  filter(passage==5 & category!="lowAbundance") %>% 
  group_by(OTU) %>% 
  summarize(Family, OTUnum, combos = n_distinct(combo), patterns = n_distinct(category)) %>% 
  unique()

# Get ASVs present in multiple mixtures with consistent patterns.
singlePatternASVs <- dataPatterns %>%
  filter(combos > 1 & patterns==1) %>% 
  select(Family, OTUnum, combos) %>% 
  left_join(data %>% 
              filter(passage==5 & category!="lowAbundance") %>% 
              select(Family, OTUnum, category)) %>% 
  unique()

# Get distribution of other colonization patterns for ASVs in each category.
categories <- c("noisy", "flat", "overcolonizing", "undercolonizing", "dose-dependent")
categoriesRBY <- c("overcolonizing","undercolonizing","dose-dependent")

numASVPatterns <- foreach(x=categories, .combine="rbind") %do% {
  # Get set of ASVs belonging to the current category.
  ASVs <- unique((data %>% filter(passage==5 & category==x))$OTU)
  numASVPatterns <- data %>% 
    filter(passage==5 & category!="lowAbundance" & OTU %in% ASVs) %>%
    filter(category %in% categoriesRBY) %>% 
    ungroup() %>% 
    select(Family, OTUnum, combo, category) %>% 
    unique() %>% 
    group_by(Family, OTUnum) %>% 
    # Count number of combos each ASV is present in.
    mutate(nCombos = n_distinct(combo)) %>% 
    ungroup() %>% 
    group_by(Family, OTUnum, category) %>%
    # Count number of combos per category.
    summarize(nCombos, nOTUs = n()) %>% 
    ungroup() %>% 
    unique() %>% 
    # Bind family coding.
    left_join(familyCoding %>% rename(Family=family), by="Family") %>% 
    # Count number of combos with current category.
    mutate(nSameSet = ifelse(category==x, nOTUs, 0),
           OTUID = paste0(code, "-", OTUnum),
           set = x) %>% 
    group_by(code, OTUnum) %>% 
    mutate(nSameSet = sum(nSameSet)) %>% 
    ungroup() %>% 
    mutate(category = fct_relevel(category, categories))
}

combos <- c("XBA-XBB", "XCA-XCB", "XDA-XDB", "XFA-XFB", "XBA-XCA", "XCA-XDA", "XDA-XFA", "XFA-XBA")

# Get simplified version of dataframe with just combos and colonization patterns.
dataColonization <- data %>% 
  filter(passage==5 & category!="lowAbundance") %>% 
  filter(category %in% c("dose-dependent", "overcolonizing", "undercolonizing")) %>% 
  ungroup() %>% 
  select(Family, OTUnum, combo, category) %>% 
  unique() %>% 
  left_join(familyCoding %>% rename(Family=family), by="Family") %>% 
  mutate(OTUID = paste0(Family, "-", OTUnum),
         OTUIDshort = paste0(code, "-", OTUnum)) %>% 
  group_by(OTUID) %>% 
  mutate(nCombos = n_distinct(combo)) %>% 
  ungroup() %>% 
  mutate(OTUID = fct_reorder(OTUID, nCombos),
         combo = fct_relevel(combo, combos))

# Plots. -----------------------------------------------------------------

# Plot number of distinct patterns per number of mixtures.
p_patternsPerMixtureNum <- dataPatterns %>% 
  ggplot() +
  geom_jitter(aes(x=combos, y=patterns, color=taxa), width=0.1, height=0.1) +
  geom_smooth(aes(x=combos, y=patterns)) +
  scale_color_manual(values=taxaPaletteList) +
  xlim(0.85,8.15) +
  ylim(0.85,5.15) +
  theme(legend.position="none") + 
  DEFAULTS.THEME_PRES
p_patternsPerMixtureNum
save_plot("out/patternsPerMixtureNum.png", p_patternsPerMixtureNum)

p_patternsPerMixtureNumPerFamily <- dataPatterns %>% 
  ggplot() +
  geom_jitter(aes(x=combos, y=patterns, color=taxa), width=0.1, height=0.1) +
  scale_color_manual(values=taxaPaletteList) +
  xlim(0.85,8.15) +
  ylim(0.85,5.15) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~Family, nrow=3)
p_patternsPerMixtureNumPerFamily
save_plot("out/patternsPerMixtureNumPerFamily.png", p_patternsPerMixtureNumPerFamily, nrow=2, ncol=3)

# Plot number of distinct patterns per family.
p_patternsPerFamily <- dataPatterns %>% 
  ggplot() +
  geom_jitter(aes(x=Family, y=patterns, color=taxa)) +
  scale_color_manual(values=taxaPaletteList) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES
p_patternsPerFamily
save_plot("out/patternsPerFamily.png", p_patternsPerFamily, nrow=1, ncol=1.5)

# Plot distribution of categories of DD ASVs.
p_numDDASVPatterns <- numASVPatterns %>% 
  filter(set=="dose-dependent") %>% 
  mutate(OTUID = fct_reorder(OTUID, -nSameSet)) %>% 
  ggplot() +
  geom_bar(aes(x=OTUID, y=nOTUs, fill=category), stat="identity") +
  scale_fill_manual(name="Category",
                    labels=c("noisy"="Noisy",
                             "flat"="Flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("noisy"="#999999",
                             "flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#3C91E6",
                             "dose-dependent"="#A63A43")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  DEFAULTS.THEME_PRES
p_numDDASVPatterns
save_plot("out/numDDASVPatternsRBY.png", p_numDDASVPatterns, nrow=1.25, ncol=1.5)

# Plot distribution of categories of OC ASVs.
p_numOCASVPatterns <- numASVPatterns %>% 
  filter(set=="overcolonizing") %>% 
  mutate(OTUID = fct_reorder(OTUID, -nSameSet)) %>% 
  ggplot() +
  geom_bar(aes(x=OTUID, y=nOTUs, fill=category), stat="identity") +
  scale_fill_manual(name="Category",
                    labels=c("noisy"="Noisy",
                             "flat"="Flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("noisy"="#999999",
                             "flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#3C91E6",
                             "dose-dependent"="#A63A43")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  DEFAULTS.THEME_PRES
p_numOCASVPatterns
save_plot("out/numOCASVPatternsRBY.png", p_numOCASVPatterns, nrow=1.25, ncol=1.5)

# Plot distribution of categories of UC ASVs.
p_numUCASVPatterns <- numASVPatterns %>% 
  filter(set=="undercolonizing") %>% 
  mutate(OTUID = fct_reorder(OTUID, -nSameSet)) %>% 
  ggplot() +
  geom_bar(aes(x=OTUID, y=nOTUs, fill=category), stat="identity") +
  scale_fill_manual(name="Category",
                    labels=c("noisy"="Noisy",
                             "flat"="Flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("noisy"="#999999",
                             "flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#3C91E6",
                             "dose-dependent"="#A63A43")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  DEFAULTS.THEME_PRES
p_numUCASVPatterns
save_plot("out/numUCASVPatternsRBY.png", p_numUCASVPatterns, nrow=1.25, ncol=1.5)

p_numAllPatterns <- numASVPatterns %>% 
  filter(set %in% c("dose-dependent","undercolonizing","overcolonizing")) %>% 
  group_by(OTUID, category) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(OTUID) %>% 
  mutate(nCategories = n_distinct(category)) %>% 
  ungroup() %>% 
  mutate(OTUID = fct_reorder(OTUID, -nCategories)) %>%
  filter(nCombos > 1) %>% 
  ggplot() +
  geom_bar(aes(x=OTUID, y=nOTUs, fill=category), stat="identity") +
  scale_fill_manual(name="Category",
                    labels=c("overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("overcolonizing"="#D5A021",
                             "undercolonizing"="#3C91E6",
                             "dose-dependent"="#A63A43")) +
  xlab("ASV") +
  ylab("# mixtures") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT
p_numAllPatterns
save_plot("out/numAllPatterns.png", p_numAllPatterns, base_width=4, base_height=2)


# Plot distribution of categories of noisy ASVs.
p_numNoisyASVPatterns <- numASVPatterns %>% 
  filter(set=="noisy") %>%
  mutate(OTUID = fct_reorder(OTUID, -nSameSet)) %>% 
  ggplot() +
  geom_bar(aes(x=OTUID, y=nOTUs, fill=category), stat="identity") +
  scale_fill_manual(name="Category",
                    labels=c("noisy"="Noisy",
                             "flat"="Flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("noisy"="#999999",
                             "flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#3C91E6",
                             "dose-dependent"="#A63A43")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  DEFAULTS.THEME_PRES
p_numNoisyASVPatterns
save_plot("out/numNoisyASVPatterns.png", p_numNoisyASVPatterns, nrow=1.25, ncol=1.5)

# Plot distribution of categories of flat ASVs.
p_numFlatASVPatterns <- numASVPatterns %>% 
  filter(set=="flat") %>% 
  mutate(OTUID = fct_reorder(OTUID, -nSameSet)) %>% 
  ggplot() +
  geom_bar(aes(x=OTUID, y=nOTUs, fill=category), stat="identity") +
  scale_fill_manual(name="Category",
                    labels=c("noisy"="Noisy",
                             "flat"="Flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("noisy"="#999999",
                             "flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#3C91E6",
                             "dose-dependent"="#A63A43")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  DEFAULTS.THEME_PRES
p_numFlatASVPatterns
save_plot("out/numFlatASVPatterns.png", p_numFlatASVPatterns, nrow=1.25, ncol=1.5)

# Plot colonization pattern of each ASV across mixtures.
p_patternTiles <- dataColonization %>%
  mutate(combo=case_when(
    combo=="XBA-XBB" ~ "XBA\nXBB",
    combo=="XCA-XCB" ~ "XCA\nXCB",
    combo=="XDA-XDB" ~ "XDA\nXDB",
    combo=="XFA-XFB" ~ "XFA\nXFB",
    combo=="XBA-XCA" ~ "XBA\nXCA",
    combo=="XCA-XDA" ~ "XCA\nXDA",
    combo=="XDA-XFA" ~ "XDA\nXFA",
    combo=="XFA-XBA" ~ "XFA\nXBA"
    ),
    combo = fct_relevel(
      combo, c("XBA\nXBB","XCA\nXCB","XDA\nXDB","XFA\nXFB","XBA\nXCA","XCA\nXDA","XDA\nXFA","XFA\nXBA"))) %>% 
  filter(nCombos > 1) %>% 
  group_by(OTUIDshort) %>% 
  mutate(diffPatterns = n_distinct(category)) %>% 
  ungroup() %>% 
  mutate(OTUIDshort = fct_reorder(OTUIDshort, nCombos)) %>% 
  ggplot() +
  geom_tile(aes(x=combo, y=OTUIDshort, fill=category), color="black") +
  scale_fill_manual(name="Category",
                    labels=c("noisy"="Noisy",
                             "flat"="Flat",
                             "overcolonizing"="Overcolonizing",
                             "undercolonizing"="Undercolonizing",
                             "dose-dependent"="Dose-dependent"),
                    values=c("noisy"="#999999",
                             "flat"="#616667",
                             "overcolonizing"="#D5A021",
                             "undercolonizing"="#3C91E6",
                             "dose-dependent"="#A63A43")) +
  xlab("Mixture") +
  ylab("ASV")+
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT
  #facet_wrap(~nCombos, ncol=1)
  #facet_wrap(~diffPatterns, scales="free")
p_patternTiles
save_plot("out/patternTilesRBY.pdf", p_patternTiles, base_width=2.75, base_height=3)
#save_plot("out/patternTilesRBY.png", p_patternTiles, nrow=2, ncol=1.5)

# Plot colonization patterns of each ASV as ASV on Y axis and pattern on X.
p_patternMatching <- dataColonization %>% 
  #filter(nCombos>1) %>% 
  mutate(present="Present") %>% 
  select(!combo) %>% 
  group_by(OTUIDshort) %>% 
  mutate(nCategories = n_distinct(category)) %>% 
  ungroup() %>% 
  unique() %>% 
  complete(nesting(OTUIDshort, nCombos, nCategories), category, fill=list(present="Absent")) %>% 
  mutate(present=fct_relevel(present, c("Present","Absent")),
         nCombosBin = case_when(nCombos==1 ~ "1 mixture",
                                nCombos==2 ~ "2 mixtures",
                                nCombos>2 ~ "3+ mixtures"),
         category = case_when(category=="dose-dependent" ~ "DD",
                              category=="overcolonizing" ~ "OC",
                              category=="undercolonizing" ~ "UC"),
         category = fct_relevel(category, c("DD","OC","UC")),
         OTUIDshort = fct_reorder(OTUIDshort, nCategories)) %>% 
  ggplot() +
  geom_tile(aes(x=category, y=OTUIDshort, fill=present), color="black") +
  scale_fill_manual(values=c("Present"="#28112B","Absent"="#7D7D7D")) +
  xlab("Pattern") +
  ylab("ASV") +
  #coord_equal() +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  facet_wrap(~nCombosBin, scales="free") +
  DEFAULTS.THEME_PRINT
p_patternMatching
save_plot("out/patternMatching.pdf", p_patternMatching, base_width=5, base_height=3)

# Plot colonization patterns of each ASV as ASV on Y axis and pattern on X.
p_patternMatching2combo <- dataColonization %>% 
  filter(nCombos==2) %>% 
  #filter(nCombos>=3)
  mutate(present="Present") %>% 
  select(!combo) %>% 
  group_by(OTUIDshort) %>% 
  mutate(nCategories = n_distinct(category)) %>% 
  ungroup() %>% 
  unique() %>% 
  complete(nesting(OTUIDshort, nCombos, nCategories), category, fill=list(present="Absent")) %>% 
  mutate(present=fct_relevel(present, c("Present","Absent")),
         nCombosBin = case_when(nCombos==1 ~ "1 mixture",
                                nCombos==2 ~ "2 mixtures",
                                nCombos>2 ~ "3+ mixtures"),
         category = case_when(category=="dose-dependent" ~ "DD",
                              category=="overcolonizing" ~ "OC",
                              category=="undercolonizing" ~ "UC"),
         category = fct_relevel(category, c("DD","OC","UC")),
         OTUIDshort = fct_reorder(OTUIDshort, nCategories)) %>% 
  ggplot() +
  geom_tile(aes(x=category, y=OTUIDshort, fill=present), color="black") +
  scale_fill_manual(values=c("Present"="#28112B","Absent"="#7D7D7D")) +
  xlab("Pattern") +
  ylab("ASV") +
  #coord_equal() +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~nCombosBin, scales="free") +
  DEFAULTS.THEME_PRINT
p_patternMatching2combo
save_plot("out/patternMatching2combo.pdf", p_patternMatching2combo, base_width=1.4, base_height=2)

patternMatchingLegend <- get_legend(p_patternMatching1combo)
save_plot("out/patternMatchingLegend.pdf", patternMatchingLegend, base_width=0.5, base_height=0.5)

# Plot tile colonization patterns of each ASV with number per pattern.
p_patternMatchingNumbers3combo <- dataColonization %>% 
  filter(nCombos>=3) %>% 
  group_by(OTUIDshort, category) %>% 
  summarize(count = n(), nCombos) %>% 
  unique() %>% 
  ungroup() %>% 
  group_by(OTUIDshort) %>% 
  mutate(nCategories = n_distinct(category)) %>%
  ungroup() %>% 
  complete(nesting(OTUIDshort, nCombos, nCategories), category, fill=list(count=0)) %>% 
  mutate(nCombosBin = case_when(nCombos==1 ~ "1 mixture",
                                nCombos==2 ~ "2 mixtures",
                                nCombos>2 ~ "3+ mixtures"),
         category = case_when(category=="dose-dependent" ~ "DD",
                              category=="overcolonizing" ~ "OC",
                              category=="undercolonizing" ~ "UC"),
         category = fct_relevel(category, c("DD","OC","UC")),
         OTUIDshort = fct_reorder(OTUIDshort, nCategories)) %>% 
  ggplot() +
  geom_tile(aes(x=category, y=OTUIDshort, fill=as.character(count)), color="black") +
  scale_fill_manual(values=c("4"="#28112B",
                             "3"="#441D49",
                             "2"="#5F2966",
                             "1"="#7A3483",
                             "0"="#7D7D7D")) +
  xlab("Pattern") +
  ylab("ASV") +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  facet_wrap(~nCombosBin, scales="free") +
  DEFAULTS.THEME_PRINT
p_patternMatchingNumbers3combo
save_plot("out/patternMatchingNumbers3combo.pdf", p_patternMatchingNumbers3combo, base_width=1.4, base_height=2)
  
patternMatchingNumbersLegend <- get_legend(p_patternMatchingNumbers3combo)
save_plot("out/patternMatchingNumbersLegend.pdf", patternMatchingNumbersLegend, base_width=0.25, base_height=0.5)

# Pattern matching heatmap with dendrogram. ------------------------------------------------

# Import phyloseq to get all ASVs.
ps <- readRDS("../../workflow/out/e0017/DADA2_output/ps_all.rds")
df_tax_raw <- as.data.frame(tax_table(ps))
allASVs <- unique(rownames(df_tax_raw))

# Import phylogeny.
tree <- read.tree("../../workflow/out/e0017/DADA2_output/dsvs_msa.tree")

# Get ASVs not in DD, OC, or UC categories.
includedASVs <- data %>% 
  filter(passage==5 & category %in% c("dose-dependent","overcolonizing","undercolonizing")) %>% 
  pull(OTU) %>% 
  unique()
excludedASVs <- allASVs[which(!(allASVs %in% includedASVs))]

# Subset tree to only DD, OC, UC ASVs.
trimmedTree <- drop.tip(tree, excludedASVs)

# Change tip labels.
treedf <- as_tibble(trimmedTree) %>% 
  left_join(
    data %>% 
      ungroup() %>% 
      filter(OTU %in% includedASVs) %>% 
      select(OTU, Family, OTUnum) %>% 
      #left_join(familyCoding %>% rename(Family=family)) %>% 
      mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
      #mutate(OTUID = paste0(code, "-", OTUnum)) %>% 
      unique() %>% 
      select(OTU, OTUID, Family) %>% 
      rename(label = OTU)
      #left_join(taxaPalette %>% select(hex, Family))
  ) %>% 
  mutate(label = OTUID) %>% 
  select(!OTUID)
trimmedLabeledTree <- as.phylo(treedf)

# Generate dendrogram for plotting.
dendro <- as.dendrogram(trimmedLabeledTree)
dData <- dendro_data(dendro, type="rectangle")

dendroLabelsHex <- label(dData) %>%
  mutate(taxashort = gsub("-.*","",label)) %>% 
  left_join(palette %>% select(hex, taxashort)) %>% 
  pull(hex)

p_dendroNew <- ggplot(segment(dData)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  scale_y_reverse(expand=c(0.2,0)) +
  ylab("Phylogenetic distance") +
  xlab("ASV") +
  theme(#axis.ticks.x=element_blank(),
        #axis.line.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.title.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(vjust=-5)) +
  DEFAULTS.THEME_PRINT
p_dendroNew
save_plot("out/DD-OC-UC-ASVsdendrogram-new.pdf", p_dendroNew, base_width=1.4, base_height=3.97)

# Plot dendrogram.
p_dendro <- ggdendrogram(data=dendro, rotate=TRUE) +
  scale_y_reverse() +
  #scale_y_discrete(position="right") +
  #theme(axis.text.y=element_text(size=6)) +
  #geom_text(data=label(dData), aes(label=label, x=x, y=-1.5, color=label(dData)$code), size=3) +
  theme(axis.text.x=element_text(hjust=1),
        axis.text.y=element_blank()) +
  DEFAULTS.THEME_PRINT
p_dendro
save_plot("out/DD-OC-UC-ASVsdendrogram.pdf", p_dendro, base_width=1.5, base_height=4)

# Get order of ASVs from dendrogram.
dendroData <- dendro_data(dendro)
dendro_order <- (dendroData$labels)$label

# Plot colonization patterns of each ASV as ASV on Y axis and pattern on X, all ASVS in dendrogram.
p_patternMatchingAllASVs <- dataColonization %>% 
  group_by(OTUID, category) %>% 
  summarize(count=n()) %>% 
  unique() %>% 
  ungroup() %>% 
  complete(OTUID, category, fill=list(count=0)) %>% 
  mutate(category = case_when(category=="dose-dependent" ~ "DD",
                              category=="overcolonizing" ~ "S",
                              category=="undercolonizing" ~ "W"),
         category = fct_relevel(category, c("DD","OC","UC")),
         OTUID = fct_relevel(OTUID, dendro_order)) %>% 
  ggplot() +
  geom_tile(aes(x=category, y=OTUID, fill=as.character(count)), color="black") +
  scale_fill_manual(values=c("4"="#0E060F",
                             "3"="#602966",
                             "2"="#BC6DC5",
                             "1"="#E3C4E8",
                             "0"="#7E7E7E")) +
  xlab("Pattern") +
  ylab("ASV") +
  #coord_equal() +
  guides(fill = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  #theme(legend.position="none") +
  theme(axis.text.y=element_text(color = dendroLabelsHex),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.line.x=element_blank()) +
  #facet_wrap(~nCombosBin, scales="free") +
  DEFAULTS.THEME_PRINT
p_patternMatchingAllASVs
save_plot("out/patternMatchingAllASVs.pdf", p_patternMatchingAllASVs, base_width=1.8, base_height=3.77)

patternMatchingLegend <- get_legend(p_patternMatchingAllASVs)
save_plot("out/patternMatchingLegend.pdf", patternMatchingLegend, base_width=0.25, base_height=0.5)
