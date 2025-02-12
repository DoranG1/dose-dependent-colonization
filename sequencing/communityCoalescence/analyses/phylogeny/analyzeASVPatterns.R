library(tidyverse)
library(phyloseq)
library(cowplot)
library(data.table)
library(foreach)
library(ape)
library(treeio)
library(ggdendro)
library(phylogram)

theme_set(theme_cowplot())
source("../../config/palettes/plotDefaults.R")

# Same data as other analyses, just read an intermediate instead of running the generateMixturesDataframe script.
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
  filter(taxashort %in% sort(unique(data %>% filter(rel_abundance > 0.001) %>% pull(Family))))
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

dataColonizationAllPatterns <- data %>% 
  filter(passage==5 & category!="lowAbundance") %>% 
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


# Pattern matching heatmap with dendrogram, all five categories. ----------


# Import phyloseq to get all ASVs.
ps <- readRDS("../../workflow/out/e0017/DADA2_output/ps_all.rds")
df_tax_raw <- as.data.frame(tax_table(ps))
allASVs <- unique(rownames(df_tax_raw))

# Import phylogeny.
tree <- read.tree("../../workflow/out/e0017/DADA2_output/dsvs_msa.tree")

# Identify ASVs not in DD, OC, UC, flat, or noisy categories.
includedASVsAllPatterns <- data %>% 
  filter(passage==5 & category!="lowAbundance") %>% 
  pull(OTU) %>% 
  unique()
excludedASVsAllPatterns <- allASVs[which(!(allASVs %in% includedASVsAllPatterns))]

# Subset tree to only included ASV patterns.
trimmedTreeAllPatterns <- drop.tip(tree, excludedASVsAllPatterns)

# Change tip labels.
treedfAllPatterns <- as_tibble(trimmedTreeAllPatterns) %>% 
  left_join(
    data %>% 
      ungroup() %>% 
      filter(OTU %in% includedASVsAllPatterns) %>% 
      select(OTU, Family, OTUnum) %>% 
      mutate(OTUID = paste0(Family, "-", OTUnum)) %>%
      unique() %>% 
      select(OTU, OTUID, Family) %>% 
      rename(label = OTU)
  ) %>% 
  mutate(label = OTUID) %>% 
  select(!OTUID)
trimmedLabeledTreeAllPatterns <- as.phylo(treedfAllPatterns)

# Generate dendrogram for plotting.
dendroAllPatterns <- as.dendrogram(trimmedLabeledTreeAllPatterns)
dDataAllPatterns <- dendro_data(dendroAllPatterns, type="rectangle")

dendroLabelsHexAllPatterns <- label(dDataAllPatterns) %>%
  mutate(taxashort = gsub("-.*","",label)) %>% 
  left_join(palette %>% select(hex, taxashort)) %>% 
  # Remove the second Peptostreptococcaceae coloring.
  filter(!(taxashort=="Peptostreptococcaceae" & hex=="#daa520")) %>% 
  pull(hex)

# Plot dendrogram.
p_dendroAllPatterns <- ggplot(segment(dDataAllPatterns)) +
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
p_dendroAllPatterns
save_plot("out/ASVsdendrogramAllPatterns.pdf", p_dendroAllPatterns, base_width=1.7, base_height=7.27)

# Get order of ASVs from dendrogram.
dendroDataAllPatterns <- dendro_data(dendroAllPatterns)
dendro_orderAllPatterns <- (dendroDataAllPatterns$labels)$label

# Plot colonization patterns of each ASV as ASV on Y axis and pattern on X, all ASVS in dendrogram.
p_patternMatchingAllASVsAllPatterns <- dataColonizationAllPatterns %>% 
  group_by(OTUID, category) %>% 
  summarize(count=n()) %>% 
  unique() %>% 
  ungroup() %>% 
  complete(OTUID, category, fill=list(count=0)) %>% 
  mutate(category = case_when(category=="dose-dependent" ~ "DD",
                              category=="overcolonizing" ~ "S",
                              category=="undercolonizing" ~ "W",
                              category=="flat" ~ "R",
                              category=="noisy" ~ "N"),
         category = fct_relevel(category, c("DD","S","W","R","N")),
         OTUID = fct_relevel(OTUID, dendro_orderAllPatterns)) %>% 
  ggplot() +
  geom_tile(aes(x=category, y=OTUID, fill=as.character(count)), color="black") +
  scale_fill_manual(values=c("8"="#0E060F",
                             "7"="#37183B",
                             "6"="#4C2151",
                             "5"="#602966",
                             "4"="#8E4B96",
                             "3"="#BC6DC5",
                             "2"="#D099D7",
                             "1"="#E3C4E8",
                             "0"="#7E7E7E")) +
  xlab("Pattern") +
  ylab("ASV") +
  #coord_equal() +
  guides(fill = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  theme(axis.text.y=element_text(color = dendroLabelsHexAllPatterns),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.line.x=element_blank()) +
  #facet_wrap(~nCombosBin, scales="free") +
  DEFAULTS.THEME_PRINT
p_patternMatchingAllASVsAllPatterns
save_plot("out/patternMatchingAllASVsAllPatterns.pdf", p_patternMatchingAllASVsAllPatterns, 
          base_width=2, base_height=6.75)

patternMatchingLegendAllPatterns <- get_legend(p_patternMatchingAllASVsAllPatterns)
save_plot("out/patternMatchingLegendAllPatterns.pdf", patternMatchingLegendAllPatterns, base_width=0.25, base_height=1)

