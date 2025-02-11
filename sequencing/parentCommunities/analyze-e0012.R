library(dada2)
library(phyloseq)
library(tidyverse)
library(ggfittext)
library(ggtext)
library(seqinr)
library(cowplot)
library(exactRankTests)
library(rlist)
library(foreach)
library(RColorBrewer)
library(patchwork)
library(data.table)
library(scales)

# Set working directory for running locally.
setwd("../")
# Set the plot theme.
theme_set(theme_cowplot())
source("config/palettes/plotDefaults.R")

plotdir <- "analysis/e0012Plots"

# Import read data.
ps <- readRDS("config/e0012-ps_all.rds")
# Import samplesheet.
ss <- read_tsv("config/20210201-e0012-samplesheet.txt")
# Import taxa color palette.
palette <- read.table("config/palettes/KCHcolors-Silva-partial.txt",
                      header=TRUE, stringsAsFactors = FALSE) %>% 
  mutate(taxashort=gsub("*.*\\.","",taxa))

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
  left_join(ss, by="filename") %>% 
  filter(count>0 & (media=="mBHImucin" | passage==0) & PCRreplicate==1) %>% 
  group_by(filename) %>% 
  mutate(rel_abundance=count/sum(count)) %>% 
  group_by(filename, Family) %>% 
  mutate(FamilyRel_abundance=sum(rel_abundance)) %>% 
  ungroup()

# Add passage 0 to all replicates.
df <- df %>% 
  rbind(df %>% filter(passage==0) %>% mutate(inoculationreplicate=2)) %>% 
  rbind(df %>% filter(passage==0) %>% mutate(inoculationreplicate=3))

dfHighAbundance <- df %>% 
  filter(FamilyRel_abundance > 0.01)
  #filter(rel_abundance > 0.01)

# Create stripped-down color palette containing only families present in the dataset.
palettee0012 <- palette %>% 
  filter(taxashort %in% sort(unique(dfHighAbundance$Family)))
# Make a named list.
palettee0012Vector <- palettee0012$hex
names(palettee0012Vector) <- palettee0012$taxashort

# Import family coding and bind to dataframe.
familyCoding <- fread("config/familyCodingFull.txt")
df <- df %>% left_join(familyCoding %>% rename(Family=family))

# Create version of palette to use with family codes.
palettee0012Codes <- left_join(palettee0012, familyCoding %>% rename(taxashort=family))
palettee0012VectorCodes <- palettee0012Codes$hex
names(palettee0012VectorCodes) <- palettee0012Codes$code

# Calculate alpha diversity metrics.
dfAlpha <- estimate_richness(ps, measures="Shannon")
dfAlpha <- dfAlpha %>% 
  mutate(sample = rownames(dfAlpha),
         sample = sub(".*_","",sample)) %>% 
  left_join(ss, by="sample") %>% 
  filter((media=="mBHImucin" | passage==0) & PCRreplicate==1)
dfAlpha <- dfAlpha %>% 
  rbind(dfAlpha %>% filter(passage==0) %>% mutate(inoculationreplicate=2)) %>% 
  rbind(dfAlpha %>% filter(passage==0) %>% mutate(inoculationreplicate=3))

# Calculate JSD.
JSD <- distance(ps, "jsd")
JSD_df <- melt(as.matrix(JSD), varnames=c("sample1","sample2")) %>% rename(JSD=value)

# Parse sample metadata from JSD df.
JSD_df_passages <- JSD_df %>% 
  mutate(sample1 = sub(".*_","",sample1),
         sample2 = sub(".*_","",sample2)) %>% 
  separate(sample1, into=c("subject1","timepoint1","media1","rep1","passage1","PCRrep1"), sep="-") %>% 
  separate(sample2, into=c("subject2","timepoint2","media2","rep2","passage2","PCRrep2"), sep="-") %>% 
  filter(subject1==subject2 & timepoint1==timepoint2 & rep1==rep2 & PCRrep1==PCRrep2 & passage1!=passage2) %>% 
  filter(((media1=="mBHImucin" & media2=="inoculum") | (media1=="inoculum" & media2=="mBHImucin") | (media1=="mBHImucin" & media2=="mBHImucin")) & 
           PCRrep1=="PCR1") %>%
  mutate(passage1 = case_when(passage1=="p0" ~ 0,
                              passage1=="p1" ~ 1,
                              passage1=="p3" ~ 3,
                              passage1=="p7" ~ 7,
                              passage1=="p15" ~ 15),
         passage2 = case_when(passage2=="p0" ~ 0,
                              passage2=="p1" ~ 1,
                              passage2=="p3" ~ 3,
                              passage2=="p7" ~ 7,
                              passage2=="p15" ~ 15)) %>% 
  mutate(comparison = paste0(passage1, "-", passage2))

# Get dataframe with passage 15 of all media + inoculum.
dfAllMedia <- left_join(df_otu, df_tax) %>% 
  mutate(filename=gsub("_.*", "", sample)) %>% 
  select(-sample) %>% 
  left_join(ss, by="filename") %>% 
  filter(count>0 & (passage==15 | media=="inoculum") & PCRreplicate==1) %>% 
  group_by(filename) %>% 
  mutate(rel_abundance=count/sum(count))

# JSD stability plots. ----------------------------------------------------

p_JSD <- JSD_df_passages %>%
  filter(comparison %in% c("0-1","1-3","3-7","7-15")) %>% 
  rbind(JSD_df_passages %>% filter(comparison=="0-1") %>% mutate(rep1="r2"),
        JSD_df_passages %>% filter(comparison=="0-1") %>% mutate(rep1="r3")) %>% 
  mutate(comparison = fct_relevel(comparison, c("0-1","1-3","3-7","7-15"))) %>% 
  ggplot() +
  geom_line(aes(x=as.character(comparison), y=JSD, group=interaction(subject1, rep1)), size=0.25) +
  ylim(0,0.65) +
  xlab("Comparison between passages") +
  ylab("JSD") +
  DEFAULTS.THEME_PRINT
p_JSD
save_plot(paste0(plotdir, "/p_JSDoverPassages.pdf"), p_JSD, base_width=1.7, base_height=1.2)

# Relative abundance barplots. --------------------------------------------

# Plot family-level relative abundance across passages.
p_familyAbundances <- df %>% 
  filter(FamilyRel_abundance>10^-2) %>% 
  ggplot() +
  geom_bar(aes(x=as.character(passage), y=rel_abundance, fill=Family), stat="identity", color="black", size=0.01) +
  scale_fill_manual(values=palettee0012Vector) +
  scale_x_discrete("Passage",limits=c("0","1","3","7","15")) +
  ylab("Relative abundance") +
  #theme(legend.position="none") +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        legend.text=element_text(size=5.5)) +
  guides(fill = guide_legend(ncol = 4)) +
  DEFAULTS.THEME_PRES +
  facet_grid(~biosample~inoculationreplicate)
p_familyAbundances
save_plot(paste0(plotdir,"/p_familyAbundances.png"), p_familyAbundances, nrow=1.75, ncol=2)

e0012Legend <- get_legend(p_familyAbundances)
save_plot(paste0(plotdir, "/e0012Legend.pdf"), e0012Legend, base_width=4.75, base_height=1)

p_familyAbundancesExampleXBA <- df %>% 
  filter(biosample=="XBA-029" & inoculationreplicate==1) %>%
  mutate(biosample = "A1") %>% 
  ggplot() +
  geom_bar(aes(x=as.character(passage), y=rel_abundance, fill=Family), stat="identity", color="black", size=0.05) +
  scale_fill_manual(values=palettee0012Vector) +
  scale_x_discrete("Passage",limits=c("0","1","3","7","15")) +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~biosample)
p_familyAbundancesExampleXBA
save_plot(paste0(plotdir,"/p_familyAbundancesExampleXBA-small.pdf"), p_familyAbundancesExampleXBA, base_width=1.75, base_height=1.5)
save_plot(paste0(plotdir,"/p_familyAbundancesExampleXBA-small.png"), p_familyAbundancesExampleXBA, base_width=1.5, base_height=1.2)

p_familyAbundancesExampleXCA <- df %>% 
  filter(biosample=="XCA-029" & inoculationreplicate==1) %>%
  mutate(biosample = "B1") %>% 
  ggplot() +
  geom_bar(aes(x=as.character(passage), y=rel_abundance, fill=Family), stat="identity", color="black", size=0.05) +
  scale_fill_manual(values=palettee0012Vector) +
  scale_x_discrete("Passage",limits=c("0","1","3","7","15")) +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~biosample)
p_familyAbundancesExampleXCA
save_plot(paste0(plotdir,"/p_familyAbundancesExampleXCA.pdf"), p_familyAbundancesExampleXCA, base_width=1.75, base_height=1.5)

p_familyAbundancesExampleXCB <- df %>% 
  filter(biosample=="XCB-029" & inoculationreplicate==1) %>%
  mutate(biosample = "B2") %>% 
  ggplot() +
  geom_bar(aes(x=as.character(passage), y=rel_abundance, fill=Family), stat="identity", color="black", size=0.05) +
  scale_fill_manual(values=palettee0012Vector) +
  scale_x_discrete("Passage",limits=c("0","1","3","7","15")) +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~biosample)
p_familyAbundancesExampleXCB
save_plot(paste0(plotdir,"/p_familyAbundancesExampleXCB.pdf"), p_familyAbundancesExampleXCB, base_width=1.75, base_height=1.5)

p_familyAbundancesExampleXDA <- df %>% 
  filter(biosample=="XDA-029" & inoculationreplicate==1) %>%
  mutate(biosample = "C1") %>% 
  ggplot() +
  geom_bar(aes(x=as.character(passage), y=rel_abundance, fill=Family), stat="identity", color="black", size=0.05) +
  scale_fill_manual(values=palettee0012Vector) +
  scale_x_discrete("Passage",limits=c("0","1","3","7","15")) +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~biosample)
p_familyAbundancesExampleXDA
save_plot(paste0(plotdir,"/p_familyAbundancesExampleXDA.pdf"), p_familyAbundancesExampleXDA, base_width=1.75, base_height=1.5)

p_familyAbundancesExampleXDB <- df %>% 
  filter(biosample=="XDB-029" & inoculationreplicate==1) %>%
  mutate(biosample = "C2") %>% 
  ggplot() +
  geom_bar(aes(x=as.character(passage), y=rel_abundance, fill=Family), stat="identity", color="black", size=0.05) +
  scale_fill_manual(values=palettee0012Vector) +
  scale_x_discrete("Passage",limits=c("0","1","3","7","15")) +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~biosample)
p_familyAbundancesExampleXDB
save_plot(paste0(plotdir,"/p_familyAbundancesExampleXDB.pdf"), p_familyAbundancesExampleXDB, base_width=1.75, base_height=1.5)

p_familyAbundancesExampleXFA <- df %>% 
  filter(biosample=="XFA-029" & inoculationreplicate==1) %>%
  mutate(biosample = "D1") %>% 
  ggplot() +
  geom_bar(aes(x=as.character(passage), y=rel_abundance, fill=Family), stat="identity", color="black", size=0.05) +
  scale_fill_manual(values=palettee0012Vector) +
  scale_x_discrete("Passage",limits=c("0","1","3","7","15")) +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~biosample)
p_familyAbundancesExampleXFA
save_plot(paste0(plotdir,"/p_familyAbundancesExampleXFA.pdf"), p_familyAbundancesExampleXFA, base_width=1.75, base_height=1.5)

p_familyAbundancesExampleXFB <- df %>% 
  filter(biosample=="XFB-029" & inoculationreplicate==1) %>%
  mutate(biosample = "D2") %>% 
  ggplot() +
  geom_bar(aes(x=as.character(passage), y=rel_abundance, fill=Family), stat="identity", color="black", size=0.05) +
  scale_fill_manual(values=palettee0012Vector) +
  scale_x_discrete("Passage",limits=c("0","1","3","7","15")) +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~biosample)
p_familyAbundancesExampleXFB
save_plot(paste0(plotdir,"/p_familyAbundancesExampleXFB.pdf"), p_familyAbundancesExampleXFB, base_width=1.75, base_height=1.5)

# Get abundant families from XBA-029.
familiesAbundantExample <- (df %>% 
                              filter(biosample=="XBA-029") %>% 
                              group_by(taxa) %>% 
                              summarize(maxFamilyAbundance=max(FamilyRel_abundance)) %>% 
                              filter(maxFamilyAbundance>=0.005))$taxa
# Export abundant families for taxa key generation.
write.table(familiesAbundantExample, file=paste0(plotdir,"/familiesAbundantExample.tsv"), quote=FALSE, row.names=FALSE)

# Get abundant families from all samples.
familiesAbundantAlle0012 <- df %>%
  group_by(taxa) %>% 
  summarize(maxFamilyAbundance = max(FamilyRel_abundance)) %>% 
  filter(maxFamilyAbundance>=0.005) %>% 
  pull(taxa)
write.table(familiesAbundantAlle0012, file=paste0(plotdir, "/familiesAbundantAlle0012.txt"), 
            quote=FALSE, row.names=FALSE, sep="\t")

# Plot abundances for one community example across media.
p_allMediaExample <- dfAllMedia %>% 
  filter(biosample=="XBA-029" & inoculationreplicate==1) %>% 
  ungroup() %>% 
  mutate(media=ifelse(media=="mBHImucin","mBHI+mucin",media),
         media=fct_relevel(media, c("inoculum","BHI","mBHI","mBHI+mucin","TYG")),
         biosample="A1")%>% 
  ggplot() +
  geom_bar(aes(x=media, y=rel_abundance, fill=Family), stat="identity", color="black", size=0.05) +
  scale_fill_manual(values=palettee0012Vector) +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~biosample)
p_allMediaExample
save_plot(paste0(plotdir, "/p_allMediaExample.pdf"), p_allMediaExample, base_width=2.5, base_height=1.75)

# Relative abundance correlations. ----------------------------------------

# Plot correlation of abundances between passages 0 and 3.
abundanceP0P3_df <- df %>% 
  # Turn dataframe into wide format for plotting.
  filter(passage==0) %>%
  select(Family, biosample, inoculationreplicate, FamilyRel_abundance) %>% 
  rename(FamilyRelAbundance_0=FamilyRel_abundance) %>% 
  unique() %>% 
  left_join(
    df %>%
      filter(passage==3) %>% 
      select(Family, biosample, inoculationreplicate, FamilyRel_abundance) %>% 
      rename(FamilyRelAbundance_3=FamilyRel_abundance) %>% 
      unique(),
    by=c("Family","biosample","inoculationreplicate")
  ) %>% 
  # Replace NA values with 10^-4 relative abundance.
  mutate(FamilyRelAbundance_0=replace(FamilyRelAbundance_0, is.na(FamilyRelAbundance_0), 10^-3),
         FamilyRelAbundance_3=replace(FamilyRelAbundance_3, is.na(FamilyRelAbundance_3), 10^-3)) %>% 
  group_by(Family, biosample) %>%
  summarize(Family, biosample, 
            logFamilyRelAbundance_0_avg=log10(mean(FamilyRelAbundance_0)),
            logFamilyRelAbundance_0_max=max(log10(FamilyRelAbundance_0)),
            logFamilyRelAbundance_0_min=min(log10(FamilyRelAbundance_0)),
            logFamilyRelAbundance_0_range=logFamilyRelAbundance_0_max-logFamilyRelAbundance_0_min,
            logFamilyRelAbundance_3_avg=log10(mean(FamilyRelAbundance_3)),
            logFamilyRelAbundance_3_max=max(log10(FamilyRelAbundance_3)),
            logFamilyRelAbundance_3_min=min(log10(FamilyRelAbundance_3)),
            logFamilyRelAbundance_3_range=logFamilyRelAbundance_3_max-logFamilyRelAbundance_3_min) %>% 
  unique()

# Plot scatterplot of relative abundances for p0 vs p3.
p_abundanceScatterP0P3 <- abundanceP0P3_df %>% 
  ggplot() +
  geom_point(aes(x=logFamilyRelAbundance_0_avg, y=logFamilyRelAbundance_3_avg, color=Family)) +
  geom_errorbar(aes(x=logFamilyRelAbundance_0_avg, 
                    ymin=logFamilyRelAbundance_3_min,
                    ymax=logFamilyRelAbundance_3_max,
                    color=Family), alpha=0.25) +
  geom_errorbarh(aes(y=logFamilyRelAbundance_3_avg, 
                     xmin=logFamilyRelAbundance_0_min,
                     xmax=logFamilyRelAbundance_0_max,
                     color=Family), alpha=0.25) +
  geom_abline(slope=1, linetype="dotted", color="gray") +
  scale_color_manual(values=palettee0012Vector) +
  xlab("Passage 0 relative abundance") +
  ylab("Passage 3 relative abundance") +
  xlim(c(-3,0)) +
  ylim(c(-3,0)) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~biosample, nrow=2)
p_abundanceScatterP0P3
save_plot(paste0(plotdir,"/p_abundanceScatterP0P3.png"), p_abundanceScatterP0P3, nrow=1.25, ncol=1.5)

# Plot scatterplot of relative abundances for p0 vs p3 example.
p_abundanceScatterP0P3Example <- abundanceP0P3_df %>% 
  filter(biosample=="XBA-029") %>% 
  mutate(biosample = "A1") %>% 
  ggplot() +
  geom_point(aes(x=logFamilyRelAbundance_0_avg, y=logFamilyRelAbundance_3_avg, color=Family), size=1) +
  geom_errorbar(aes(x=logFamilyRelAbundance_0_avg, 
                    ymin=logFamilyRelAbundance_3_min,
                    ymax=logFamilyRelAbundance_3_max,
                    color=Family), alpha=0.25, size=0.5) +
  geom_errorbarh(aes(y=logFamilyRelAbundance_3_avg, 
                     xmin=logFamilyRelAbundance_0_min,
                     xmax=logFamilyRelAbundance_0_max,
                     color=Family), alpha=0.25, size=0.5) +
  geom_abline(slope=1, linetype="dashed", color="gray", size=0.5) +
  scale_color_manual(values=palettee0012Vector) +
  scale_x_continuous(name="Passage 0\nrelative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_y_continuous(name="Passage 3\nrelative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  #theme(axis.text.x=element_text(angle=90)) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~biosample, nrow=2)
p_abundanceScatterP0P3Example
save_plot(paste0(plotdir,"/p_abundanceScatterP0P3Example.png"), p_abundanceScatterP0P3Example, nrow=1, ncol=.7)

# Plot correlation of abundances between passages 3 and 7.
abundanceP3P7_df <- df %>% 
  # Turn dataframe into wide format for plotting.
  filter(passage==3) %>%
  select(Family, biosample, inoculationreplicate, FamilyRel_abundance) %>% 
  rename(FamilyRelAbundance_3=FamilyRel_abundance) %>% 
  unique() %>% 
  left_join(
    df %>%
      filter(passage==7) %>% 
      select(Family, biosample, inoculationreplicate, FamilyRel_abundance) %>% 
      rename(FamilyRelAbundance_7=FamilyRel_abundance) %>% 
      unique(),
    by=c("Family","biosample","inoculationreplicate")
  ) %>% 
  # Replace NA values with 10^-4 relative abundance.
  mutate(FamilyRelAbundance_3=replace(FamilyRelAbundance_3, is.na(FamilyRelAbundance_3), 10^-3),
         FamilyRelAbundance_7=replace(FamilyRelAbundance_7, is.na(FamilyRelAbundance_7), 10^-3)) %>% 
  group_by(Family, biosample) %>%
  summarize(Family, biosample, 
            logFamilyRelAbundance_3_avg=log10(mean(FamilyRelAbundance_3)),
            logFamilyRelAbundance_3_max=max(log10(FamilyRelAbundance_3)),
            logFamilyRelAbundance_3_min=min(log10(FamilyRelAbundance_3)),
            logFamilyRelAbundance_3_range=logFamilyRelAbundance_3_max-logFamilyRelAbundance_3_min,
            logFamilyRelAbundance_7_avg=log10(mean(FamilyRelAbundance_7)),
            logFamilyRelAbundance_7_max=max(log10(FamilyRelAbundance_7)),
            logFamilyRelAbundance_7_min=min(log10(FamilyRelAbundance_7)),
            logFamilyRelAbundance_7_range=logFamilyRelAbundance_7_max-logFamilyRelAbundance_7_min) %>% 
  unique()

# Plot scatterplot of relative abundances for p3 vs p7.
p_abundanceScatterP3P7 <- abundanceP3P7_df %>% 
  ggplot() +
  geom_point(aes(x=logFamilyRelAbundance_3_avg, y=logFamilyRelAbundance_7_avg, color=Family)) +
  geom_errorbar(aes(x=logFamilyRelAbundance_3_avg, 
                    ymin=logFamilyRelAbundance_7_min,
                    ymax=logFamilyRelAbundance_7_max,
                    color=Family), alpha=0.25) +
  geom_errorbarh(aes(y=logFamilyRelAbundance_7_avg, 
                     xmin=logFamilyRelAbundance_3_min,
                     xmax=logFamilyRelAbundance_3_max,
                     color=Family), alpha=0.25) +
  geom_abline(slope=1, linetype="dotted", color="gray") +
  scale_color_manual(values=palettee0012Vector) +
  xlab("Passage 3 relative abundance") +
  ylab("Passage 7 relative abundance") +
  xlim(c(-3,0)) +
  ylim(c(-3,0)) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~biosample, nrow=2)
p_abundanceScatterP3P7
save_plot(paste0(plotdir,"/p_abundanceScatterP3P7.png"), p_abundanceScatterP3P7, nrow=1.25, ncol=1.5)

# Plot scatterplot of relative abundances for p3 vs p7 example.
p_abundanceScatterP3P7Example <- abundanceP3P7_df %>% 
  filter(biosample=="XBA-029") %>% 
  mutate(biosample = "A1") %>% 
  ggplot() +
  geom_point(aes(x=logFamilyRelAbundance_3_avg, y=logFamilyRelAbundance_7_avg, color=Family), size=1) +
  geom_errorbar(aes(x=logFamilyRelAbundance_3_avg, 
                    ymin=logFamilyRelAbundance_7_min,
                    ymax=logFamilyRelAbundance_7_max,
                    color=Family), alpha=0.25, size=0.5) +
  geom_errorbarh(aes(y=logFamilyRelAbundance_7_avg, 
                     xmin=logFamilyRelAbundance_3_min,
                     xmax=logFamilyRelAbundance_3_max,
                     color=Family), alpha=0.25, size=0.5) +
  geom_abline(slope=1, linetype="dashed", color="gray", size=0.5) +
  scale_color_manual(values=palettee0012Vector) +
  scale_x_continuous(name="Passage 3\nrelative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_y_continuous(name="Passage 7\nrelative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~biosample, nrow=2)
p_abundanceScatterP3P7Example
save_plot(paste0(plotdir,"/p_abundanceScatterP3P7Example.png"), p_abundanceScatterP3P7Example, nrow=1, ncol=.7)

save_plot(paste0(plotdir,"/p_abundanceScatterCombinedExampleP0P3P7.pdf"), 
          p_abundanceScatterP0P3Example + p_abundanceScatterP3P7Example,
          base_width=3.2,base_height=1.8)

# Calculate correlation of relative abundances between p0 and p3 / p3 and p7.
p0p3corr <- cor.test(abundanceP0P3_df %>% filter(biosample=="XBA-029") %>% pull(logFamilyRelAbundance_3_avg),
                     abundanceP0P3_df %>% filter(biosample=="XBA-029") %>% pull(logFamilyRelAbundance_0_avg),
                     method = "pearson")

p3p7corr <- cor.test(abundanceP3P7_df %>% filter(biosample=="XBA-029") %>% pull(logFamilyRelAbundance_7_avg),
                     abundanceP3P7_df %>% filter(biosample=="XBA-029") %>% pull(logFamilyRelAbundance_3_avg),
                     method="pearson")

# ASV diversity. ----------------------------------------------------------

# Count ASVs and family-level diversity above 10^-3 relative abundance.
diversitydf <- df %>% 
  filter(rel_abundance>=10^-3) %>% 
  group_by(biosample, inoculationreplicate, passage) %>% 
  summarize(ASVcount=n(),
            familyDiversity=n_distinct(Family))

p_ASVcounts <- diversitydf %>% 
  ggplot() +
  geom_line(aes(x=as.character(passage), y=ASVcount, group=interaction(biosample, inoculationreplicate)),
            size=0.25, color="#E41A1C") +
  scale_x_discrete("Passage", limits=c("0","1","3","7","15")) +
  ylab("# ASVs") +
  ylim(c(0,125)) +
  DEFAULTS.THEME_PRINT
p_ASVcounts
save_plot(paste0(plotdir,"/p_ASVcounts.pdf"), p_ASVcounts, base_width=1.75, base_height=1.2)

p_familyCounts <- diversitydf %>% 
  ggplot() +
  geom_line(aes(x=as.character(passage), y=familyDiversity, group=interaction(biosample, inoculationreplicate)),
            size=0.25, color="#377EB8") +
  scale_x_discrete("Passage", limits=c("0","1","3","7","15")) +
  ylab("# families") +
  ylim(c(0,25)) +
  DEFAULTS.THEME_PRINT
p_familyCounts
save_plot(paste0(plotdir, "/p_familyCounts.pdf"), p_familyCounts, base_width=1.75, base_height=1.2)

# Shanon diversity. -------------------------------------------------------

# Plot Shannon effective # species over passages.
p_shannon <- dfAlpha %>% 
  ggplot() +
  geom_line(aes(x=as.character(passage), y=exp(Shannon), group=interaction(biosample, inoculationreplicate)),
            size=0.25, color="#E41A1C") +
  scale_x_discrete("Passage", limits=c("0","1","3","7","15")) +
  ylab("Effective # species") +
  ylim(c(0,80)) +
  DEFAULTS.THEME_PRINT
p_shannon
save_plot(paste0(plotdir, "/p_shannon.pdf"), p_shannon, base_width=1.5, base_height=1.2)

# Recapitulated diversity. ------------------------------------------------

# Calculate % passage 0 relative abundance recapitulated at different passages for both ASV and family.
recapAbundance_df <- df %>% 
  #filter(rel_abundance>=10^-3) %>%
  filter(passage==0 | rel_abundance>=10^-3) %>% 
  mutate(P0 = passage==0) %>% 
  group_by(OTU, biosample, inoculationreplicate) %>% 
  mutate(inP0 = TRUE %in% P0) %>% 
  group_by(Family, biosample, inoculationreplicate) %>% 
  summarize(OTU, passage, inP0, FamilyinP0 = TRUE %in% P0, rel_abundance, FamilyRel_abundance) %>% 
  filter(passage!=0 & (inP0==TRUE | FamilyinP0==TRUE)) %>% 
  mutate(rel_abundance=replace(rel_abundance, inP0==FALSE, NA)) %>% 
  group_by(biosample, inoculationreplicate, passage) %>% 
  summarize(FamilyRel_abundance, ASV=sum(rel_abundance, na.rm=TRUE)) %>%
  unique() %>% 
  group_by(biosample, inoculationreplicate, passage) %>% 
  summarize(ASV, Family=sum(FamilyRel_abundance)) %>% 
  unique() %>% 
  pivot_longer(cols=c(ASV, Family), names_to="type", values_to="recapAbundance")

p_recapAbundance <- recapAbundance_df %>% 
  ggplot() +
  geom_jitter(aes(x=as.character(passage), y=recapAbundance, color=type), width=0.15, size=0.25)+
  geom_boxplot(aes(x=as.character(passage), y=recapAbundance), alpha=0, size=0.25) +
  scale_color_brewer(palette="Set1", name="Taxon rank") +
  scale_x_discrete("Passage", limits=c("1","3","7","15")) +
  ylab("Relative abundance \nof inoculum") +
  ylim(c(0,1)) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~type)
p_recapAbundance
save_plot(paste0(plotdir,"/p_recapAbundance.pdf"), p_recapAbundance, base_width=2.5, base_height=1.5)