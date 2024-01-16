library(tidyverse)
library(phyloseq)
library(cowplot)
library(data.table)
library(foreach)
library(FAVA)
library(data.table)
library(ape)
library(ggrepel)

theme_set(theme_cowplot())
source("../../config/palettes/plotDefaults.R")

# Import data.
OTU_table <- as.data.frame(t(read.table("../OTUtableTheoretical.txt", header=TRUE, sep="\t")))
# Convert counts to relative abundances.
OTU_table <- OTU_table/rowSums(OTU_table)
# Add sample names as column.
OTU_table <- OTU_table %>% 
  mutate(sample = rownames(OTU_table),
         sample = gsub("\\.","-",sample),
         sample = sub("1-0","1:0",sample),
         sample = sub("1000-1","1000:1",sample),
         sample = sub("100-1","100:1",sample),
         sample = sub("10-1","10:1",sample),
         sample = sub("1-1","1:1",sample),
         sample = sub("1-10",'1:10',sample),
         sample = sub("1-100","1:100",sample),
         sample = sub("1-1000","1:1000",sample),
         sample = sub("0-1","0:1",sample))
# Import sample metadata.
sam_table <- fread("../samTableTheoretical.txt")

# Combine OTU and sample data.
df <- left_join(sam_table, OTU_table)
# Remove poorly sequenced samples.
dfFiltered <- df %>% 
  filter(!(combo %in% c("XBA-XBB", "XBA-XCA", "XFA-XBA") & passage==3 & line=="theoretical-1")) %>% 
  filter(!(passage==5 & line=="actual-1" & label=="XFA-100:1-XBA"))

# Read in phylogenetic tree.
tree <- read.tree("../../workflow/out/e0017/DADA2_output/dsvs_msa.tree")
# Import phyloseq to get all ASVs.
ps <- readRDS("../../workflow/out/e0017/DADA2_output/ps_all.rds")
df_tax_raw <- as.data.frame(tax_table(ps))
allASVs <- unique(rownames(df_tax_raw))

# Get ASVs not in data.
includedASVs <- colnames(OTU_table[,1:1362]) 
excludedASVs <- allASVs[which(!(allASVs %in% includedASVs))]

# Subset tree to only ASVs in data.
trimmedTree <- drop.tip(tree, excludedASVs)

# Create phylogenetic distance matrix between ASVs from tree.
distMatrix <- cophenetic.phylo(trimmedTree)

# Convert distance matrix to a similarity matrix.
simMatrix <- exp(-distMatrix)

# Calculate FAVA for each grouping of variables.
dfFAVA <- dfFiltered %>% 
  mutate(sample = paste0(line,"_",combo,"_",passage)) 
FAVA <- fava(dfFAVA, group="sample", K=1362) %>% 
  separate(sample, into=c("line","combo","passage"), sep="_", remove=FALSE) %>% 
  mutate(mixtureType = sub("-.*","",line))
FAVASummarized <- FAVA %>% 
  group_by(combo, passage, mixtureType) %>% 
  summarize(avgFAVA = mean(FAVA), maxFAVA = max(FAVA), minFAVA = min(FAVA))

# Import full dataframe.
data <- fread("../mixtureDataframe.txt")

# Get total relative abundance of dose-dependent ASVs in data at 1:1 mixture ratio, average across replicates.
DDRelAb <- data %>% 
  filter(category=="dose-dependent" & ratio=="1:1" & mixtureType=="actual") %>% 
  group_by(passage, combo, line) %>% 
  summarize(totalAbundance = sum(rel_abundance_old, na.rm=TRUE))
DDRelAbSummarized <- DDRelAb %>% 
  ungroup() %>% 
  group_by(passage, combo) %>% 
  summarize(maxAbundance = max(totalAbundance), minAbundance = min(totalAbundance), totalAbundance = mean(totalAbundance))

# Plot FAVA and total DD relative abundance of actual mixtures, average across replicates.
DDFAVASummarized <- DDRelAbSummarized %>% 
  mutate(passage = as.character(passage)) %>% 
  left_join(FAVASummarized %>% filter(mixtureType=="actual"))

p_DDFAVASummarized <- DDFAVASummarized %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_point(aes(x=totalAbundance, y=avgFAVA), alpha=0.6, size=0.75) +
  geom_errorbar(aes(x=totalAbundance, ymax=maxFAVA, ymin=minFAVA), width=0, alpha=0.6, size=0.5) +
  geom_errorbarh(aes(xmax=maxAbundance, xmin=minAbundance, y=avgFAVA), height=0, alpha=0.6, size=0.5) +
  #geom_text_repel(aes(x=totalAbundance, y=avgFAVA, label=combo), size=1.5, 
  #                box.padding=0.35, seed=0, nudge_x=-0.025, nudge_y=0.0005) +
  #geom_smooth(aes(x=totalAbundance, y=avgFAVA), method="lm", color="black", linetype="dotted", se=FALSE) +
  #geom_text(aes(x=totalAbundance, y=avgFAVA, label=combo), size=2, nudge_x=-0.04, nudge_y=0.0015) +
  xlab("Total DD relative abundance") +
  ylab("Variability across mix. ratios (FST)") +
  xlim(0,0.3) +
  ylim(0,0.05) +
  DEFAULTS.THEME_PRINT
p_DDFAVASummarized
save_plot("out/DDFAVAsummarized.pdf", p_DDFAVASummarized, base_width=3, base_height=2)

# Calculate correlation of FAVA with DD ASVs, average across replicates.
DDFAVAlm <- cor.test(DDFAVASummarized %>% filter(passage==5) %>% pull(totalAbundance),
                     DDFAVASummarized %>% filter(passage==5) %>% pull(avgFAVA), 
                     method = "pearson")

DDFAVAlmFiltered <- cor.test(DDFAVASummarized %>% filter(passage==5 & combo!="XCA-XDA") %>% pull(totalAbundance),
                             DDFAVASummarized %>% filter(passage==5 & combo!="XCA-XDA") %>% pull(avgFAVA), 
                             method = "pearson")

# Plot FAVA and total DD relative abundance of actual mixtures for each replicate separately.
DDFAVA <- DDRelAb %>% 
  mutate(passage = as.character(passage)) %>% 
  left_join(FAVA %>% filter(mixtureType=="actual"))

p_DDFAVA <- DDFAVA %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_point(aes(x=totalAbundance, y=FAVA)) +
  geom_smooth(aes(x=totalAbundance, y=FAVA), method="lm", color="black", linetype="dotted", se=FALSE) +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~line)
p_DDFAVA

# Calculate correlation of FAVA with DD ASVs for each replicate separately.
DDFAVAlm1 <- lm(DDFAVA %>% filter(line=="actual-1") %>% pull(FAVA) ~ 
                  DDFAVA %>% filter(line=="actual-1") %>% pull(totalAbundance))
DDFAVAlm2 <- lm(DDFAVA %>% filter(line=="actual-2") %>% pull(FAVA) ~ 
                  DDFAVA %>% filter(line=="actual-2") %>% pull(totalAbundance))
DDFAVAlm3 <- lm(DDFAVA %>% filter(line=="actual-3") %>% pull(FAVA) ~ 
                  DDFAVA %>% filter(line=="actual-3") %>% pull(totalAbundance))

# Compute FAVA in sliding windows.
FAVAwindows <- window_fava(dfFAVA, group="sample", K=1362, window_size=4, window_step=1)
FAVAwindowsDf <- FAVAwindows$window_data %>% 
  separate(group, into=c("line","combo","passage"), sep="_", remove=FALSE) %>% 
  mutate(mixtureType = sub("-.*","",line))

# Plot FAVA across windows for each mixture individually.
p_FAVAwindows <- FAVAwindowsDf %>% 
  filter(passage==5 & mixtureType=="actual") %>% 
  ggplot() +
  geom_point(aes(x=window_index, y=FAVA)) +
  geom_line(aes(x=window_index, y=FAVA, group=line)) +
  facet_wrap(~combo, nrow=2) +
  DEFAULTS.THEME_PRINT
p_FAVAwindows  

# Test built-in sliding window plot function.
FAVAwindows$window_plot +
  facet_wrap(~group) +
  theme(legend.position="none")

# Test built-in relative abundance plot function.
plot_relabund(dfFAVA, group="sample", K=1362)
