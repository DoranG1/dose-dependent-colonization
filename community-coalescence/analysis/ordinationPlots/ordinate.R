library(dada2)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggfittext)
library(seqinr)
library(cowplot)
library(exactRankTests)
library(reshape2)
library(forcats)
library(foreach)

# Import data and sample sheet. ---------------------------------------------

# Set working directory for running locally.
setwd("/oak/stanford/groups/relman/users/dorang/210827-e0017-16S")
source("analysis/generateMixturesDataframe.R")
plotdir <- "analysis/ordinationPlots"

# Select just extraction replicate 1.
ps <- prune_samples(sample_data(ps)$extractionReplicate==1, ps)

# Prune samples without any reads.
ps <- prune_samples(sample_sums(ps)>=1, ps)

# Ordination --------------------------------------------------------------

# Preprocess the data.
# Remove OTUs that do not appear more than 5 times in more than half the samples.
wh0 <- genefilter_sample(ps, filterfun_sample(function(x) x>5), 
                         A=0.5*nsamples(ps))
ps1 <- prune_taxa(wh0, ps)
# Transform to even sampling depth.
ps1 <- transform_sample_counts(ps1, function(x) 1E6 * x/sum(x))
# Prune samples with fewer than 100 total reads.
ps1 <- prune_samples(sample_sums(ps1)>=100, ps1)

# Create a dataframe from certain sample data variables for easier manipulation.
filename <- sample_data(ps1)$filename
sample1 <- gsub("-.*","",sample_data(ps1)$biosample1)
sample2 <- gsub("-.*","",sample_data(ps1)$biosample2)
ratio <- sample_data(ps1)$ratio
metadata_df <- data.frame(filename, biosample1=sample1, biosample2=sample2, ratio)
# Dummy dataframes to record original sample order.
order_df <- data.frame(filename=sample_data(ps1)$filename)
controls_df <- metadata_df %>%
  filter(is.na(biosample2)) %>%
  mutate(combo=biosample1, label=biosample1)
  
# Add combo and labels to sample metadata
combos <- c("XBA-XBB", "XCA-XCB", "XDA-XDB", "XFA-XFB", "XBA-XCA", "XCA-XDA", "XDA-XFA", "XFA-XBA")
metadata_df <- foreach(x=combos, .combine="rbind") %do% {
  selectsample1 <- gsub("-.*", "", x)
  selectsample2 <- gsub(".*-", "", x)
  df <- metadata_df %>%
    filter(biosample1 %in% c(selectsample1, selectsample2) &
             biosample2 %in% c(selectsample1, selectsample2)) %>%
    mutate(combo=x) %>%
    mutate(ratio1=gsub(":.*", "", ratio), ratio2=gsub(".*:", "", ratio)) %>%
    mutate(label = case_when(
      biosample1==selectsample1 & !is.na(biosample2) ~ 
        paste(selectsample1, ratio, selectsample2, sep="-"),
      biosample1!=selectsample1 & !is.na(biosample2) ~ 
        paste(selectsample1, paste0(ratio2, ":", ratio1), selectsample2, sep="-"))) %>%
    select(-c(ratio1, ratio2))
}

# Bind labeled combo df to controls and rearrange to original order, add back to sample metadata.
metadata_df <- rbind(metadata_df, controls_df)
metadata_df <- left_join(order_df, metadata_df)
sample_data(ps1)$combo <- metadata_df$combo
sample_data(ps1)$label <- metadata_df$label

# Filtered versions of dataframe for each passage.
ps1_P3 <- prune_samples(sample_data(ps1)$passage==3, ps1)
ps1_P5 <- prune_samples(sample_data(ps1)$passage==5, ps1)

# Ordinate.
ps1_ord <- ordinate(ps1, method="NMDS")
ps1_ord_P3 <- ordinate(ps1_P3, method="NMDS")
ps1_ord_P5 <- ordinate(ps1_P5, method="NMDS")

# Plot ordination.

p_ps1_ord <- plot_ordination(ps1, ps1_ord, color="combo")
save_plot(paste0(plotdir, "/allMixturesOrdination.png"), p_ps1_ord)

p_ps1_P3_ord <- plot_ordination(ps1_P3, ps1_ord_P3, color="combo")
save_plot(paste0(plotdir, "/allMixturesOrdination_P3.png"), p_ps1_P3_ord)

p_ps1_P5_ord <- plot_ordination(ps1_P5, ps1_ord_P5, color="combo")
save_plot(paste0(plotdir, "/allMixturesOrdination_P5.png"), p_ps1_P5_ord)

# Split ordination into each set of controls/combo to see more clearly.
foreach(x=combos) %do% {
  sample1 <- gsub("-.*","",x)
  sample2 <- gsub(".*-","",x)
  # Plot both passages.
  psCombo <- subset_samples(ps1,sample_data(ps1)$combo %in% c(x, sample1, sample2))
  ordCombo <- ordinate(psCombo, method="NMDS")
  ordComboPlot <- plot_ordination(psCombo, ordCombo, color="combo")
  plotName <- paste0(x, "-NMDS")
  save_plot(paste0(plotdir, "/", plotName, ".png"), ordComboPlot)
  # Plot passage 3.
  psCombo_P3 <- subset_samples(ps1_P3,sample_data(ps1_P3)$combo %in% c(x, sample1, sample2))
  ordCombo_P3 <- ordinate(psCombo_P3, method="NMDS")
  ordCombo_P3_Plot <- plot_ordination(psCombo_P3, ordCombo_P3, color="combo")
  plotName <- paste0(x, "-NMDS_P3")
  save_plot(paste0(plotdir, "/", plotName, ".png"), ordCombo_P3_Plot)
  # Plot passage 5.
  psCombo_P5 <- subset_samples(ps1_P5,sample_data(ps1_P5)$combo %in% c(x, sample1, sample2))
  ordCombo_P5 <- ordinate(psCombo_P5, method="NMDS")
  ordCombo_P5_Plot <- plot_ordination(psCombo_P5, ordCombo_P5, color="combo")
  plotName <- paste0(x, "-NMDS_P5")
  save_plot(paste0(plotdir, "/", plotName, ".png"), ordCombo_P5_Plot)
}
