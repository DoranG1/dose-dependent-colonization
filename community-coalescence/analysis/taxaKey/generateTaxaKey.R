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
library(data.table)

# Import data and sample sheet. ---------------------------------------------

# Set working directory for running locally.
setwd("../")
plotdir <- "analysis/taxaKey"

# Set the plot theme.
theme_set(theme_cowplot())
source("config/palettes/plotDefaults.R")

# Import taxa color palette.
palette <- read.table("config/palettes/KCHcolors-Silva-partial.txt",
                      header=TRUE, stringsAsFactors = FALSE)

# Import read summary statistics.
outdir <- "workflow/out/e0017"
output_path <- file.path(outdir,"DADA2_output")
data <- read.table(file.path(output_path, "L5_summary.txt"), header=TRUE, stringsAsFactors=FALSE)
taxa <- colnames(data)
data <- rownames_to_column(data) %>% 
  dplyr::rename(sample=rowname) %>%
  separate(sample, into=c("filename","stem"), sep="_")

# Import sample sheet.
ss <- read_tsv("config/210819-e0015-e0020-16S-samplesheet.txt") %>%
  filter(group=="e0017") %>%
  select(-c(sample, group, plate, media, 
            cultureID, stockID, bioTyperDate, bioTyperClassification))
data <- left_join(data, ss, by=c("filename"))

# Tidy data.
data <- data %>% gather(key="taxa",value="count",taxa)

# Calculate relative abundance.
data <- data %>%
  group_by(filename) %>%
  mutate(relAbundance=count/sum(count),
         relAbundance=ifelse(is.na(relAbundance),0,relAbundance))

# Create short version of taxa names.
data <- data %>%
  mutate(taxaShort=gsub(".*\\.","",taxa))

# Bind colors.
data <- left_join(data, 
                  palette %>% dplyr::select(taxa, hex) %>% filter(taxa %in% data$taxa), 
                  by=c("taxa"))

# Extract color palette. Convert unknown taxa to dark gray.
taxaPalette <- data %>% ungroup() %>%
  dplyr::select(taxa, taxaShort, hex) %>%
  unique() %>% mutate(hex=ifelse(is.na(hex),"#615c5c", hex)) %>%
  arrange(taxa) %>% 
  mutate(taxaShort=fct_reorder(taxaShort, taxa))
taxaPaletteList <- list(taxaPalette$hex)[[1]]
names(taxaPaletteList) <- taxaPalette$taxa
# Create a version of the color palette with shortened taxa names.
taxaPalette <- taxaPalette %>%
  mutate(taxaShort=gsub(".*\\.","",taxa))
taxaShortPaletteList <- list(taxaPalette$hex)[[1]]
names(taxaShortPaletteList) <- taxaPalette$taxaShort

# Extract list of more prevalent taxa.
abundantTaxa <- (data %>% ungroup() %>% group_by(taxa) %>%
                   summarize(maxRelAbundance=max(relAbundance)) %>%
                   filter(maxRelAbundance>0.03))$taxa

# Get abundant families for example plot (XBA-029).
abundantTaxaExample <- (read_tsv("analysis/e0012plots/familiesAbundantExample.tsv", col_names=TRUE))$x

# Get abundant families from all e0012.
abundantTaxaAlle0012 <- fread("analysis/e0012plots/familiesAbundantAlle0012.txt")$x

# Plot key of abundant taxa.
taxaKeyAlle0012 <- taxaPalette %>%
  filter(taxa %in% abundantTaxaAlle0012) %>% 
  ggplot() +
  geom_tile(aes(x=0, y=factor(taxaShort), fill=taxaShort), color="black", size=0.05) +
  scale_fill_manual(values=taxaShortPaletteList) +
  scale_y_discrete(limits=
                     ((taxaPalette %>% 
                         filter(taxa %in% abundantTaxaAlle0012))$taxaShort)
                   [order(rev((taxaPalette %>% 
                                 filter(taxa %in% abundantTaxaAlle0012))$taxa))]) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(hjust=1, size=5.5))
taxaKeyAlle0012
save_plot(paste0(plotdir, "/taxaKeyAlle0012.pdf"), taxaKeyAlle0012, base_width=1.3, base_height=1.6)
  