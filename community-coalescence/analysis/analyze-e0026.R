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

# Set working directory for running locally.
setwd("../")

# Import read data.
ps <- readRDS("config/e0026-ps_all.rds")
# Import samplesheet.
ss <- read_tsv("config/220118-e0023-e0026-e0029-e0030-16S-samplesheet.txt")

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
  filter(count>0 & experiment=="e0026") %>% 
  group_by(filename) %>% 
  mutate(rel_abundance=count/sum(count))

# Check OTU presence in XFA-XFB.
# OTU #7 from XFA.
dfXFA <- df %>% 
  filter(biosample1 %in% c("XFA-029", "XFB-029") &
           OTU=="TACGTAGGTCCCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGTGGTTTCTTAAGTCTGATGTAAAAGGCAGTGGCTCAACCATTGTGTGCATTGGAAACTGGGAGACTTGAGTGCAGGAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGAGGCGAAAGCGGCTCTCTGGCCTGTAACTGACACTGAGGCTCGAAAGCGTGGGGAGCAAACAGG")

# OTU #17 from XFB.
dfXFB <- df %>% 
  filter(biosample1 %in% c("XFB-029", "XFA-029") &
           OTU=="TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGTGGTTTATTAAGTCTGGTGTAAAAGGCAGTGGCTCAACCATTGTATGCATTGGAAACTGGTAGACTTGAGTGCAGGAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCCTGTAACTGACACTGAGGCTCGAAAGCGTGGGGAGCAAACAGG")
