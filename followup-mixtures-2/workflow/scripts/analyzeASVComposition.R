library(phyloseq)
library(ape)
library(tidyverse)
library(cowplot)
library(foreach)

# List snakemake parameters required for running script.
# snakemake@params[["outdir"]] Path to output directory for plots and small output files
# snakemake@params[["metadata"]] Path to sample metadata
OUTDIR <- snakemake@params[["outdir"]]
ifelse(!dir.exists(OUTDIR),dir.create(OUTDIR, recursive=TRUE),FALSE)
METADATA <- snakemake@params[["metadata"]]

# Set output directories.
outdir <- OUTDIR
output_path <- file.path(outdir,"DADA2_output")
if(!file_test("-d", output_path)) dir.create(output_path)


# Generate full phyloseq object -------------------------------------------

# Import phyloseq object with OTU table and taxonomy table.
# The OTU table is generated after filtering out bimeras.
ps_taxa <- readRDS(file.path(output_path,"ps_taxa.rds"))

# Import the tree and add it to the phyloseq object.
tree <- read.tree(file.path(output_path,'dsvs_msa.tree'))

# Import sample metadata.
meta <- read.table(METADATA,
                   header=TRUE, stringsAsFactors = FALSE)
meta <- meta %>%
  mutate(fullSample=paste0(filename,"_",sample))
# Filter to include only samples that were retained through dada2.
meta <- meta %>%
  filter(fullSample %in% rownames(otu_table(ps_taxa)))
row.names(meta) <- meta$fullSample

# Merge the OTU table, taxonomy table, sample data, and tree.
ps <- merge_phyloseq(ps_taxa, tree)
sample_data(ps) <- meta

# Export the complete phyloseq object.
saveRDS(ps, file.path(output_path,"ps_all.rds"))

