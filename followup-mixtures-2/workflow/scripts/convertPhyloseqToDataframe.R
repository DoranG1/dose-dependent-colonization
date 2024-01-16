library(phyloseq)
library(tidyverse)
library(ape)

# List snakemake parameters required for running script.
# snakemake@params[["outdir"]] Path to output directory for plots and small output files
OUTDIR <- snakemake@params[["outdir"]]
ifelse(!dir.exists(OUTDIR),dir.create(OUTDIR, recursive=TRUE),FALSE)

# Set output directories.
outdir <- OUTDIR
output_path <- file.path(outdir,"DADA2_output")
if(!file_test("-d", output_path)) dir.create(output_path)

# Import full phyloseq object with OTU table, taxonomy table, tree, and metadata.
ps <- readRDS(file.path(output_path,"ps_all.rds"))

# Extract each component of the phyloseq object.
# Convert OTU table to data frame.
# Method described here: https://github.com/joey711/phyloseq/issues/613
OTU1 <- as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
dataOTU <- as.data.frame(OTU1)
# Tidy the data frame.
dataOTU$sample <- rownames(otu_table(ps))
dataOTUtidy <- dataOTU %>%
  pivot_longer(-sample, names_to="OTU", values_to="count")

# Convert samplesheet to data frame.
sampleData1 <- as(sample_data(ps), "matrix")
if(taxa_are_rows(ps)){sampleData1 <- t(sampleData1)}
# Coerce to data.frame
dataSamplesheet <- as.data.frame(sampleData1)
dataSamplesheet$sample <- rownames(sample_data(ps))

# Merge the OTU counts with the sample metadata.
data <- left_join(dataOTUtidy, dataSamplesheet, by=c("sample"))
rm(dataOTUtidy)

# Filter out OTUs that are not observed in a particular sample.
data <- data %>%
  filter(count!=0)

# Calculate the relative abundance of each OTU.
data <- data %>% group_by(sample) %>%
  mutate(relAbundance=count/sum(count))

# Merge the OTU taxonomy with the OTU counts and samplesheet.
# Extract the OTU taxonomy.
dataTaxonomy <- as.data.frame(tax_table(ps))
dataTaxonomy$OTU <- row.names(dataTaxonomy)

# Annotate the taxonomy of each OTU with information about its abundance.
data <- left_join(data, dataTaxonomy, by=c("OTU"))

# Export all data.
write.table(data, gzfile(file.path(output_path,"ps_all.txt.gz")),
            row.names=FALSE, quote=FALSE)
