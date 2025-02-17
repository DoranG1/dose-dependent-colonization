library(phyloseq)
library(ape)
library(tidyverse)
library(cowplot)
library(foreach)

# Set output directories.
outdir <- snakemake@params[["outdir"]]
output_path <- file.path(outdir,"DADA2_output")
if(!file_test("-d", output_path)) dir.create(output_path)

# Set the plot theme.
theme_set(theme_cowplot())
source("../../config/palettes/plotDefaults.R")
plotdir <- "workflow/out/plots/ASVComposition"


# Generate full phyloseq object -------------------------------------------

if(!file_test("-f",file.path(output_path,"ps_all.rds"))){
  # Import phyloseq object with OTU table and taxonomy table.
  # The OTU table is generated after filtering out bimeras.
  ps_taxa <- readRDS(file.path(output_path,"ps_taxa.rds"))
  
  # Import the tree and add it to the phyloseq object.
  tree <- read.tree(file.path(output_path,'dsvs_msa.tree'))
  
  # Import sample metadata.
  meta <- read.table("../../config/20210201-e0012-samplesheet.txt",
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
} else {
  ps <- readRDS(file.path(output_path,"ps_all.rds"))
}


# Plot alpha diversity for each sample ------------------------------------

# Calculate alpha diversity for each sample.
alpha <- estimate_richness(ps)
# Append the alpha diversity to the sample metadata.
alpha <- cbind(sample_data(ps), alpha)

# Arrange the media for plotting.
mediaOrder <- c("inoculum","BHI","mBHI","mBHImucin","TYG")

# Plot the species richness over time for each media and community.
save_plot(file.path(plotdir,"speciesRichness-byPassage.png"),
          alpha %>%
            filter(PCRreplicate==1) %>%
            ggplot() +
            geom_point(aes(x=passage, y=Observed, color=fct_relevel(media, levels=mediaOrder))) +
            geom_line(aes(x=passage, y=Observed, color=fct_relevel(media, levels=mediaOrder), 
                          group=interaction(media,inoculationreplicate))) +
            facet_wrap(~biosample, scales="free", ncol=4) +
            scale_color_brewer(name="Media", palette="Set1") +
            xlab("Passage") + ylab("Number of ASVs") +
            ylim(0,max(alpha$Observed)+25) +
            DEFAULTS.THEME_ALL, ncol=1.5)
# Plot the species richness over time for a single set of replicate communities.
save_plot(file.path(plotdir,"speciesRichness-byPassage-XBA.png"),
          alpha %>%
            filter(PCRreplicate==1, biosample=="XBA-029", 
                   media %in% c("inoculum","mBHImucin")) %>%
            ggplot() +
            geom_point(aes(x=passage, y=Observed, color=fct_relevel(media, levels=mediaOrder))) +
            geom_line(aes(x=passage, y=Observed, color=fct_relevel(media, levels=mediaOrder), 
                          group=interaction(media,inoculationreplicate))) +
            scale_color_manual(name="", values=c("firebrick3","black")) +
            xlab("Passage") + ylab("Number of ASVs") +
            ylim(0,150) + theme(legend.position="top") +
            DEFAULTS.THEME_ALL, ncol=0.4, nrow=0.5)

# Plot the species richness at the final timepoint for each biological sample.
save_plot(file.path(plotdir,"speciesRichness-final-bySample.png"),
          alpha %>%
            filter(PCRreplicate==1, passage %in% c(0,15)) %>%
            ggplot() +
            geom_point(aes(x=biosample, y=Observed, 
                           color=fct_relevel(media, levels=mediaOrder),
                           group=fct_relevel(media, levels=mediaOrder))) +
            scale_color_brewer(name="Media", palette="Set1") +
            xlab("Sample") + ylab("Number of ASVs, passage 15") +
            theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
            DEFAULTS.THEME_ALL, ncol=0.8)


# Plot the species richness at the final timepoint for each media condition.
save_plot(file.path(plotdir,"speciesRichness-final-byMedia.png"),
          alpha %>%
            filter(PCRreplicate==1, passage %in% c(0,15)) %>%
            ggplot() +
            geom_point(aes(x=media, y=Observed, 
                           color=factor(biosample))) +
            scale_color_brewer(name="Sample", palette="Paired") +
            xlab("Media") + ylab("Number of ASVs, passage 15") +
            theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
            DEFAULTS.THEME_ALL, ncol=0.8)


# Plot ordinations --------------------------------------------------------

# Preprocess the data.
# Remove OTUs that do not appear more than 5 times in more than half the samples.
wh0 <- genefilter_sample(ps, filterfun_sample(function(x) x>5), 
                        A=0.5*nsamples(ps))
ps1 <- prune_taxa(wh0, ps)
# Transform to even sampling depth.
ps1 <- transform_sample_counts(ps1, function(x) 1E6 * x/sum(x))
# Prune samples with fewer than 100 total reads.
ps1 <- prune_samples(sample_sums(ps1)>=100, ps1)


# Ordinate the data using a variety of distance metrics.
distances <- c("jaccard","bray","jsd","unifrac","wunifrac")
ords <- foreach(x=distances, .combine="c") %do% {
  # Ordinate the data.
  ps1.ord <- ordinate(ps1,"NMDS",x)
  # Plot the ordination.
  save_plot(paste0(plotdir,"/ordination/",x,"-all.png"),
            plot_ordination(ps1, ps1.ord, type="samples", 
                            color="biosample", shape="media") +
              scale_color_brewer(name="Sample", palette="Paired") +
              ggtitle(x) +
              DEFAULTS.THEME_ALL)
  # Plot the ordination, splitting panels by timepoint and media.
  save_plot(paste0(plotdir,"/ordination/",x,"-split.png"),
            plot_ordination(ps1, ps1.ord, type="samples", 
                            color="biosample") +
              geom_polygon(aes(fill=biosample), alpha=0.2) +
              facet_grid(fct_relevel(media, levels=mediaOrder)~passage) +
              scale_color_brewer(name="Sample", palette="Paired") +
              scale_fill_brewer(name="Sample", palette="Paired") +
              ggtitle(x) +
              DEFAULTS.THEME_ALL, ncol=1.3, nrow=1.3)
  ps1.ord
}
