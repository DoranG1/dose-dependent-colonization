library(dada2); packageVersion("dada2")
library(tidyverse)

# EDIT THE FOLLOWING FOUR LINES
outdir <- snakemake@params[["outdir"]]
per_sample_folder_fp <- snakemake@params[["trimdir"]] #THIS SHOULD BE THE FOLDER CONTAINING THE FWD AND REV_READS FOLDERS
truncLenF <- 200
truncLenR <- 130

fwd_path <- file.path(per_sample_folder_fp,"R1")
rev_path <- file.path(per_sample_folder_fp,"R2")

fwd_filt_path <- file.path(outdir,"filtered/R1")
rev_filt_path <- file.path(outdir,"filtered/R2")

fwd_fastq <- sort(list.files(fwd_path, pattern="fastq.gz"))
rev_fastq <- sort(list.files(rev_path, pattern="fastq.gz"))


out <- filterAndTrim(fwd=file.path(fwd_path, fwd_fastq), filt=file.path(fwd_filt_path, fwd_fastq), 
                     rev=file.path(rev_path, rev_fastq), filt.rev=file.path(rev_filt_path, rev_fastq), 
                     truncLen=c(truncLenF, truncLenR), maxEE=c(2,2), truncQ=2,
                     maxN=0, rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=FALSE)

# Export summary of reads analyzed and filtered.
output_path <- file.path(outdir,"DADA2_output")
if(!file_test("-d", output_path)) dir.create(output_path)
write.table(as.data.frame(out) %>% rownames_to_column("sample"),
            file.path(output_path,"read_summary.txt"), row.names=FALSE, quote=FALSE)

filtpathF <- fwd_filt_path
filtpathR <- rev_filt_path

filtFs <- list.files(filtpathF, pattern=".fastq", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern=".fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "\\."), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "\\."), `[`, 1)

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

# Learn forward error rates
errF <- learnErrors(filtFs, nread=2e6, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nread=2e6, multithread=TRUE)


# plot error rates
plotErrors(errF, nominalQ=TRUE)
ggsave(filename = paste(output_path, "error_plot_F.pdf", sep = "/"))
plotErrors(errR, nominalQ=TRUE)
ggsave(filename = paste(output_path, "error_plot_R.pdf", sep = "/"))

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))

names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file.path(output_path,"seqtab.rds"))

seqtab_nobim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
write.table(seqtab_nobim, file.path(output_path,'all_dsvs.no_bimera.txt'), sep='\t', row.names=TRUE, col.names=TRUE)
saveRDS(seqtab_nobim, file.path(output_path,"seqtab_all.rds"))

