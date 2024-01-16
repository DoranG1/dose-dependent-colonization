library(dada2)
library(msa)
library(phangorn)
library(phyloseq)

#CHANGE THE FOLLOWING LINE
outdir <- snakemake@params[["outdir"]]
output_path <- file.path(outdir,"DADA2_output") #THIS SHOULD BE THE DADA2_output folder
seqtab_nobim <- readRDS( file.path(output_path,"seqtab_all.rds"))

seqs <- getSequences(seqtab_nobim)
names(seqs) <- seqs

mult <- msa(seqs, method="ClustalW", type="dna", order="input")
saveRDS(mult, file.path(output_path,"mult.rds"))

phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
saveRDS(fit, file.path(output_path,"fit.rds"))

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
saveRDS(fitGTR, file.path(output_path,"mult.rds"))

write.tree(phyloseq(phy_tree(fitGTR$tree)),file.path(output_path,'dsvs_msa.tree'))
