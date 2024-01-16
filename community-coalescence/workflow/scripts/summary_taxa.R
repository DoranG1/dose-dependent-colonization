library(dada2)
library(phyloseq)

#CHANGE THE NEXT LINES
outdir <- snakemake@params[["outdir"]]
output_path <- file.path(outdir,"DADA2_output") #THIS SHOULD BE THE DADA2_output folder
seqtab_nobim <- readRDS( file.path(output_path,"seqtab_all.rds"))
ref_fasta <- "workflow/out/Sil138_BAhq_V4_train_set.fa.gz"
taxa <- assignTaxonomy(seqtab_nobim,refFasta = ref_fasta,multithread=FALSE)

tax_tab <- tax_table(taxa)

pseq <- phyloseq(otu_table(seqtab_nobim,taxa_are_rows = F),tax_tab)
saveRDS(pseq, file.path(output_path,"ps_taxa.rds"))

write.table(tax_tab,  file.path(output_path,"tax_table.txt"), sep='\t', row.names=TRUE, col.names=TRUE)

otut<-otu_table(tax_glom(pseq,taxrank = "Phylum"))
taxt<-tax_table(tax_glom(pseq,taxrank = "Phylum"))
write.table(taxt,  file.path(output_path,"L2_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
a <- read.table(file.path(output_path,"L2_taxa.txt"))

a1<-as.character(a[,1]) 
a2<-as.character(a[,2])
tax_ph<- paste(a1,a2,sep=";")
write.table(otut, file.path(output_path,"L2_summary.txt"), sep='\t', row.names=TRUE, col.names=tax_ph)



otut<-otu_table(tax_glom(pseq,taxrank = "Class"))
taxt<-tax_table(tax_glom(pseq,taxrank = "Class"))
write.table(taxt,  file.path(output_path,"L3_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
a <- read.table(file.path(output_path,"L3_taxa.txt"))

a1<-as.character(a[,1]) 
a2<-as.character(a[,2])
a3<-as.character(a[,3])
tax_ph<- paste(a1,a2,a3,sep=";")
write.table(otut, file.path(output_path,"L3_summary.txt"), sep='\t', row.names=TRUE, col.names=tax_ph)



otut<-otu_table(tax_glom(pseq,taxrank = "Order"))
taxt<-tax_table(tax_glom(pseq,taxrank = "Order"))
write.table(taxt,  file.path(output_path,"L4_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
a <- read.table(file.path(output_path,"L4_taxa.txt"))

a1<-as.character(a[,1]) 
a2<-as.character(a[,2])
a3<-as.character(a[,3])
a4<-as.character(a[,4])
tax_ph<- paste(a1,a2,a3,a4,sep=";")
write.table(otut, file.path(output_path,"L4_summary.txt"), sep='\t', row.names=TRUE, col.names=tax_ph)



otut<-otu_table(tax_glom(pseq,taxrank = "Family"))
taxt<-tax_table(tax_glom(pseq,taxrank = "Family"))
write.table(taxt,  file.path(output_path,"L5_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
a <- read.table(file.path(output_path,"L5_taxa.txt"))

a1<-as.character(a[,1]) 
a2<-as.character(a[,2])
a3<-as.character(a[,3])
a4<-as.character(a[,4])
a5<-as.character(a[,5])
tax_ph<- paste(a1,a2,a3,a4,a5,sep=";")
write.table(otut, file.path(output_path,"L5_summary.txt"), sep='\t', row.names=TRUE, col.names=tax_ph)



otut<-otu_table(tax_glom(pseq,taxrank = "Genus"))
taxt<-tax_table(tax_glom(pseq,taxrank = "Genus"))
write.table(taxt,  file.path(output_path,"L6_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
a <- read.table(file.path(output_path,"L6_taxa.txt"))

a1<-as.character(a[,1]) 
a2<-as.character(a[,2])
a3<-as.character(a[,3])
a4<-as.character(a[,4])
a5<-as.character(a[,5])
a6<-as.character(a[,6])
tax_ph<- paste(a1,a2,a3,a4,a5,a6,sep=";")
write.table(otut, file.path(output_path,"L6_summary.txt"), sep='\t', row.names=TRUE, col.names=tax_ph)
