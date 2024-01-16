## Script to run on cluster for permutations.

library(phyloseq)
library(tidyverse)
library(foreach)
library(ape)
library(doParallel)
library(data.table)

registerDoParallel(cores=16)

# Read in phylogenetic tree.
tree <- read.tree("../../workflow/out/e0017/DADA2_output/dsvs_msa.tree")

# Create phylogenetic distance matrix between ASVs from tree.
distMatrix <- cophenetic.phylo(tree)

# Turn distance matrix into dataframe.
distDf <- data.frame(distMatrix) %>% 
  mutate(OTU1 = row.names(distMatrix)) %>% 
  pivot_longer(cols=!OTU1, names_to="OTU2", values_to="distance")

combos <- c("XBA-XBB", "XCA-XCB", "XDA-XDB", "XFA-XFB", "XBA-XCA", "XCA-XDA", "XDA-XFA", "XFA-XBA")
passages <- c("3","5")
comboPassages <- expand.grid(combos,passages) %>% mutate(comboPassage=paste(Var1,Var2,sep="_"))
comboPassages <- comboPassages$comboPassage

theoretical_mixtures_df_categorized <- fread("../mixtureDataframe.txt")

# Analyze closest neighbors for undercolonizing ASVs. ----------------------

# Permute UC ASVs and calculate null distribution of closest phylogenetic distances.
UCASVDistancesPermuted <- foreach(x=comboPassages, .combine="rbind") %dopar% {
  currCombo <- sub("_.*","",x)
  currPassage <- sub(".*_","",x)
  currASVsDf <- theoretical_mixtures_df_categorized %>% 
    filter(passage==currPassage & combo==currCombo & category!="lowAbundance") %>% 
    ungroup() %>% 
    select(OTU, Family, OTUnum, category) %>% 
    unique()
  numASVs <- n_distinct(currASVsDf %>% filter(category=="undercolonizing") %>% pull(OTU))
  if (numASVs>1) {
    UCASVDistancesPermuted <- foreach(y=seq(1:1000), .combine="rbind") %dopar% {
      UCASVsDf <- currASVsDf %>% 
        slice_sample(n=numASVs)
      allCurrUCASVs <- UCASVsDf %>% pull(OTU) %>% unique()
      distanceDf <- foreach(z=allCurrUCASVs, .combine="rbind") %dopar% {
        currFamily <- UCASVsDf %>% filter(OTU==z) %>% pull(Family)
        currOTUnum <- UCASVsDf %>% filter(OTU==z) %>% pull(OTUnum)
        currDf <- UCASVsDf %>% 
          filter(OTU!=z)
        currUCASVs <- unique(currDf$OTU)
        distanceDf <- distDf %>% 
          filter(OTU1==z & OTU2 %in% currUCASVs) %>% 
          left_join(currDf %>% rename(OTU2=OTU, Family2=Family, OTUnum2=OTUnum) %>% select(!category)) %>% 
          mutate(combo=currCombo, passage=currPassage, Family1=currFamily, OTUnum1=currOTUnum)
      }
      distanceDf <- distanceDf %>% mutate(permutationReplicate=y)
    }
  }
}

write.table(UCASVDistancesPermuted, "UCASVDistancesPermuted.txt", quote=FALSE, row.names=FALSE, sep="\t")