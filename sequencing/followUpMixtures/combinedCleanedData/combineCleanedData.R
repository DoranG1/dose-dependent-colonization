library(data.table)
library(tidyverse)
library(cowplot)
library(forcats)
library(R.utils)
library(foreach)
library(scales)

data4C <- fread("../4C/data/intermediates/e0063MixtureDataframe.txt") %>% 
  select(!plate)
data4D <- fread("../4D/data/intermediates/e0065MixtureDataframe.txt") %>% 
  select(!plate)
data4F <- fread("../4F-5D/data/intermediates/e0043MixtureDataframe.txt") %>% 
  mutate(combo = ifelse(combo=="D1-Strep17","XFA-Strep17",combo))
data4E <- fread("../4E-4G-5C-5E-5F-5G/data/intermediates/e0049e0050MixtureDataframe.txt")
dataS16A <- fread("../S16A/data/intermediates/XFA-Strep17-p3-mixtureDataframe.txt")

dataCombined <- rbind(data4C, data4D, data4F, data4E, dataS16A)

write.table(dataCombined, "combinedCleanedMixtureDataframe.txt", quote=FALSE, row.names=FALSE, sep="\t")
