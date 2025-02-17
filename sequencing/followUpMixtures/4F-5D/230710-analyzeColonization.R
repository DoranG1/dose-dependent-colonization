library(data.table)
library(tidyverse)
library(cowplot)
library(forcats)
library(R.utils)
library(foreach)
library(data.table)
library(scales)

theme_set(theme_cowplot())
source("../../../config/palettes/plotDefaults.R")

# Import taxa color palette.
palette <- read.table("../../../config/palettes/KCHcolors-Silva-partial.txt",
                      header=TRUE, stringsAsFactors = FALSE) %>% 
  mutate(taxashort=gsub(".*\\.","",taxa))

# Import coverage data.
coverage <- fread("../../../config/read_summary.txt")

OTU7 <- "TACGTAGGTCCCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGTGGTTTCTTAAGTCTGATGTAAAAGGCAGTGGCTCAACCATTGTGTGCATTGGAAACTGGGAGACTTGAGTGCAGGAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGAGGCGAAAGCGGCTCTCTGGCCTGTAACTGACACTGAGGCTCGAAAGCGTGGGGAGCAAACAGG"
OTU17 <- "TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGTGGTTTATTAAGTCTGGTGTAAAAGGCAGTGGCTCAACCATTGTATGCATTGGAAACTGGTAGACTTGAGTGCAGGAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCCTGTAACTGACACTGAGGCTCGAAAGCGTGGGGAGCAAACAGG"

ratios <- c("1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000")
ratiosFull <- c("1:0", "1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000", "0:1")

# Read in data and add sample variables.
data <- fread("../../../data/ps_all.txt.gz") %>% 
  mutate(row=substr(well, 1, 1), col=as.numeric(substr(well, 2, length(well))),
         communityType=ifelse(row=="A", "parent", "mixture"),
         community=case_when(communityType=="parent" & col %in% c(1,2,3) ~ "7",
                             communityType=="parent" & col %in% c(4,5,6) ~ "17",
                             communityType=="parent" & col %in% c(7,8,9) ~ "XFA",
                             communityType=="parent" & col %in% c(10,11,12) ~ "XFB",
                             communityType=="mixture" & col %in% c(1,2,3) ~ "blank",
                             communityType=="mixture" & col %in% c(4,5,6) ~ "7-17",
                             communityType=="mixture" & col %in% c(7,8,9) ~ "XFA-17",
                             communityType=="mixture" & col %in% c(10,11,12) ~ "7-XFB"),
         ratio=case_when(row=="B" ~ "1000:1",
                         row=="C" ~ "100:1",
                         row=="D" ~ "10:1",
                         row=="E" ~ "1:1",
                         row=="F" ~ "1:10",
                         row=="G" ~ "1:100",
                         row=="H" ~ "1:1000"),
         logAbundance = log10(relAbundance))

dataCoverage <- left_join(data,
                          coverage %>% 
                            mutate(fullSample = sub(".fastq.gz","",sample)) %>% 
                            select(!sample)) %>%
  filter(community!="blank") %>% 
  group_by(filename) %>% 
  slice(1)

# Calculate coverage stats for paper.
sum(dataCoverage %>% filter(!community %in% c("XFB","7-XFB")) %>% pull(reads.out))

# Reformat parent communities to plot with mixtures.
dataPlotting <- data %>% 
  select(-c(sample, filename, round2plate, round1index, well, plate, group, fullSample, row)) 
dataPlotting <- dataPlotting %>% 
  mutate(combo=ifelse(communityType=="mixture", community, NA)) %>% 
  filter(!is.na(combo)) %>% 
  rbind(dataPlotting %>% filter(community=="7") %>% mutate(combo="7-17", ratio="1:0"), 
        dataPlotting %>% filter(community=="17") %>% mutate(combo="7-17", ratio="0:1"),
        dataPlotting %>% filter(community=="XFA") %>% mutate(combo="XFA-17", ratio="1:0"),
        dataPlotting %>% filter(community=="17") %>% mutate(combo="XFA-17", ratio="0:1"),
        dataPlotting %>% filter(community=="7") %>% mutate(combo="7-XFB", ratio="1:0"),
        dataPlotting %>% filter(community=="XFB") %>% mutate(combo="7-XFB", ratio="0:1")) %>% 
  mutate(rep=ifelse(col%%3!=0, col%%3, 3),
         relAbundanceOld = relAbundance,
         relAbundance = ifelse(relAbundance >= 10^-3, relAbundance, 10^-3)) %>% 
  group_by(combo) %>%
  complete(nesting(OTU, Family),
           ratio=ratiosFull, rep=c(1,2,3),
           fill=list(logAbundance=-3, relAbundance=10^-3)) %>% 
  mutate(combo=case_when(
    "7-17" %in% combo ~ "7-17",
    "XFA-17" %in% combo ~ "XFA-17",
    "7-XFB" %in% combo ~ "7-XFB",
    "blank" %in% combo ~ "blank"),
    ratio=fct_relevel(ratio, ratiosFull)) %>% 
  # Remove poorly sequenced replicates.
  #filter(!(combo=="7-XFB" & rep==3) & !(combo=="XFA-17" & rep==3)) %>% 
  # Bind taxa colors by family.
  left_join(palette %>% select(taxashort, taxa, hex) %>%
              filter(taxashort %in% dataPlotting$Family) %>% 
              rename(Family=taxashort) %>% 
              group_by(Family) %>% 
              mutate(hexnum=row_number()) %>% 
              filter(hexnum==1) %>% 
              select(!hexnum),
            by="Family")

# Extract color palette. Convert unknown taxa to dark gray.
taxaPalette <- dataPlotting %>% ungroup() %>%
  select(Family, taxa, hex) %>%
  unique() %>%
  mutate(hex=ifelse(is.na(hex),"#615c5c", hex)) %>% 
  arrange(taxa)
taxaPaletteList <- list(taxaPalette$hex)[[1]]
names(taxaPaletteList) <- taxaPalette$taxa

# Bind OTU numbers from e0017.
datae0017 <- fread("../../../data/mixtureDataframe.txt")
dataPlotting <- dataPlotting %>% 
  left_join(fread("../../../data/mixtureDataframe.txt") %>% 
              ungroup() %>% 
              select(Family, OTU, OTUnum) %>% 
              unique(),
            by=c("Family","OTU"))

# Get Entero ASVs.
EnteC2 <- (dataPlotting %>% filter(Family=="Enterococcaceae" & OTUnum==2))$OTU[1]
EnteC3 <- (dataPlotting %>% filter(Family=="Enterococcaceae" & OTUnum==3))$OTU[1]

# Remove contamination from 7-17 mixtures.
dataMixture1 <- dataPlotting %>% 
  filter(combo=="7-17" & OTU %in% c(OTU7,OTU17)) %>%
  group_by(ratio, rep) %>%
  summarize(OTU, relAbundanceOld = relAbundanceOld/sum(relAbundanceOld, na.rm=TRUE))

# Add dose difference stat to dataframe.
dataStats <- dataPlotting %>% 
  filter(ratio %in% ratios) %>% 
  #filter(!(combo=="XFA-17" & ratio=="1:1000")) %>%
  filter(!(combo=="XFA-17" & rep==3)) %>% 
  # Add dose distance/difference statistics.
  group_by(combo, OTU, rep) %>% 
  arrange(combo, OTU, rep, ratio) %>% 
  mutate(doseDist=abs(logAbundance-lag(logAbundance)),
         doseDistTotal=mean(doseDist, na.rm=TRUE),
         doseDiff=logAbundance-lag(logAbundance),
         doseDiffTotal=mean(doseDiff, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(OTU, combo) %>% 
  mutate(doseDistance=mean(doseDistTotal), 
         doseDifference=mean(doseDiffTotal))

# Generate theoretical curves. --------------------------------------------

# Theoretical curve for whole mixtures.
dataPlottingTheoretical <- dataPlotting %>%
  mutate(mixtureType = "actual") %>% 
  rbind(dataPlotting %>%
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(combo, rep, OTU) %>%
          filter(ratio %in% c("1:0","1000:1","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/1000) %>% 
          filter(ratio=="1000:1") %>% 
          ungroup() %>% 
          group_by(combo, rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(combo, rep, OTU) %>%
          filter(ratio %in% c("1:0","100:1","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/100) %>% 
          filter(ratio=="100:1") %>% 
          ungroup() %>% 
          group_by(combo, rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(combo, rep, OTU) %>%
          filter(ratio %in% c("1:0","10:1","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/10) %>% 
          filter(ratio=="10:1") %>% 
          ungroup() %>% 
          group_by(combo, rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(combo, rep, OTU) %>%
          filter(ratio %in% c("1:0","1:1","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]) %>% 
          filter(ratio=="1:1") %>% 
          ungroup() %>% 
          group_by(combo, rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(combo, rep, OTU) %>%
          filter(ratio %in% c("1:0","1:10","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/10 + relAbundanceOld[ratio=="0:1"]) %>% 
          filter(ratio=="1:10") %>% 
          ungroup() %>% 
          group_by(combo, rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(combo, rep, OTU) %>%
          filter(ratio %in% c("1:0","1:100","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/100 + relAbundanceOld[ratio=="0:1"]) %>% 
          filter(ratio=="1:100") %>% 
          ungroup() %>% 
          group_by(combo, rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(combo, rep, OTU) %>%
          filter(ratio %in% c("1:0","1:1000","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/1000 + relAbundanceOld[ratio=="0:1"]) %>% 
          filter(ratio=="1:1000") %>% 
          ungroup() %>% 
          group_by(combo, rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
  ) %>% 
  mutate(relAbundance = ifelse(is.na(relAbundanceOld) | relAbundanceOld<10^-3, 10^-3, relAbundanceOld),
         logAbundance = log10(relAbundance))
  
# Theoretical curve for mix 1 with contaminants removed.
dataMixture1Theoretical <- dataMixture1 %>%
  mutate(mixtureType = "actual") %>% 
  rbind(dataMixture1 %>%
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(rep, OTU) %>%
          filter(ratio %in% c("1:0","1000:1","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/1000) %>% 
          filter(ratio=="1000:1") %>% 
          ungroup() %>% 
          group_by(rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataMixture1 %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(rep, OTU) %>%
          filter(ratio %in% c("1:0","100:1","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/100) %>% 
          filter(ratio=="100:1") %>% 
          ungroup() %>% 
          group_by(rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataMixture1 %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(rep, OTU) %>%
          filter(ratio %in% c("1:0","10:1","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/10) %>% 
          filter(ratio=="10:1") %>% 
          ungroup() %>% 
          group_by(rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataMixture1 %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(rep, OTU) %>%
          filter(ratio %in% c("1:0","1:1","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]) %>% 
          filter(ratio=="1:1") %>% 
          ungroup() %>% 
          group_by(rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataMixture1 %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(rep, OTU) %>%
          filter(ratio %in% c("1:0","1:10","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/10 + relAbundanceOld[ratio=="0:1"]) %>% 
          filter(ratio=="1:10") %>% 
          ungroup() %>% 
          group_by(rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataMixture1 %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(rep, OTU) %>%
          filter(ratio %in% c("1:0","1:100","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/100 + relAbundanceOld[ratio=="0:1"]) %>% 
          filter(ratio=="1:100") %>% 
          ungroup() %>% 
          group_by(rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataMixture1 %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          group_by(rep, OTU) %>%
          filter(ratio %in% c("1:0","1:1000","0:1")) %>% 
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/1000 + relAbundanceOld[ratio=="0:1"]) %>% 
          filter(ratio=="1:1000") %>% 
          ungroup() %>% 
          group_by(rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
  ) %>% 
  mutate(relAbundance = ifelse(is.na(relAbundanceOld) | relAbundanceOld<10^-3, 10^-3, relAbundanceOld),
         logAbundance = log10(relAbundance))

# Add resequenced passage data to mix 1.
S7S17p3 <- fread("../../../data/S7S17p3.txt")
dataMixture1Theoretical <- dataMixture1Theoretical %>% 
  mutate(passage = 5) %>% 
  rbind(S7S17p3 %>% select(-c(combo, Family, count, plate, communityType, community, 
                              Kingdom, Phylum, Class, Order, Genus, relAbundanceOld, ))) %>% 
  ungroup() %>% 
  mutate(ratio = fct_relevel(ratio, ratiosFull))

# Export cleaned intermediate file containing just datasets shown in paper (but including contaminants).
dataPaper <- dataPlottingTheoretical %>% 
  filter(combo %in% c("7-17","XFA-17")) %>% 
  mutate(passage = 5) %>% 
  select(-c(taxa,hex,OTUnum,col)) %>% 
  rbind(S7S17p3 %>% select(!plate)) %>% 
  mutate(combo = case_when(combo=="7-17" ~ "Strep7-Strep17",
                           combo=="XFA-17" ~ "D1-Strep17",
                           combo=="Strep7-Strep17" ~ "Strep7-Strep17"))

write.table(dataPaper, "e0043MixtureDataframe.txt", quote=FALSE, row.names=FALSE, sep="\t")

# Calculate dose difference stat on mixture 1.
dataMix1StatsOld <- dataMixture1Theoretical %>% 
  filter(mixtureType=="actual" & ratio %in% ratios) %>%
  filter(!(rep==3 & passage==5)) %>% 
  # Add dose distance/difference statistics.
  group_by(passage, rep, OTU) %>% 
  arrange(passage, OTU, rep, ratio) %>% 
  mutate(doseDist=abs(logAbundance-lag(logAbundance)),
         doseDistTotal=mean(doseDist, na.rm=TRUE),
         doseDiff=logAbundance-lag(logAbundance),
         doseDiffTotal=mean(doseDiff, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(passage, OTU) %>% 
  summarize(doseDistance=mean(doseDistTotal), 
         doseDifference=mean(doseDiffTotal)) %>% 
  unique() %>% 
  mutate(OTUID = case_when(OTU==OTU7 ~ "Strep0007",
                           OTU==OTU17 ~ "Strep0017",
                           .default=NA)) %>% 
  filter(!is.na(OTUID))

# Calculate model dose dependence statistic.
dataMix1Stats <- dataMixture1Theoretical %>%
  filter(mixtureType=="actual" & ratio %in% ratios) %>%
  # Remove points that are obviously due to contamination.
  filter(!(rep==3 & passage==5)) %>% 
  group_by(passage, rep, OTU) %>% 
  arrange(passage, OTU, rep, ratio) %>% 
  mutate(doseDist=abs(logAbundance-lag(logAbundance)),
         doseDistTotal=mean(doseDist, na.rm=TRUE),
         doseDiff=logAbundance-lag(logAbundance),
         doseDiffTotal=mean(doseDiff, na.rm=TRUE), 
         doseDependence = ifelse(
           relAbundance[ratio=="1:1000"]/relAbundance[ratio=="1000:1"] > 1,
           relAbundance[ratio=="1:1000"]/relAbundance[ratio=="1000:1"],
           relAbundance[ratio=="1000:1"]/relAbundance[ratio=="1:1000"])) %>% 
  ungroup() %>% 
  group_by(passage, OTU) %>% 
  summarize(doseDistance=mean(doseDistTotal), 
            doseDifference=mean(doseDiffTotal),
            doseDependence = mean(doseDependence)) %>% 
  unique() %>% 
  mutate(OTUID = case_when(OTU==OTU7 ~ "Strep0007",
                           OTU==OTU17 ~ "Strep0017",
                           .default=NA)) %>% 
  filter(!is.na(OTUID))
write.table(dataMix1Stats, "dataMix1Stats.txt", quote=FALSE, row.names=FALSE, sep="\t")

# Plots --------------------------------------------------

# Trajectory plots of 7-17 mixture, for paper.
p_mixture1Paper <- dataMixture1Theoretical %>% 
  ungroup() %>% 
  filter(OTU %in% c(OTU7, OTU17)) %>% 
  filter(!(rep==3 & passage==5)) %>% 
  filter(!(passage==3 & ratio %in% c("1:0","0:1"))) %>% 
  mutate(OTUID = ifelse(OTU==OTU7, "L. garvieae", "L. lactis"),
         OTUID = fct_relevel(OTUID, c("L. garvieae","L. lactis"))) %>%
  mutate(lineGroup = case_when(ratio %in% c("1:0","0:1") ~ "parent",
                               ratio %in% ratios & rep==1 & mixtureType=="theoretical" ~ "theoretical-1",
                               ratio %in% ratios & rep==2 & mixtureType=="theoretical" ~ "theoretical-2",
                               ratio %in% ratios & rep==3 & mixtureType=="theoretical" ~ "theoretical-3",
                               ratio %in% ratios & rep==1 & mixtureType=="actual" ~ "actual-1",
                               ratio %in% ratios & rep==2 & mixtureType=="actual" ~ "actual-2",
                               ratio %in% ratios & rep==3 & mixtureType=="actual" ~ "actual-3"),
         parentGroup = case_when(ratio %in% ratios ~ "mixture",
                                 ratio=="1:0" & rep==1 ~ "leftParent1",
                                 ratio=="1:0" & rep==2 ~ "leftParent2",
                                 ratio=="0:1" & rep==1 ~ "rightParent1",
                                 ratio=="0:1" & rep==2 ~ "rightParent2"),
         lineGroupStrain = case_when(lineGroup == "parent" ~ "parent", 
                                     lineGroup %in% c("theoretical-1","theoretical-2","theoretical-3") ~ "theoretical",
                                     lineGroup %in% c("actual-1","actual-2","actual-3") ~ paste0(lineGroup, "_", OTUID)),
         lineShape = ifelse(passage==3, "dotted", "solid"), 
         pointShape = ifelse(passage==3, "open", "closed")) %>%
  #mutate(OTUID = ifelse(OTUID=="L. garvieae","L. garvieae\n","L. lactis\n")) %>% 
  # Remove extra copies of theoretical line.
  filter(mixtureType == "actual" | rep==1) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, 
                group=interaction(lineGroup, parentGroup, passage), linetype=lineShape), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroupStrain, shape=pointShape), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical"="Theoretical",
                              "actual-1_EnteC-2"="Actual 1 EnteC2",
                              "actual-2_EnteC-2"="Actual 2 EnteC2",
                              "actual-3_EnteC-2"="Actual 3 EnteC2",
                              "actual-1_EnteC-3"="Actual 1 EnteC3",
                              "actual-2_EnteC-3"="Actual 2 EnteC3",
                              "actual-3_EnteC-3"="Actual 3 EnteC3",
                              "actual-1_Strep-7"="Actual 1 Strep7",
                              "actual-2_Strep-7"="Actual 2 Strep7",
                              "actual-3_Strep-7"="Actual 3 Strep7",
                              "actual-1_Strep-17"="Actual 1 Strep17",
                              "actual-2_Strep-17"="Actual 2 Strep17",
                              "actual-3_Strep-17"="Actual 3 Strep17"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical"="#808080",
                              "actual-1_E. faecalis"="#D3BC4A", 
                              "actual-2_E. faecalis"="#C1A92F",
                              "actual-3_E. faecalis"="#A48F28",
                              "actual-1_E. casseliflavus"="#D38D4A",
                              "actual-2_E. casseliflavus"="#BF732E",
                              "actual-3_E. casseliflavus"="#A46428",
                              "actual-1_L. garvieae"="#7EAAB4",
                              "actual-2_L. garvieae"="#6398A4",
                              "actual-3_L. garvieae"="#53838D",
                              "actual-1_L. lactis"="#4044B5",
                              "actual-2_L. lactis"="#343795",
                              "actual-3_L. lactis"="#2A2D79"), guide="none") +
  scale_linetype_manual(values=c("dotted"=2, "solid"=1), guide="none") +
  scale_shape_manual(values=c("open"=5, "closed"=19), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        strip.text=element_text(face="bold.italic")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~OTUID)
p_mixture1Paper
save_plot("out/mixture1PaperUpdated.pdf", p_mixture1Paper, base_width=2.5, base_height=1.6)

# Trajectory plots of XFA-17 mixture with both EnteC species, for paper.
p_mixture2PaperEnteC <- dataPlottingTheoretical %>% 
  ungroup() %>% 
  filter(combo=="XFA-17") %>% 
  filter(OTU %in% c(OTU7,OTU17,EnteC2,EnteC3)) %>%
  filter(rep!=3) %>% 
  mutate(OTUID = case_when(OTU==EnteC2 ~ "E. faecalis",
                           OTU==EnteC3 ~ "E. casseliflavus",
                           OTU==OTU7 ~ "L. garvieae",
                           OTU==OTU17 ~ "L. lactis"),
         OTUID = fct_relevel(OTUID, c("E. faecalis","E. casseliflavus","L. garvieae","L. lactis"))) %>%
  mutate(lineGroup = case_when(ratio %in% c("1:0","0:1") ~ "parent",
                               ratio %in% ratios & rep==1 & mixtureType=="theoretical" ~ "theoretical-1",
                               ratio %in% ratios & rep==2 & mixtureType=="theoretical" ~ "theoretical-2", 
                               ratio %in% ratios & rep==1 & mixtureType=="actual" ~ "actual-1",
                               ratio %in% ratios & rep==2 & mixtureType=="actual" ~ "actual-2"),
         parentGroup = case_when(ratio %in% ratios ~ "mixture",
                                 ratio=="1:0" & rep==1 ~ "leftParent1",
                                 ratio=="1:0" & rep==2 ~ "leftParent2",
                                 ratio=="0:1" & rep==1 ~ "rightParent1",
                                 ratio=="0:1" & rep==2 ~ "rightParent2"),
         lineGroupStrain = case_when(lineGroup == "parent" ~ "parent", 
                                     lineGroup == "theoretical-1" ~ "theoretical-1",
                                     lineGroup == "theoretical-2" ~ "theoretical-2",
                                     lineGroup %in% c("actual-1","actual-2") ~ paste0(lineGroup, "_", OTUID))) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, group=interaction(lineGroupStrain, parentGroup)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroupStrain), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical"="Theoretical",
                              "actual-1_EnteC-2"="Actual 1 EnteC2",
                              "actual-2_EnteC-2"="Actual 2 EnteC2",
                              "actual-3_EnteC-2"="Actual 3 EnteC2",
                              "actual-1_EnteC-3"="Actual 1 EnteC3",
                              "actual-2_EnteC-3"="Actual 2 EnteC3",
                              "actual-3_EnteC-3"="Actual 3 EnteC3",
                              "actual-1_Strep-7"="Actual 1 Strep7",
                              "actual-2_Strep-7"="Actual 2 Strep7",
                              "actual-3_Strep-7"="Actual 3 Strep7",
                              "actual-1_Strep-17"="Actual 1 Strep17",
                              "actual-2_Strep-17"="Actual 2 Strep17",
                              "actual-3_Strep-17"="Actual 3 Strep17"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080",
                              "theoretical-2"="#808080",
                              "theoretical-3"="#808080",
                              "actual-1_E. faecalis"="#D3BC4A", 
                              "actual-2_E. faecalis"="#C1A92F",
                              "actual-3_E. faecalis"="#A48F28",
                              "actual-1_E. casseliflavus"="#D38D4A",
                              "actual-2_E. casseliflavus"="#BF732E",
                              "actual-3_E. casseliflavus"="#A46428",
                              "actual-1_L. garvieae"="#7EAAB4",
                              "actual-2_L. garvieae"="#6398A4",
                              "actual-3_L. garvieae"="#53838D",
                              "actual-1_L. lactis"="#4044B5",
                              "actual-2_L. lactis"="#343795",
                              "actual-3_L. lactis"="#2A2D79"), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        strip.text=element_text(face="bold.italic")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~OTUID, nrow=1)
p_mixture2PaperEnteC
save_plot("out/mixture2PaperEnteC.pdf", p_mixture2PaperEnteC, base_width=4.45, base_height=1.6)
save_plot("out/mixture2PaperEnteC.png", p_mixture2PaperEnteC, base_width=4.45, base_height=1.6)