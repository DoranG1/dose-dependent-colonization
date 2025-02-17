library(data.table)
library(tidyverse)
library(cowplot)
library(forcats)
library(scales)

theme_set(theme_cowplot())
source("../../config/palettes/plotDefaults.R")

# Import data.
data <- fread("../../workflow/out/e0065/DADA2_output/ps_all.txt.gz")

# Import coverage.
coverage <- fread("../../workflow/out/e0065/DADA2_output/read_summary.txt")

# Bind coverage file to dataframe to get averages.
data <- left_join(data,
                  coverage %>% 
                    mutate(fullSample = sub(".fastq.gz","",sample)) %>% 
                    select(!sample))

ratios <- c("1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000")
ratiosFull <- c("1:0", "1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000", "0:1")

Bacte0126 <- "TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGACTGGTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGTCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGATGCTCGAAAGTGTGGGTATCAAACAGG"
Bacte0003 <- "TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGACTGGTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGTCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGATGCTCGAAAGTGTGGGTATCAAACAGG"
Tanne0007 <- "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGTGGTTAATTAAGTCAGCGGTGAAAGTTTGTGGCTCAACCATAAAATTGCCGTTGAAACTGGTTGACTTGAGTATATTTGAGGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCTTACTAAACTATAACTGACACTGAAGCACGAAAGCGTGGGGATCAAACAGG"

# Reformat parent communities to plot with mixtures.
dataPlotting <- data %>% 
  select(-c(sample, filename, round2plate, round1index, well, group, fullSample, reads.in, reads.out)) %>% 
  mutate(logAbundance = log10(relAbundance))
dataPlotting <- dataPlotting %>% 
  mutate(combo = ifelse(communityType=="mixture", community, NA)) %>% 
  filter(!is.na(combo)) %>% 
  rbind(
    # Bacte0126/Bacte0003
    dataPlotting %>% filter(community=="Bacte0126") %>% 
      mutate(combo = "Bacte0126-Bacte0003", ratio="1:0"),
    dataPlotting %>% filter(community=="Bacte0003") %>% 
      mutate(combo = "Bacte0126-Bacte0003", ratio="0:1"),
    # Bacte0126/Tanne0007
    dataPlotting %>% filter(community=="Bacte0126") %>% 
      mutate(combo = "Bacte0126-Tanne0007", ratio="1:0"),
    dataPlotting %>% filter(community=="Tanne0007") %>% 
      mutate(combo = "Bacte0126-Tanne0007", ratio="0:1"),
    # Bacte0003/Tanne0007
    dataPlotting %>% filter(community=="Bacte0003") %>% 
      mutate(combo = "Bacte0003-Tanne0007", ratio="1:0"),
    dataPlotting %>% filter(community=="Tanne0007") %>% 
      mutate(combo = "Bacte0003-Tanne0007", ratio="0:1")
    ) %>% 
  mutate(relAbundanceOld = relAbundance,
         relAbundance = ifelse(relAbundance < 10^-3, 10^-3, relAbundance),
         logAbundance = ifelse(logAbundance < -3, -3, logAbundance)) %>%  
  group_by(passage, combo) %>% 
  complete(nesting(OTU, Family),
           ratio = ratiosFull, rep = c(1,2,3),
           fill = list(logAbundance = -3, relAbundance = 10^-3)) %>% 
  mutate(combo = case_when(
    "Bacte0126-Bacte0003" %in% combo ~ "Bacte0126-Bacte0003",
    "Bacte0126-Tanne0007" %in% combo ~ "Bacte0126-Tanne0007",
    "Bacte0003-Tanne0007" %in% combo ~ "Bacte0003-Tanne0007")) %>% 
  ungroup() %>% 
  mutate(ratio = fct_relevel(ratio, ratiosFull),
         mixtureType = "actual")
  
dataTheoretical <- dataPlotting %>% 
  filter(ratio %in% ratios) %>% 
  rbind(data.frame(passage = 1, combo = "Bacte0126-Tanne0007",
                   OTU = c(rep(Bacte0126, 9), rep(Tanne0007, 9)),
                   Family = NA, ratio = rep(ratiosFull, 2), rep = 1,
                   count = NA, plate = NA, communityType = NA, community = NA, 
                   relAbundance = c(c(1,0.99900100,0.99009901,0.90909091,0.5,0.09090909,0.00990099,0.00100000,0),
                                    rev(c(1,0.99900100,0.99009901,0.90909091,0.5,0.09090909,0.00990099,0.00100000,0))),
                   Kingdom = NA, Phylum = NA, Class = NA, Order = NA, Genus = NA, 
                   relAbundanceOld = NA) %>% 
          mutate(logAbundance = ifelse(relAbundance!=0, log10(relAbundance), -3), mixtureType = ifelse(ratio %in% ratios, "theoretical", "actual"))) %>% 
  rbind(data.frame(passage = 1, combo = "Bacte0003-Tanne0007",
                   OTU = c(rep(Bacte0003, 9), rep(Tanne0007, 9)),
                   Family = NA, ratio = rep(ratiosFull, 2), rep = 1,
                   count = NA, plate = NA, communityType = NA, community = NA, 
                   relAbundance = c(c(1,0.99900100,0.99009901,0.90909091,0.5,0.09090909,0.00990099,0.00100000,0),
                                    rev(c(1,0.99900100,0.99009901,0.90909091,0.5,0.09090909,0.00990099,0.00100000,0))),
                   Kingdom = NA, Phylum = NA, Class = NA, Order = NA, Genus = NA, 
                   relAbundanceOld = NA) %>% 
          mutate(logAbundance = ifelse(relAbundance!=0, log10(relAbundance), -3), mixtureType = ifelse(ratio %in% ratios, "theoretical", "actual")))

# Export intermediate annotated data file.
dataPaper <- dataTheoretical %>% 
  filter(community %in% c("Bacte0126","Tanne0007","Bacte0126-Tanne0007") | (mixtureType=="theoretical" & combo=="Bacte0126-Tanne0007"))
write.table(dataPaper, "e0065MixtureDataframe.txt", quote=FALSE, row.names=FALSE, sep="\t")

dataStats <- dataTheoretical %>%
  filter(mixtureType=="actual" & ratio %in% ratios) %>% 
  group_by(combo, passage, rep, OTU) %>% 
  arrange(combo, passage, OTU, rep, ratio) %>% 
  mutate(doseDist=abs(logAbundance-lag(logAbundance)),
         doseDistTotal=mean(doseDist, na.rm=TRUE),
         doseDiff=logAbundance-lag(logAbundance),
         doseDiffTotal=mean(doseDiff, na.rm=TRUE),
         doseDependence = ifelse(relAbundance[ratio=="1:1000"]/relAbundance[ratio=="1000:1"] > 1,
                                 relAbundance[ratio=="1:1000"]/relAbundance[ratio=="1000:1"],
                                 relAbundance[ratio=="1000:1"]/relAbundance[ratio=="1:1000"])) %>% 
  ungroup() %>% 
  group_by(combo, passage, OTU) %>% 
  summarize(doseDistance=mean(doseDistTotal), 
            doseDifference=mean(doseDiffTotal),
            doseDependence = mean(doseDependence)) %>% 
  unique() %>% 
  mutate(OTUID = case_when(OTU==Bacte0003 ~ "Bacte0003",
                           OTU==Bacte0126 ~ "Bacte0126",
                           OTU==Tanne0007 ~ "Tanne0007",
                           .default=NA)) %>% 
  filter(!is.na(OTUID))
write.table(dataStats, "dataStats.txt", quote=FALSE, row.names=FALSE, sep="\t")

# Individual trajectory plots. --------------------------------------------------

# Trajectory plot for Bacte0126-Tanne0007.
p_B126T <- dataTheoretical %>% 
  filter(combo=="Bacte0126-Tanne0007" & ((mixtureType=="actual" & ratio %in% ratios) | passage==1)) %>% 
  mutate(passage = ifelse(mixtureType=="theoretical" | ratio %in% c("1:0","0:1"), 5, passage)) %>% 
  filter(OTU %in% c(Bacte0126, Tanne0007)) %>% 
  mutate(OTUID = case_when(OTU==Bacte0126 ~ "B. fragilis",
                           OTU==Tanne0007 ~ "P. goldsteinii")) %>%
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
                                 ratio=="1:0" & rep==3 ~ "leftParent3",
                                 ratio=="0:1" & rep==1 ~ "rightParent1",
                                 ratio=="0:1" & rep==2 ~ "rightParent2",
                                 ratio=="0:1" & rep==3 ~ "rightParent3"),
         lineGroupStrain = case_when(lineGroup == "parent" ~ "parent", 
                                     lineGroup %in% c("theoretical-1","theoretical-2","theoretical-3") ~ "theoretical",
                                     lineGroup %in% c("actual-1","actual-2","actual-3") ~ paste0(lineGroup, "_", OTUID))) %>%
  #mutate(OTUID = ifelse(OTUID=="B. fragilis (D1)", "B. fragilis (D1)\n", "P. goldsteinii\n")) %>% 
  ggplot() + 
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, 
                group=interaction(parentGroup,lineGroupStrain, passage),
                linetype=factor(passage)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroupStrain, shape=factor(passage)), alpha=0.8, size=0.5) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical"="Theoretical",
                              "actual-1_B. fragilis"="Actual 1 B. fragilis",
                              "actual-2_B. fragilis"="Actual 2 B. fragilis",
                              "actual-3_B. fragilis"="Actual 3 B. fragilis",
                              "actual-1_P. goldsteinii"="Actual 1 P. goldsteinii",
                              "actual-2_P. goldsteinii"="Actual 2 P. goldsteinii",
                              "actual-3_P. goldsteinii"="Actual 3 P. goldsteinii"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical"="#808080",
                              "actual-1_B. fragilis"="#943874", 
                              "actual-2_B. fragilis"="#762D5D",
                              "actual-3_B. fragilis"="#592245",
                              "actual-1_P. goldsteinii"="#DA2F35",
                              "actual-2_P. goldsteinii"="#B62025",
                              "actual-3_P. goldsteinii"="#9C1C20"), guide="none") +
  scale_linetype_manual(values=c("1"="dotted", "3"="dashed", "5"="solid")) +
  scale_shape_manual(values=c("1"=1, "3"=5, "5"=19)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        strip.text=element_text(face="bold.italic")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~OTUID)
p_B126T
save_plot("out/B126T.pdf", p_B126T, base_width=2.5, base_height=1.6)
save_plot("out/B126T.png", p_B126T, base_width=2.5, base_height=1.7)
