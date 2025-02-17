library(data.table)
library(tidyverse)
library(cowplot)
library(forcats)
library(scales)

theme_set(theme_cowplot())
source("../../config/palettes/plotDefaults.R")

# Import data.
data <- fread("../../workflow/out/e0043/DADA2_output/ps_all.txt.gz")

# Import coverage.
coverage <- fread("../../workflow/out/e0043/DADA2_output/read_summary.txt")

# Bind coverage file to dataframe to get averages.
data <- left_join(data,
                  coverage %>% 
                    mutate(fullSample = sub(".fastq.gz","",sample)) %>% 
                    select(!sample))

ratios <- c("1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000")
ratiosFull <- c("1:0", "1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000", "0:1")

Strep7 <- "TACGTAGGTCCCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGTGGTTTCTTAAGTCTGATGTAAAAGGCAGTGGCTCAACCATTGTGTGCATTGGAAACTGGGAGACTTGAGTGCAGGAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGAGGCGAAAGCGGCTCTCTGGCCTGTAACTGACACTGAGGCTCGAAAGCGTGGGGAGCAAACAGG"
Strep17 <- "TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGTGGTTTATTAAGTCTGGTGTAAAAGGCAGTGGCTCAACCATTGTATGCATTGGAAACTGGTAGACTTGAGTGCAGGAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCCTGTAACTGACACTGAGGCTCGAAAGCGTGGGGAGCAAACAGG"
EnteC2 <- "TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGAGACTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGG"
EnteC3 <- "TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGAGACTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCGAACAGG"

# Reformat parent communities to plot with mixtures.
dataPlotting <- data %>% 
  select(-c(sample, filename, round2plate, round1index, well, group, fullSample, reads.in, reads.out)) %>% 
  mutate(logAbundance = log10(relAbundance))
dataPlotting <- dataPlotting %>% 
  mutate(combo = ifelse(communityType=="mixture", community, NA)) %>% 
  filter(!is.na(combo)) %>% 
  rbind(
    # Strep7/Strep17
    dataPlotting %>% filter(community=="Strep7") %>% 
      mutate(combo = "Strep7-Strep17", ratio="1:0"),
    dataPlotting %>% filter(community=="Strep17") %>% 
      mutate(combo = "Strep7-Strep17", ratio="0:1"),
    # D1/Strep17
    dataPlotting %>% filter(community=="D1") %>% 
      mutate(combo = "D1-Strep17", ratio="1:0"),
    dataPlotting %>% filter(community=="Strep17") %>% 
      mutate(combo = "D1-Strep17", ratio="0:1")) %>% 
  mutate(relAbundanceOld = relAbundance,
         relAbundance = ifelse(relAbundance < 10^-3, 10^-3, relAbundance),
         logAbundance = ifelse(logAbundance < -3, -3, logAbundance)) %>%  
  group_by(passage, combo) %>% 
  complete(nesting(OTU, Family),
           ratio = ratiosFull, rep = c(1,2,3),
           fill = list(logAbundance = -3, relAbundance = 10^-3)) %>% 
  mutate(combo = case_when(
    "Strep7-Strep17" %in% combo ~ "Strep7-Strep17",
    "D1-Strep17" %in% combo ~ "D1-Strep17")) %>% 
  ungroup() %>% 
  mutate(ratio = fct_relevel(ratio, ratiosFull),
         mixtureType = "actual") %>% 
  # Add parent community 0:1 for Strep17 in D1-Strep17 mixture.
  mutate(relAbundanceOld = ifelse(combo=="D1-Strep17" & OTU==Strep17 & ratio=="0:1", 1, relAbundanceOld))
  
dataTheoretical <- dataPlotting %>% 
  filter(combo=="D1-Strep17") %>% 
  rbind(dataPlotting %>% 
          filter(combo=="D1-Strep17") %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          group_by(rep, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/1000,
                 ratio = "1000:1") %>% 
          group_by(rep, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm = TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          filter(combo=="D1-Strep17") %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          group_by(rep, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/100,
                 ratio = "100:1") %>% 
          group_by(combo, rep, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          filter(combo=="D1-Strep17") %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          group_by(combo, rep, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/10,
                 ratio = "10:1") %>% 
          group_by(combo, rep, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          filter(combo=="D1-Strep17") %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>%
          group_by(combo, rep, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"],
                 ratio = "1:1") %>% 
          group_by(combo, rep, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep) %>%  
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          filter(combo=="D1-Strep17") %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>%
          group_by(combo, rep, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/10 + relAbundanceOld[ratio=="0:1"],
                 ratio = "1:10") %>% 
          group_by(combo, rep, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          filter(combo=="D1-Strep17") %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          group_by(combo, rep, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/100 + relAbundanceOld[ratio=="0:1"],
                 ratio = "1:100") %>% 
          group_by(combo, rep, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataPlotting %>% 
          filter(combo=="D1-Strep17") %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          group_by(combo, rep, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/1000 + relAbundanceOld[ratio=="0:1"],
                 ratio = "1:1000") %>% 
          group_by(combo, rep, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical")
  ) %>% 
  ungroup() %>% 
  mutate(relAbundance = ifelse(is.na(relAbundanceOld) | relAbundanceOld<10^-3, 10^-3, relAbundanceOld),
         logAbundance = log10(relAbundance),
         ratio = fct_relevel(ratio, ratiosFull))

# Export S7-S17 mixture data for replotting as part of figs. 4F-5D.
dataS7S17 <- dataPlotting %>% 
  filter(combo=="Strep7-Strep17")
write.table(dataS7S17, "S7S17p3.txt", row.names=FALSE, quote=FALSE, sep="\t")

# Export XFA-S17 intermediate annotated data file for public repo (s7-S17 data is included in e0043 repo.)
dataPaper <- dataPlotting %>% 
  filter(combo=="D1-Strep17") %>% 
  select(!plate) %>% 
  mutate(combo = "XFA-Strep17")
write.table(dataPaper, "XFA-Strep17-p3-mixtureDataframe.txt", quote=FALSE, row.names=FALSE, sep="\t")

# Individual trajectory plots. --------------------------------------------------

# Trajectory plot for D1-Strep17, passage 3.
p_D1S17 <- dataTheoretical %>% 
  filter(combo=="D1-Strep17") %>% 
  filter(OTU %in% c(EnteC2, EnteC3, Strep7, Strep17)) %>% 
  mutate(OTUID = case_when(OTU==EnteC2 ~ "E. faecalis",
                           OTU==EnteC3 ~ "E. casseliflavus",
                           OTU==Strep7 ~ "L. garvieae",
                           OTU==Strep17 ~ "L. lactis")) %>%
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
                                     lineGroup == "theoretical-1" ~ "theoretical-1",
                                     lineGroup == "theoretical-2" ~ "theoretical-2",
                                     lineGroup == "theoretical-3" ~ "theoretical-3",
                                     lineGroup %in% c("actual-1","actual-2","actual-3") ~ paste0(lineGroup, "_", OTUID)),
         OTUID = fct_relevel(OTUID, c("E. faecalis", "E. casseliflavus", "L. garvieae", "L. lactis"))) %>%
  ggplot() + 
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, 
                group=interaction(parentGroup,lineGroupStrain, passage)), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroupStrain), alpha=0.8, size=0.5) +
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
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        strip.text=element_text(face="bold.italic")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_grid(~OTUID)
p_D1S17
save_plot("out/D1S17.pdf", p_D1S17, base_width=4.45, base_height=1.6)
save_plot("out/D1S17.png", p_D1S17, base_width=4.45, base_height=1.6)
