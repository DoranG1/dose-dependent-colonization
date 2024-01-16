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
  filter(!(combo=="7-XFB" & rep==3) & !(combo=="XFA-17" & rep==3)) %>% 
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
  filter(!(combo=="XFA-17" & ratio=="1:1000")) %>%
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

# Add CASEU data to mix 1.
Strep7Strep17p3 <- fread("../../../../e0032-strain-isolation-flow-cytometry/analysis/231109-DG-e0043-HKX802-Elim824768/strep7strep17p3.txt")

Strep7Strep17p3Clean <- Strep7Strep17p3 %>% 
  filter(R2>0.9) %>% 
  select(!R2) %>% 
  pivot_longer(col=c(Strep7,Strep17), names_to="OTU", values_to="relAbundanceOld") %>% 
  mutate(mixtureType = "actual",
         relAbundance = ifelse(relAbundanceOld > 10^-3, relAbundanceOld, 10^-3),
         logAbundance = log10(relAbundance),
         OTU = ifelse(OTU=="Strep7",OTU7,OTU17))

dataMixture1Theoretical <- dataMixture1Theoretical %>% 
  mutate(passage = 5) %>% 
  rbind(Strep7Strep17p3Clean %>% mutate(passage =3)) %>%
  ungroup() %>% 
  mutate(ratio = fct_relevel(ratio, ratiosFull))

# Assess complementation for mixture 2. -----------------------------------

# Get set of OTUs in e0043 XFA community.
OTUsInXFA <- unique((dataPlotting %>% filter(combo=="XFA-17" & ratio=="1:0"))$OTU)

# Assess complementation of OTU 17 with XFA ASVs.
dataComplementation17 <- foreach(x = OTUsInXFA, .combine="rbind") %do% {
  df <- dataStats %>% 
    filter(combo=="XFA-17" & OTU==x) %>% 
    ungroup() %>% 
    select(Family, OTU, OTUnum, doseDifference) %>%
    unique() %>% 
    mutate(OTU17doseDifference = (dataStats %>% 
                                    filter(OTU==OTU17 & combo=="XFA-17"))$doseDifference[1],
           complementationValue = abs(doseDifference + OTU17doseDifference))
}

# Compare XFB e0043 to e0017. ---------------------------------------------

# Get set of OTUs in e0043 XFB community (which do not compete with OTU 7).
OTUsInXFB <- unique((dataXFBAnnotated %>% filter(inXFB))$OTU)

# Read in data from main e0017 experiment.
datae0017 <- fread("../../../data/mixtureDataframe.txt") %>% 
  filter(combo=="XFA-XFB" & ratio=="0:1" & category!="lowAbundance")

# Get set of OTUs in e0017 XFB community which compete with OTU 7
OTUsInXFBe0017 <- unique(datae0017$OTU)

# Get set of OTUs in e0017 but not in e0043, one of which may be competing with OTU 7.
OTUsNotIne0043 <- OTUsInXFBe0017[which(!(OTUsInXFBe0017 %in% OTUsInXFB))]
# Get taxonomy of OTUs.
OTUsNotIne0043Tax <- datae0017 %>%
  filter(passage==5 & OTU %in% OTUsNotIne0043) %>% 
  filter(!(is.na(Phylum) & is.na(Class) & is.na(Order) & is.na(Genus))) %>% 
  ungroup() %>% 
  select(combo, Phylum, Class, Order, Family, Genus, OTU, OTUnum, category) %>% 
  unique()

# Assess complementation of OTU 7 in e0017 with XFB ASVs missing from e0043.
dataComplementation7 <- foreach(x = OTUsNotIne0043, .combine="rbind") %do% {
    df <- datae0017 %>% 
      filter(passage==5 & OTU==x) %>% 
      ungroup() %>% 
      select(OTUnum, Family, a_doseDifference) %>%
      unique() %>% 
      mutate(OTU7doseDifference = (datae0017 %>% 
                                     filter(passage==5 & OTU==OTU7) %>% 
                                     ungroup() %>% 
                                     select(a_doseDifference) %>% 
                                     unique())$a_doseDifference,
             complementationValue = abs(a_doseDifference + OTU7doseDifference))
  }

# Combine e0017 and e0043 relative abundance data.
dataRelAbundance <- rbind(
  data %>% 
    filter(community=="XFB") %>% 
    select(Family, OTU, relAbundance, col) %>% 
    mutate(rep=ifelse(col%%3!=0, col%%3, 3)) %>% 
    select(!col) %>% 
    mutate(experiment="e0043"),
  fread("../../../data/mixtureDataframe.txt") %>% 
    filter(combo=="XFA-XFB" & ratio=="0:1" & passage==5) %>% 
    select(Family, OTU, rel_abundance_old, inoculationReplicate) %>% 
    rename(relAbundance=rel_abundance_old, rep=inoculationReplicate) %>% 
    mutate(experiment="e0017")
) %>% 
  left_join(data %>% select(Family, OTU, taxa)) %>% 
  unique()

# Plots --------------------------------------------------

# Trajectory plot for 7-17 mixture.
p_mixture1 <- dataPlotting %>% 
  filter(combo == "7-17") %>%
  filter(OTU %in% c(OTU7,OTU17)) %>%
  #filter(Genus=="Lactococcus") %>% 
  mutate(OTUID = ifelse(OTU==OTU7, "Strep-7", "Strep-17")) %>% 
  filter(rep!=3) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, group=interaction(OTU, rep), color=OTUID)) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
  scale_color_manual(name="ASV", values=c("#BBD686","#6D2E46")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  #ggtitle("7-17") +
  DEFAULTS.THEME_PRINT
  #facet_wrap(~col)
p_mixture1
save_plot("out/mixture1.png", p_mixture1, base_width=3.5, base_height=1.95)

# Trajectory plot for XFA-17 mixture.
p_mixture2 <- dataPlotting %>% 
  filter(combo == "XFA-17") %>%
  filter(OTU %in% c(OTU7,OTU17)) %>%
  #filter(Genus=="Lactococcus") %>% 
  mutate(OTUID = ifelse(OTU==OTU7, "Strep-7", "Strep-17")) %>% 
  mutate(group = ifelse(ratio %in% c("1:0","0:1"), "parent", "mixture")) %>% 
  filter(rep!=3) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, group=interaction(OTU, rep, group), color=OTUID), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=OTUID), alpha=0.8, size=0.75) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
  scale_color_manual(name="ASV", values=c("#BBD686","#6D2E46")) +
  #ggtitle("XFA-17")
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  DEFAULTS.THEME_PRINT
  #facet_wrap(~col)
p_mixture2
save_plot("out/mixture2.png", p_mixture2, base_width=2.55, base_height=1.5)

# Trajectory plot for 7-XFB mixture.
p_mixture3 <- dataPlotting %>% 
  filter(combo == "7-XFB") %>%
  filter(OTU %in% c(OTU7,OTU17)) %>%
  #filter(Genus=="Lactococcus") %>% 
  filter(logAbundance >= -3) %>% 
  filter(rep!=3) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, group=interaction(OTU, rep), color=OTU)) +
  ylim(-3,0) +
  scale_color_manual(name="OTU", values=c("dark red", "dark blue"), labels=c("7","17")) +
  ggtitle("7-XFB")
  #theme(legend.position="none")
  #facet_wrap(~col)
p_mixture3
save_plot("out/mixture3.png", p_mixture3, nrow=1, ncol=1)

# Trajectory plot for 7-17 mixture with contaminants removed.
p_mixture1Corrected <- dataMixture1 %>% 
  filter(OTU %in% c(OTU7,OTU17)) %>%
  mutate(OTUID = ifelse(OTU==OTU7, "Strep-7", "Strep-17")) %>%
  mutate(group = ifelse(ratio %in% c("1:0","0:1"), "parent", "mixture")) %>% 
  #filter(Genus=="Lactococcus") %>% 
  filter(rep!=3) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, group=interaction(OTU, rep, group), color=OTUID), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=OTUID), alpha=0.8, size=0.75) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
  scale_color_manual(name="ASV", values=c("#BBD686","#6D2E46")) +
  #ggtitle("7-17") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  DEFAULTS.THEME_PRINT
  #facet_wrap(~col)
p_mixture1Corrected
save_plot("out/mixture1Corrected.png", p_mixture1Corrected, base_width=2.55, base_height=1.52)

# Trajectory plots of 7-17 mixture, for paper.
p_mixture1Paper <- dataMixture1Theoretical %>% 
  ungroup() %>% 
  filter(OTU %in% c(OTU7, OTU17)) %>% 
  filter(!(rep==3 & passage==5)) %>% 
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
  # Remove extra copies of theoretical line.
  filter(mixtureType == "actual" | rep==1) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, 
                group=interaction(lineGroup, parentGroup, passage), linetype=lineShape), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroupStrain, shape=pointShape), alpha=0.8, size=0.5) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
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
  scale_linetype_manual(values=c("dotted"=3, "solid"=1), guide="none") +
  scale_shape_manual(values=c("open"=1, "closed"=16), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        strip.text=element_text(face="bold.italic")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~OTUID)
p_mixture1Paper
save_plot("out/mixture1Paper.pdf", p_mixture1Paper, base_width=2.5, base_height=1.6)

# Trajectory plots of XFA-17 mixture, for paper.
p_mixture2Paper <- dataPlottingTheoretical %>% 
  ungroup() %>% 
  filter(combo=="XFA-17") %>% 
  filter(OTU %in% c(OTU7,OTU17)) %>%
  filter(rep!=3) %>% 
  mutate(OTUID = ifelse(OTU==OTU7, "Streptococcaceae-7", "Streptococcaceae-17"),
         OTUID = fct_relevel(OTUID, c("Streptococcaceae-7","Streptococcaceae-17"))) %>%
  mutate(lineGroup = case_when(ratio %in% c("1:0","0:1") ~ "parent",
                               ratio %in% ratios & rep==1 & mixtureType=="theoretical" ~ "theoretical-1",
                               ratio %in% ratios & rep==2 & mixtureType=="theoretical" ~ "theoretical-2", 
                               ratio %in% ratios & rep==1 & mixtureType=="actual" ~ "actual-1",
                               ratio %in% ratios & rep==2 & mixtureType=="actual" ~ "actual-2"),
         parentGroup = case_when(ratio %in% ratios ~ "mixture",
                                 ratio=="1:0" & rep==1 ~ "leftParent1",
                                 ratio=="1:0" & rep==2 ~ "leftParent2",
                                 ratio=="0:1" & rep==1 ~ "rightParent1",
                                 ratio=="0:1" & rep==2 ~ "rightParent2")) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroup, group=interaction(lineGroup, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroup), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
  scale_color_manual(name="Mixture type",
                     labels=c("parent"="Parent",
                              "theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2"), 
                     values=c("parent"="#4DAF4A",
                              "theoretical-1"="#808080", 
                              "theoretical-2"="#707070", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~OTUID)
p_mixture2Paper
save_plot("out/mixture2Paper.pdf", p_mixture2Paper, base_width=2.5, base_height=1.6)

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
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, group=interaction(lineGroupStrain, parentGroup)), alpha=0.6, size=0.5) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroupStrain), alpha=0.8, size=0.75) +
  #geom_hline(yintercept=c(0, -0.5, -1, -1.5, -2, -2.5, -3), linetype="dotted", alpha=0.6) +
  xlab("Ratio") +
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

# Trajectory plot for 7-XFB mixture with all ASVs.
p_mixture3AllASVs <- dataXFBAnnotated %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, group=interaction(OTU, rep), color=OTU)) +
  ylim(-3,0) +
  #scale_color_manual(name="OTU", values=c("dark red", "dark blue"), labels=c("7","17")) +
  theme(legend.position="none") +
  ggtitle("7-XFB")
p_mixture3AllASVs
save_plot("out/mixture3AllASVs.png", p_mixture3AllASVs, nrow=1, ncol=1)

# Trajectory plot for 7-XFB mixture with all ASVs that change comparably.
p_mixture3AllASVsChanging <- dataXFBAnnotated %>% 
  filter(abs(doseDifference) > 0.05) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, group=interaction(OTU, rep), color=OTU)) +
  ylim(-3,0) +
  #scale_color_manual(name="OTU", values=c("dark red", "dark blue"), labels=c("7","17")) +
  theme(legend.position="none") +
  ggtitle("7-XFB")
p_mixture3AllASVsChanging
save_plot("out/mixture3AllASVsChanging.png", p_mixture3AllASVsChanging, nrow=1, ncol=1)

# Plot family level composition of XFA and XFB
p_communityCompXFAXFB <- data %>% 
  filter(community %in% c("XFA","XFB")) %>% 
  ggplot() +
  geom_bar(aes(x=col, y=relAbundance, fill=taxa), stat="identity", size=0.1, color="black") +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~community, scales="free")
p_communityCompXFAXFB
save_plot("out/communityCompXFAXFB.png", p_communityCompXFAXFB, nrow=1, ncol=1)

# Plot family level composition of XFB in e0017 and e0043.
p_XFBComparison <- dataRelAbundance %>% 
  ggplot() +
  geom_bar(aes(x=rep, y=relAbundance, fill=taxa), stat="identity", size=0.1, color="black") +
  scale_fill_manual(values=taxaPaletteList) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~experiment, scales="free")
p_XFBComparison
save_plot("out/XFBcomparison.png", p_XFBComparison, nrow=1, ncol=1)

# Plot all ASVs in mixture 2 on same plot.
p_mixture2allASVs <- dataPlotting %>% 
  filter(combo == "XFA-17") %>%
  filter(rep!=3) %>% 
  filter(ratio %in% ratios) %>%
  filter(ratio!="1:1000") %>% 
  mutate(OTUID = paste0(Family, "-", OTUnum)) %>% 
  group_by(OTUID) %>% 
  filter(mean(logAbundance, na.rm=TRUE)!=-3) %>% 
  ggplot() +
  geom_line(aes(x=ratio, y=logAbundance, group=rep, color=taxa)) +
  geom_point(aes(x=ratio, y=logAbundance, color=taxa)) +
  xlab("Ratio") +
  scale_y_continuous(name="Relative abundance", 
                     limits=c(-3,0), breaks=c(-3,-2,-1,0), 
                     labels=c("1e-03", "1e-02", "1e-01", "1e+00")) +
  scale_color_manual(values=taxaPaletteList) +
  #ggtitle("XFA-17")
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  DEFAULTS.THEME_PRINT +
  theme(legend.position="none") +
  facet_wrap(~OTUID)
p_mixture2allASVs
save_plot("out/mixture2allASVs.png", p_mixture2allASVs, base_width=9, base_height=9)
