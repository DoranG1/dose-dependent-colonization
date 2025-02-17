library(data.table)
library(tidyverse)
library(cowplot)
library(forcats)
library(R.utils)
library(foreach)
library(scales)

theme_set(theme_cowplot())
source("../../../config/palettes/plotDefaults.R")

# Import taxa color palette.
palette <- read.table("../../../config/palettes/KCHcolors-Silva-partial.txt",
                      header=TRUE, stringsAsFactors = FALSE) %>% 
  mutate(taxaShort=gsub("*.*\\.","",taxa))

# Import coverage.
coverageFile <- fread("../../../config/read_summary.txt")

Strep7 <- "TACGTAGGTCCCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGTGGTTTCTTAAGTCTGATGTAAAAGGCAGTGGCTCAACCATTGTGTGCATTGGAAACTGGGAGACTTGAGTGCAGGAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGAGGCGAAAGCGGCTCTCTGGCCTGTAACTGACACTGAGGCTCGAAAGCGTGGGGAGCAAACAGG"
Strep17 <- "TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGTGGTTTATTAAGTCTGGTGTAAAAGGCAGTGGCTCAACCATTGTATGCATTGGAAACTGGTAGACTTGAGTGCAGGAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCCTGTAACTGACACTGAGGCTCGAAAGCGTGGGGAGCAAACAGG"
EnteC2 <- "TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGAGACTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGG"
EnteC3 <- "TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGAGACTTGAGTGCAGAAGAGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCGAACAGG"


ratios <- c("1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000")
ratiosFull <- c("1:0", "1000:1", "100:1", "10:1", "1:1", "1:10", "1:100", "1:1000", "0:1")

# Read in data and add sample variables.
data <- fread("../../out/e0049-e0050/DADA2_output/ps_all.txt.gz") %>% 
  mutate(row=substr(well, 1, 1), col=as.numeric(substr(well, 2, length(well))))

# Create stripped-down color palette containing only families present in dataset.
palette <- palette %>% 
  filter(taxaShort %in% sort(unique(data %>% filter(relAbundance > 0.01) %>% pull(Family))))
# Make a named list.
paletteVector <- palette$hex
names(paletteVector) <- palette$taxaShort

dataAnnotated <- data %>% 
  mutate(communityType = case_when(col %in% c("4","8","9") ~ "blank",
                                   row == "A" & !(col %in% c("4","8","9")) ~ "parent",
                                   row != "A" & !(col %in% c("4","8","9")) ~ "mixture"), 
         community = case_when(# Blanks.
                               communityType == "blank" ~ "blank",
                               # Parents. 
                               plate %in%c("e0049-3","e0049-5","e0049-8") & row == "A" &
                                 col%in% c("1","2","3") ~ "EnteC2",
                               plate %in%c("e0049-3","e0049-5","e0049-8") & row == "A" &
                                 col%in% c("5","6","7") ~ "EnteC3",
                               plate %in%c("e0049-3","e0049-5","e0049-8") & row == "A" &
                                 col%in% c("10","11","12") ~ "Strep17",
                               plate %in%c("e0050-A-3","e0050-A-5") & row == "A" &
                                 col%in% c("1","2","3") ~ "EnteC2",
                               plate %in%c("e0050-A-3","e0050-A-5") & row == "A" &
                                 col%in% c("5","6","7") ~ "EnteC3",
                               plate %in%c("e0050-A-3","e0050-A-5") & row == "A" &
                                 col%in% c("10","11","12") ~ "XFB",
                               plate %in%c("e0050-B-3","e0050-B-5") & row == "A" &
                                 col%in% c("1","2","3") ~ "XFA",
                               plate %in%c("e0050-B-3","e0050-B-5") & row == "A" &
                                 col%in% c("5","6","7") ~ "Strep7",
                               plate %in%c("e0050-B-3","e0050-B-5") & row == "A" &
                                 col%in% c("10","11","12") ~ "Strep17",
                               plate == "e0050-8" & row == "A" &
                                 col %in% c("1","2","3") ~ "XFA",
                               plate == "e0050-8" & row == "A" &
                                 col %in% c("5","6","7") ~ "blank",
                               plate == "e0050-8" & row == "A" &
                                 col %in% c("10","11","12") ~ "XFB",
                               # Mixtures.
                               plate %in% c("e0049-3","e0049-5","e0049-8") & row!="A" &
                                 col %in% c("1","2","3") ~ "EnteC2-Strep17",
                               plate %in% c("e0049-3","e0049-5","e0049-8") & row!="A" &
                                col %in% c("5","6","7") ~ "EnteC3-Strep17",
                              plate %in% c("e0049-3","e0049-5","e0049-8") & row!="A" &
                                 col %in% c("10","11","12") ~ "EnteC-Strep17",
                               plate %in% c("e0050-A-3", "e0050-A-5") & row!="A" &
                                 col %in% c("1","2","3") ~ "EnteC2-XFB",
                               plate %in% c("e0050-A-3", "e0050-A-5") & row!="A" &
                                 col %in% c("5","6","7") ~ "EnteC3-XFB",
                               plate %in% c("e0050-A-3", "e0050-A-5") & row!="A" &
                                 col %in% c("10","11","12") ~ "EnteC-XFB",
                               plate %in% c("e0050-B-3", "e0050-B-5") & row!="A" &
                                 col %in% c("1","2","3") ~ "XFA-XFB",
                               plate %in% c("e0050-B-3", "e0050-B-5") & row!="A" &
                                 col %in% c("5","6","7") ~ "Strep7-XFB",
                               plate %in% c("e0050-B-3", "e0050-B-5") & row!="A" &
                                 col %in% c("10","11","12") ~ "EnteC3-Strep17",
                               plate == "e0050-8" & row!="A" &
                                 col %in% c("1","2","3") ~ "EnteC2-XFB",
                               plate == "e0050-8" & row!="A" &
                                 col %in% c("5","6","7") ~ "EnteC3-XFB",
                               plate == "e0050-8" & row!="A" &
                                 col %in% c("10","11","12") ~ "Strep7-XFB"),
         ratio=case_when(row=="B" ~ "1000:1",
                         row=="C" ~ "100:1",
                         row=="D" ~ "10:1",
                         row=="E" ~ "1:1",
                         row=="F" ~ "1:10",
                         row=="G" ~ "1:100",
                         row=="H" ~ "1:1000"),
         rep = case_when(col %in% c("1","5","10") ~ 1,
                         col %in% c("2","6","11") ~ 2,
                         col %in% c("3","7","12") ~ 3),
         mixtureRep = ifelse(plate %in% c("e0050-B-3", "e0050-B-5") & row!="A" & col %in% c("10","11","12"),
                             2, 1),
         logAbundance = log10(relAbundance))

# Check coverage of all samples.
coverage <- dataAnnotated %>% 
  group_by(passage, row, col, plate, community, ratio, rep, mixtureRep) %>% 
  summarize(sumCounts = sum(count)) %>% 
  unique()

# Bind coverage file to dataframe to get averages.
dataCoverage <- left_join(dataAnnotated,
                          coverageFile %>% 
                            mutate(fullSample = sub(".fastq.gz","",sample)) %>% 
                            select(!sample)) %>%
  filter(community!="blank") %>% 
  group_by(filename) %>% 
  slice(1)

# Counts reads in dataset used in paper.
sum(dataCoverage %>% filter(!(community %in% c("EnteC-XFB","XFA","XFA-XFB")) & mixtureRep==1) %>% pull(reads.out))

# Reformat parent communities to plot with mixtures.
dataPlotting <- dataAnnotated %>%
  left_join(coverage) %>% 
  select(-c(sample, filename, round2plate, round1index, well, group, fullSample, row, col))
dataPlotting <- dataPlotting %>% 
  mutate(combo = ifelse(communityType == "mixture", community, NA)) %>% 
  filter(!is.na(combo)) %>% 
  rbind(# e0049.3
        dataPlotting %>% filter(plate == "e0049-3" & community == "EnteC2") %>% 
          mutate(combo = "EnteC2-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-3" & community == "EnteC2") %>% 
          mutate(combo = "EnteC-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-3" & community == "EnteC3") %>% 
          mutate(combo = "EnteC3-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-3" & community == "EnteC3") %>% 
          mutate(combo = "EnteC-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-3" & community == "Strep17") %>% 
          mutate(combo = "EnteC2-Strep17", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0049-3" & community == "Strep17") %>% 
          mutate(combo = "EnteC3-Strep17", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0049-3" & community == "Strep17") %>% 
          mutate(combo = "EnteC-Strep17", ratio = "0:1"),
        # e0049.5
        dataPlotting %>% filter(plate == "e0049-5" & community == "EnteC2") %>% 
          mutate(combo = "EnteC2-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-5" & community == "EnteC2") %>% 
          mutate(combo = "EnteC-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-5" & community == "EnteC3") %>% 
          mutate(combo = "EnteC3-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-5" & community == "EnteC3") %>% 
          mutate(combo = "EnteC-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-5" & community == "Strep17") %>% 
          mutate(combo = "EnteC2-Strep17", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0049-5" & community == "Strep17") %>% 
          mutate(combo = "EnteC3-Strep17", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0049-5" & community == "Strep17") %>% 
          mutate(combo = "EnteC-Strep17", ratio = "0:1"),
        # e0049.8
        dataPlotting %>% filter(plate == "e0049-8" & community == "EnteC2") %>% 
          mutate(combo = "EnteC2-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-8" & community == "EnteC2") %>% 
          mutate(combo = "EnteC-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-8" & community == "EnteC3") %>% 
          mutate(combo = "EnteC3-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-8" & community == "EnteC3") %>% 
          mutate(combo = "EnteC-Strep17", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0049-8" & community == "Strep17") %>% 
          mutate(combo = "EnteC2-Strep17", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0049-8" & community == "Strep17") %>% 
          mutate(combo = "EnteC3-Strep17", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0049-8" & community == "Strep17") %>% 
          mutate(combo = "EnteC-Strep17", ratio = "0:1"),
        # e0050.A.3
        dataPlotting %>% filter(plate == "e0050-A-3" & community == "EnteC2") %>% 
          mutate(combo = "EnteC2-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-A-3" & community == "EnteC2") %>% 
          mutate(combo = "EnteC-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-A-3" & community == "EnteC3") %>% 
          mutate(combo = "EnteC3-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-A-3" & community == "EnteC3") %>% 
          mutate(combo = "EnteC-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-A-3" & community == "XFB") %>% 
          mutate(combo = "EnteC2-XFB", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0050-A-3" & community == "XFB") %>% 
          mutate(combo = "EnteC3-XFB", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0050-A-3" & community == "XFB") %>% 
          mutate(combo = "EnteC-XFB", ratio = "0:1"),
        # e0050.A.5
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "EnteC2") %>% 
          mutate(combo = "EnteC2-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "EnteC2") %>% 
          mutate(combo = "EnteC-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "EnteC3") %>% 
          mutate(combo = "EnteC3-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "EnteC3") %>% 
          mutate(combo = "EnteC-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "XFB") %>% 
          mutate(combo = "EnteC2-XFB", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "XFB") %>% 
          mutate(combo = "EnteC3-XFB", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "XFB") %>% 
          mutate(combo = "EnteC-XFB", ratio = "0:1"),
        # e0050.B.3
        dataPlotting %>% filter(plate == "e0050-B-3" & community == "XFA") %>% 
          mutate(combo = "XFA-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-B-3" & community == "Strep7") %>% 
          mutate(combo = "Strep7-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-B-3" & community == "Strep17") %>% 
          mutate(combo = "EnteC3-Strep17", ratio = "0:1", mixtureRep = 2),
        dataPlotting %>% filter(plate == "e0050-A-3" & community == "XFB") %>% 
          mutate(combo = "XFA-XFB", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0050-A-3" & community == "XFB") %>% 
          mutate(combo = "Strep7-XFB", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0050-A-3" & community == "EnteC3") %>% 
          mutate(combo = "EnteC3-Strep17", ratio = "1:0", mixtureRep = 2),
        # e0050.B.5
        dataPlotting %>% filter(plate == "e0050-B-5" & community == "XFA") %>% 
          mutate(combo = "XFA-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-B-5" & community == "Strep7") %>% 
          mutate(combo = "Strep7-XFB", ratio = "1:0"),
        dataPlotting %>% filter(plate == "e0050-B-5" & community == "Strep17") %>% 
          mutate(combo = "EnteC3-Strep17", ratio = "0:1", mixtureRep = 2),
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "XFB") %>% 
          mutate(combo = "XFA-XFB", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "XFB") %>% 
          mutate(combo = "Strep7-XFB", ratio = "0:1"),
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "EnteC3") %>% 
          mutate(combo = "EnteC3-Strep17", ratio = "1:0", mixtureRep = 2),
        # e0050.8
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "EnteC2") %>% 
          mutate(combo = "EnteC2-XFB", ratio = "1:0", passage = 8),
        dataPlotting %>% filter(plate == "e0050-A-5" & community == "EnteC3") %>% 
          mutate(combo = "EnteC3-XFB", ratio = "1:0", passage = 8),
        dataPlotting %>% filter(plate == "e0050-B-5" & community == "Strep7") %>% 
          mutate(combo = "Strep7-XFB", ratio = "1:0", passage = 8),
        dataPlotting %>% filter(plate == "e0050-8" & community == "XFB") %>% 
          mutate(combo = "EnteC2-XFB", ratio = "0:1", passage = 8),
        dataPlotting %>% filter(plate == "e0050-8" & community == "XFB") %>% 
          mutate(combo = "EnteC3-XFB", ratio = "0:1", passage = 8),
        dataPlotting %>% filter(plate == "e0050-8" & community == "XFB") %>% 
          mutate(combo = "Strep7-XFB", ratio = "0:1", passage = 8)) %>% 
  mutate(relAbundanceOld = relAbundance,
         relAbundance = ifelse(relAbundance < 10^-3, 10^-3, relAbundance),
         logAbundance = ifelse(logAbundance < -3, -3, logAbundance)) %>%  
  group_by(passage, combo, mixtureRep) %>% 
  complete(nesting(OTU, Family),
           ratio = ratiosFull, rep = c(1,2,3),
           fill = list(logAbundance = -3, relAbundance = 10^-3)) %>% 
  mutate(combo = case_when(
    "EnteC2-Strep17" %in% combo ~ "EnteC2-Strep17",
    "EnteC3-Strep17" %in% combo ~ "EnteC3-Strep17",
    "EnteC-Strep17" %in% combo ~ "EnteC-Strep17",
    "EnteC2-XFB" %in% combo ~ "EnteC2-XFB",
    "EnteC3-XFB" %in% combo ~ "EnteC3-XFB",
    "EnteC-XFB" %in% combo ~ "EnteC-XFB",
    "XFA-XFB" %in% combo ~ "XFA-XFB",
    "Strep7-XFB" %in% combo ~ "Strep7-XFB")) %>% 
  ungroup() %>% 
  mutate(ratio = fct_relevel(ratio, ratiosFull))
  
# Get list of low coverage samples.
lowCoverageSamples <- dataPlotting %>%
  filter(!((combo %in% c("EnteC2-Strep17","EnteC3-Strep17","EnteC-Strep17") & 
              ratio %in% ratios & sumCounts >= 100) |
             (combo %in% c("EnteC2-Strep17","EnteC3-Strep17","EnteC-Strep17") & 
                ratio %in% c("1:0","0:1")) |
           (combo %in% c("EnteC2-XFB","EnteC3-XFB","EnteC-XFB","Strep7-XFB") & 
              ratio %in% c(ratios,"0:1") & sumCounts >= 1000) |
           (combo %in% c("EnteC2-XFB","EnteC3-XFB","EnteC-XFB","Strep7-XFB") &
              ratio == "1:0") |
           (combo == "XFA-XFB" & sumCounts >= 1000))) %>% 
  group_by(passage, combo, ratio, mixtureRep, rep, sumCounts) %>% 
  slice(1) %>% 
  select(c(passage, combo, ratio, mixtureRep, rep, sumCounts))
  
# Bind OTU numbers from e0017.
datae0017 <- fread("../../../config/mixtureDataframe.txt")
dataPlotting <- dataPlotting %>% 
  left_join(datae0017 %>% 
              ungroup() %>% 
              select(Family, OTU, OTUnum) %>% 
              unique(),
            by=c("Family","OTU"))

# Bind family coding.
familyCoding <- fread("../../../config/familyCodingFull.txt")
dataPlotting <- dataPlotting %>% 
  left_join(familyCoding %>% 
              rename(Family=family))

# Add dose difference stat.
dataStats <- dataPlotting %>% 
  filter(!(combo=="Strep7-XFB" & passage==5 & rep==1)) %>% 
  #filter(!(combo=="EnteC-XFB" & passage %in% c(3,5) & ratio=="1:1000")) %>% 
  filter(ratio %in% ratios) %>% 
  group_by(combo, rep, mixtureRep, passage, OTU) %>% 
  arrange(combo, rep, mixtureRep, passage, OTU, ratio) %>% 
  mutate(doseDist=abs(logAbundance-lag(logAbundance)),
         doseDistTotal=mean(doseDist, na.rm=TRUE),
         doseDiff=logAbundance-lag(logAbundance),
         doseDiffTotal=mean(doseDiff, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(combo, mixtureRep, passage, OTU) %>% 
  mutate(doseDistance=mean(doseDistTotal), 
         doseDifference=mean(doseDiffTotal))

# Generate theoretical mixtures.
dataTheoretical <- dataPlotting %>% 
  # Remove low coverage samples.
  left_join(lowCoverageSamples %>% 
              ungroup() %>% 
              select(-sumCounts) %>% 
              mutate(filterLowCoverage = TRUE),
            by=c("passage","combo","ratio","mixtureRep","rep")) %>% 
  filter(is.na(filterLowCoverage)) %>% 
  ungroup() %>% 
  group_by(passage, combo, ratio, rep, mixtureRep) %>% 
  filter(!is.nan(mean(sumCounts, na.rm=TRUE))) %>% 
  # Bind to theoretical.
  mutate(mixtureType = "actual")
  
# Alternative theoretical dataframe, do not remove low coverage samples (for public data export.)
dataTheoretical <- dataPlotting %>% 
  mutate(mixtureType = "actual")

dataTheoretical <- dataTheoretical %>% 
  rbind(dataTheoretical %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          group_by(combo, rep, mixtureRep, passage) %>% 
          filter("1:0" %in% ratio & "0:1" %in% ratio) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/1000,
                 ratio = "1000:1") %>% 
          group_by(combo, rep, mixtureRep, passage, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm = TRUE),
                 mixtureType = "theoretical"),
        dataTheoretical %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          group_by(combo, rep, mixtureRep, passage) %>% 
          filter("1:0" %in% ratio & "0:1" %in% ratio) %>% 
          ungroup() %>%
          group_by(combo, rep, mixtureRep, passage, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/100,
                 ratio = "100:1") %>% 
          group_by(combo, rep, mixtureRep, passage, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataTheoretical %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          group_by(combo, rep, mixtureRep, passage) %>%  
          filter("1:0" %in% ratio & "0:1" %in% ratio) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"]/10,
                 ratio = "10:1") %>% 
          group_by(combo, rep, mixtureRep, passage, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataTheoretical %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>%
          group_by(combo, rep, mixtureRep, passage) %>% 
          filter("1:0" %in% ratio & "0:1" %in% ratio) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"] + relAbundanceOld[ratio=="0:1"],
                 ratio = "1:1") %>% 
          group_by(combo, rep, mixtureRep, passage, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage) %>%  
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataTheoretical %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>%
          group_by(combo, rep, mixtureRep, passage) %>%  
          filter("1:0" %in% ratio & "0:1" %in% ratio) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/10 + relAbundanceOld[ratio=="0:1"],
                 ratio = "1:10") %>% 
          group_by(combo, rep, mixtureRep, passage, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataTheoretical %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          group_by(combo, rep, mixtureRep, passage) %>% 
          filter("1:0" %in% ratio & "0:1" %in% ratio) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/100 + relAbundanceOld[ratio=="0:1"],
                 ratio = "1:100") %>% 
          group_by(combo, rep, mixtureRep, passage, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical"),
        dataTheoretical %>% 
          mutate(relAbundanceOld = ifelse(is.na(relAbundanceOld), 0, relAbundanceOld)) %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          group_by(combo, rep, mixtureRep, passage) %>% 
          filter("1:0" %in% ratio & "0:1" %in% ratio) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage, OTU) %>%
          mutate(relAbundanceOld = relAbundanceOld[ratio=="1:0"]/1000 + relAbundanceOld[ratio=="0:1"],
                 ratio = "1:1000") %>% 
          group_by(combo, rep, mixtureRep, passage, OTU, ratio) %>% 
          slice(1) %>% 
          ungroup() %>% 
          group_by(combo, rep, mixtureRep, passage) %>% 
          mutate(relAbundanceOld = relAbundanceOld / sum(relAbundanceOld, na.rm=TRUE),
                 mixtureType = "theoretical")
  ) %>% 
  ungroup() %>% 
  mutate(relAbundance = ifelse(is.na(relAbundanceOld) | relAbundanceOld<10^-3, 10^-3, relAbundanceOld),
         logAbundance = log10(relAbundance),
         ratio = fct_relevel(ratio, ratiosFull))

# Export cleaned intermediate file containing just datasets shown in paper (but including contaminants).
dataPaper <- dataTheoretical %>% 
  filter(combo %in% c("EnteC-Strep17","EnteC2-Strep17","EnteC3-Strep17","EnteC2-XFB","EnteC3-XFB","Strep7-XFB")) %>% 
  filter(mixtureRep == 1) %>% 
  select(-c(plate,mixtureRep,OTUnum,code,sumCounts))
write.table(dataPaper, "e0049e0050MixtureDataframe.txt", quote=FALSE, row.names=FALSE, sep="\t")

# Calculate new DD statistics.
dataStatsModel <- dataTheoretical %>% 
  filter(mixtureRep==1) %>% 
  filter(mixtureType=="actual" & ratio %in% ratios) %>%
  filter(combo %in% c("EnteC2-Strep17","EnteC3-Strep17") & passage %in% c(3,5)) %>% 
  group_by(combo, passage, rep, OTU) %>% 
  arrange(combo, passage, OTU, rep, ratio) %>% 
  mutate(doseDist=abs(logAbundance-lag(logAbundance)),
         doseDistTotal=mean(doseDist, na.rm=TRUE),
         doseDiff=logAbundance-lag(logAbundance),
         doseDiffTotal=mean(doseDiff, na.rm=TRUE),
         doseDependence = ifelse(
           relAbundance[ratio=="1:1000"]/relAbundance[ratio=="1000:1"] > 1,
           relAbundance[ratio=="1:1000"]/relAbundance[ratio=="1000:1"],
           relAbundance[ratio=="1000:1"]/relAbundance[ratio=="1:1000"])) %>% 
  ungroup() %>% 
  group_by(combo, passage, OTU) %>% 
  summarize(doseDistance=mean(doseDistTotal), 
            doseDifference=mean(doseDiffTotal),
            doseDependence = mean(doseDependence)) %>% 
  unique() %>% 
  mutate(OTUID = case_when(OTU==EnteC2 ~ "EnteC0002",
                           OTU==EnteC3 ~ "EnteC0003",
                           OTU==Strep17 ~ "Strep0017",
                           .default=NA)) %>% 
  filter(!is.na(OTUID)) %>% 
  filter(combo=="EnteC2-Strep17" & OTUID %in% c("EnteC0002","Strep0017") |
           combo=="EnteC3-Strep17" & OTUID %in% c("EnteC0003","Strep0017"))
write.table(dataStatsModel, "dataStatsModel.txt", quote=FALSE, row.names=FALSE, sep="\t")

# Individual trajectory plots. --------------------------------------------------

# Plot EnteC2-Strep17 with both passages on same plot.
p_EnteC2Strep17BothPassages <- dataTheoretical %>% 
  #filter(passage==5 | mixtureType=="theoretical") %>% 
  filter(combo=="EnteC2-Strep17" & mixtureRep==1 & passage %in% c(3,5)) %>% 
  filter(OTU %in% c(EnteC2, Strep17)) %>% 
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
                                       lineGroup %in% c("theoretical-1","theoretical-2","theoretical-3") ~ "theoretical",
                                       lineGroup %in% c("actual-1","actual-2","actual-3") ~ paste0(lineGroup, "_", OTUID)),
         lineShape = ifelse(passage==3, "dotted", "solid"),
         pointShape = ifelse(passage==3, "open", "closed")) %>% 
  #mutate(OTUID = ifelse(OTUID=="E. faecalis", "E. faecalis\n", "L. lactis\n")) %>% 
  # Remove extra copies of theoretical line.
  filter((mixtureType == "actual" & ratio %in% ratios) | (rep==1 & passage==5)) %>% 
  ggplot() + 
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, 
                group=interaction(lineGroupStrain, parentGroup, passage), linetype=lineShape), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroupStrain, shape=lineShape), alpha=0.8, size=0.5) +
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
  scale_linetype_manual(name = "Passage", 
                        labels = c("dotted" = "Passage 3",
                                   "solid" = "Passage 5"),
                        values=c("dotted"=2, "solid"=1)) +
  scale_shape_manual(name = "Passage", 
                     labels = c("dotted" = "Passage 3",
                                "solid" = "Passage 5"),
                     values = c("dotted"=5, "solid"=19)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        strip.text=element_text(face="bold.italic")) +
  guides(linetype=guide_legend(ncol=2)) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~OTUID)
p_EnteC2Strep17BothPassages
save_plot("out/EnteC2Strep17BothPassages.pdf", p_EnteC2Strep17BothPassages, base_width=2.5, base_height=1.6)
save_plot("out/EnteC2Strep17BothPassages.png", p_EnteC2Strep17BothPassages, base_width=2.5, base_height=1.6)

# Plot EnteC3-Strep17 with both passages on same plot.
p_EnteC3Strep17BothPassages <- dataTheoretical %>% 
  filter(combo=="EnteC3-Strep17" & mixtureRep==1 & passage %in% c(3,5)) %>% 
  filter(OTU %in% c(EnteC3, Strep17)) %>% 
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
                                     lineGroup %in% c("theoretical-1","theoretical-2","theoretical-3") ~ "theoretical",
                                     lineGroup %in% c("actual-1","actual-2","actual-3") ~ paste0(lineGroup, "_", OTUID)),
         lineShape = ifelse(passage==3, "dotted", "solid"), 
         pointShape = ifelse(passage==3, "open", "closed")) %>%
  #mutate(OTUID = ifelse(OTUID=="E. casseliflavus","E. casseliflavus\n", "L. lactis\n")) %>% 
  # Remove extra copies of theoretical line.
  filter((mixtureType == "actual" & ratio %in% ratios) | (rep==1 & passage==5)) %>%
  ggplot() + 
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, 
                group=interaction(lineGroupStrain, parentGroup, passage), linetype=lineShape), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroupStrain, shape=pointShape), alpha=0.8, size=0.5) +
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
p_EnteC3Strep17BothPassages
save_plot("out/EnteC3Strep17BothPassages.pdf", p_EnteC3Strep17BothPassages, base_width=2.5, base_height=1.6)
save_plot("out/EnteC3Strep17BothPassages.png", p_EnteC3Strep17BothPassages, base_width=2.5, base_height=1.6)

# Plot EnteC2-EnteC3-Strep17 with both passages on same plot.
p_EnteCStrep17BothPassages <- dataTheoretical %>% 
  #filter(passage==3 | mixtureType=="theoretical") %>% 
  filter(combo=="EnteC-Strep17" & mixtureRep==1 & passage %in% c(3,5)) %>%
  filter(OTU %in% c(EnteC2, EnteC3, Strep17)) %>%
  mutate(OTUID = case_when(OTU==EnteC2 ~ "E. faecalis",
                           OTU==EnteC3 ~ "E. casseliflavus",
                           OTU==Strep7 ~ "L. garvieae",
                           OTU==Strep17 ~ "L. lactis"),
         OTUID = fct_relevel(OTUID, c("E. faecalis","E. casseliflavus","L. garvieae","L. lactis"))) %>%
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
                                     lineGroup %in% c("actual-1","actual-2","actual-3") ~ paste0(lineGroup, "_", OTUID)),
         lineShape = ifelse(passage==3, "dotted", "solid"), 
         pointShape = ifelse(passage==3, "open", "closed")) %>%
  # Remove extra copies of theoretical line.
  filter(mixtureType == "actual" | (rep==1 & passage==5)) %>% 
  ggplot() + 
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, 
                group=interaction(lineGroupStrain, parentGroup, passage), linetype=lineShape), alpha=0.6, size=0.4) +
  geom_point(aes(x=ratio, y=logAbundance, color=lineGroupStrain, shape=pointShape), alpha=0.8, size=0.5) +
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
  scale_linetype_manual(values=c("dotted"=3, "solid"=1), guide="none") +
  scale_shape_manual(values=c("open"=1, "closed"=19), guide="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines"),
        strip.text=element_text(face="bold.italic")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~OTUID)
p_EnteCStrep17BothPassages
save_plot("out/EnteCStrep17BothPassages.pdf", p_EnteCStrep17BothPassages, base_width=3.5, base_height=1.6)
save_plot("out/EnteCStrep17BothPassages.png", p_EnteCStrep17BothPassages, base_width=2.5, base_height=1.6)

# Trajectory plots for Strep7-XFB.
p_Strep7XFBp8 <- dataTheoretical %>% 
  #filter(combo=="Strep7-XFB" & passage==3) %>% 
  filter(combo=="Strep7-XFB" & passage==8) %>% 
  filter(OTU %in% c(EnteC2, EnteC3, Strep7, Strep17)) %>%
  mutate(OTUID = case_when(OTU==EnteC2 ~ "E. faecalis",
                           OTU==EnteC3 ~ "E. casseliflavus",
                           OTU==Strep7 ~ "L. garvieae",
                           OTU==Strep17 ~ "L. lactis"),
         OTUID = fct_relevel(OTUID, c("E. faecalis","E. casseliflavus","L. garvieae","L. lactis"))) %>%
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
                                     lineGroup %in% c("actual-1","actual-2","actual-3") ~ paste0(lineGroup, "_", OTUID))) %>%
  # Filter noisy replicate in p5.
  #filter(rep!=1) %>%
  # Filter out noise in p8.
  filter(ratio!="1000:1") %>% 
  ggplot() + 
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, group=interaction(lineGroupStrain, parentGroup)), alpha=0.6, size=0.4) +
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
p_Strep7XFBp8
save_plot("out/Strep7XFBp8.pdf", p_Strep7XFBp8, base_width=4.45, base_height=1.6)
save_plot("out/Strep7XFBp8.png", p_Strep7XFBp8, base_width=4.45, base_height=1.6)

# Trajectory plots for EnteC2-XFB.
p_EnteC2XFBp8 <- dataTheoretical %>% 
  #filter(combo=="EnteC2-XFB" & passage==3) %>% 
  #filter(combo=="EnteC2-XFB" & passage==5) %>%
  filter(combo=="EnteC2-XFB" & passage==8 & rep!=1) %>%
  filter(OTU %in% c(EnteC2, EnteC3, Strep7, Strep17)) %>%
  mutate(OTUID = case_when(OTU==EnteC2 ~ "E. faecalis",
                           OTU==EnteC3 ~ "E. casseliflavus",
                           OTU==Strep7 ~ "L. garvieae",
                           OTU==Strep17 ~ "L. lactis"),
         OTUID = fct_relevel(OTUID, c("E. faecalis","E. casseliflavus","L. garvieae","L. lactis"))) %>%
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
                                     lineGroup %in% c("actual-1","actual-2","actual-3") ~ paste0(lineGroup, "_", OTUID))) %>%
  ggplot() + 
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, group=interaction(lineGroupStrain, parentGroup)), alpha=0.6, size=0.4) +
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
                              "theoretical-1"="#808080",
                              "theoretical-2"="#808080",
                              "theoretical-2"="#808080",
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
  facet_wrap(~OTUID)
p_EnteC2XFBp8
save_plot("out/EnteC2XFBp8.pdf", p_EnteC2XFBp8, base_width=3.5, base_height=1.6)
save_plot("out/EnteC2XFBp8.png", p_EnteC2XFBp8, base_width=3.5, base_height=1.6)

# Trajectory plots for EnteC3-XFB.
p_EnteC3XFBp8 <- dataTheoretical %>% 
  #filter(combo=="EnteC3-XFB" & passage==5) %>%
  filter(combo=="EnteC3-XFB" & passage==8) %>% 
  filter(OTU %in% c(EnteC2, EnteC3, Strep7, Strep17)) %>%
  mutate(OTUID = case_when(OTU==EnteC2 ~ "E. faecalis",
                           OTU==EnteC3 ~ "E. casseliflavus",
                           OTU==Strep7 ~ "L. garvieae",
                           OTU==Strep17 ~ "L. lactis"),
         OTUID = fct_relevel(OTUID, c("E. faecalis","E. casseliflavus","L. garvieae","L. lactis"))) %>%
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
                                     lineGroup %in% c("actual-1","actual-2","actual-3") ~ paste0(lineGroup, "_", OTUID))) %>%
  # Remove noise in p5.
  #filter(rep!=2) %>%
  # Remove noise in p8.
  filter(rep!=1) %>% 
  ggplot() + 
  geom_line(aes(x=ratio, y=logAbundance, color=lineGroupStrain, group=interaction(lineGroupStrain, parentGroup)), alpha=0.6, size=0.4) +
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
                              "theoretical-1"="#808080",
                              "theoretical-2"="#808080",
                              "theoretical-2"="#808080",
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
  facet_wrap(~OTUID)
p_EnteC3XFBp8
save_plot("out/EnteC3XFBp8.pdf", p_EnteC3XFBp8, base_width=3.5, base_height=1.6)
save_plot("out/EnteC3XFBp8.png", p_EnteC3XFBp8, base_width=3.5, base_height=1.6)