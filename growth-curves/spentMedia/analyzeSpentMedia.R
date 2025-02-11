library(tidyverse)
library(cowplot)
library(foreach)
library(RColorBrewer)
library(gcplyr)
library(data.table)
library(scales)
library(viridis)

theme_set(theme_cowplot())
source("../../../config/palettes/plotDefaults.R")

# Import spent media data from replicate 1.
rep1Data <- fread("../../../data/e0069/e0069.A1-cleaned.txt") %>% 
  select(!`T° 600`) %>% 
  separate(Time, into=c("hour", "minute", "second"), sep=":") %>% 
  mutate(time = as.numeric(as.integer(hour)*60 + as.integer(minute) + as.integer(second)/60)) %>% 
  select(-c(hour, minute, second)) %>% 
  pivot_longer(!time, names_to="well", values_to="OD600") %>% 
  mutate(row = substr(well, 1, 1), col = as.numeric(substr(well, 2, length(well))))
# Get minimum blank value from plate.
rep1Min <- min(rep1Data$OD600)

# Import spent media data from replicate 2.
rep2Data <- fread("../../../data/e0069/e0069.B1-cleaned.txt") %>% 
  select(!`T° 600`) %>% 
  separate(Time, into=c("hour", "minute", "second"), sep=":") %>% 
  mutate(time = as.numeric(as.integer(hour)*60 + as.integer(minute) + as.integer(second)/60)) %>% 
  select(-c(hour, minute, second)) %>% 
  pivot_longer(!time, names_to="well", values_to="OD600") %>% 
  mutate(row = substr(well, 1, 1), col = as.numeric(substr(well, 2, length(well))))
# Get minimum blank value from plate.
rep2Min <- min(rep2Data$OD600)

# Combine replicates and subtract blanks.
dataCombined <- rbind(
  rep1Data %>% mutate(rep = 1),
  rep2Data %>% mutate(rep = 2)
  ) %>% 
  # left_join(rbind(rep1BlankLastTimepoints %>% mutate(rep = 1) %>% rename(blankOD = finalOD),
  #                 rep2BlankLastTimepoints %>% mutate(rep = 2) %>% rename(blankOD = finalOD))) %>% 
  # mutate(blankSubtractedOD = OD600-blankOD,
  #        blankSubtractedOD = ifelse(blankSubtractedOD<0.005, 0.005, blankSubtractedOD))
  # mutate(blankSubtractedOD = ifelse(rep==1, OD600 - rep1MedBlank, OD600 - rep2MedBlank),
  #        blankSubtractedOD = ifelse(blankSubtractedOD<0.005, 0.005, blankSubtractedOD))
  # Blank-subtract by subtracting average of first three timepoints from each well.
  group_by(rep, well) %>% 
  mutate(blank = mean(OD600[time %in% unique(time)[1:3]]),
         blankSubtractedOD = OD600 - blank,
         blankSubtractedOD = ifelse(blankSubtractedOD<0.005, 0.005, blankSubtractedOD))

# Annotate strains and spent media in combined dataframe.
dataCombined <- dataCombined %>% 
  mutate(strain = case_when(row %in% c("A","H") | col %in% c(1,11,12) ~ "blank",
                            row=="B" & col %in% c(2,3,4,5,6,7,8,9,10) ~ "B. fragilis",
                            row=="C" & col %in% c(2,3,4,5,6,7,8,9,10) ~ "P. goldsteinii",
                            row=="D" & col %in% c(2,3,4,5,6,7,8,9,10) ~ "E. faecalis",
                            row=="E" & col %in% c(2,3,4,5,6,7,8,9,10) ~ "E. casseliflavus",
                            row=="F" & col %in% c(2,3,4,5,6,7,8,9,10) ~ "L. garvieae",
                            row=="G" & col %in% c(2,3,4,5,6,7,8,9,10) ~ "L. lactis"),
         spentMedia = case_when((row %in% c("A","H") & col %in% c(1,2,3,4,5,6,7,8,9,10,12)) |
                                  col %in% c(1,10,12) ~ "fresh",
                                (col==2 & row %in% c("B","C","D","E","F","G")) | 
                                  well=="B11" ~ "B. fragilis",
                                (col==3 & row %in% c("B","C","D","E","F","G")) | 
                                  well=="C11" ~ "P. goldsteinii",
                                (col==4 & row %in% c("B","C","D","E","F","G")) | 
                                  well=="D11" ~ "E. faecalis",
                                (col==5 & row %in% c("B","C","D","E","F","G")) | 
                                  well=="E11" ~ "E. casseliflavus",
                                (col==6 & row %in% c("B","C","D","E","F","G")) | 
                                  well=="F11" ~ "L. garvieae",
                                (col==7 & row %in% c("B","C","D","E","F","G")) | 
                                  well=="G11" ~ "L. lactis",
                                (col==8 & row %in% c("B","C","D","E","F","G")) | 
                                  well=="A11" ~ "Community D1",
                                (col==9 & row %in% c("B","C","D","E","F","G")) | 
                                  well=="H11" ~ "Community D2"))

# Calculate fold-changes separately for each replicate of the strain and the fresh media and then average.
nicheOverlapsSeparateReps <- dataCombined %>% 
  filter(strain!="blank" & spentMedia!="fresh") %>% 
  filter(time<=1440) %>% 
  group_by(strain, spentMedia, rep) %>% 
  summarize(maxOD = max(blankSubtractedOD)) %>% 
  left_join(dataCombined %>% 
              filter(strain!="blank" & spentMedia=="fresh") %>% 
              group_by(strain, rep) %>% 
              summarize(freshOD = max(blankSubtractedOD))) %>% 
  mutate(foldChange = maxOD/freshOD,
         nicheOverlap = 1-foldChange) %>% 
  ungroup() %>% 
  group_by(strain, spentMedia) %>% 
  summarize(foldChange = mean(foldChange),
            nicheOverlap = mean(nicheOverlap))

# Plots. ------------------------------------------------------------------

# Plot a schematic curve of E. faecalis in the L. lactis spent media.
faecalisLactisExample <- dataCombined %>% 
  ungroup %>% 
  filter(time <= 1440) %>% 
  filter(strain=="E. faecalis" & spentMedia %in% c("L. lactis", "fresh")) %>%
  mutate(spentMedia = ifelse(spentMedia=="fresh", "Fresh mBHI", "Strain 2")) %>% 
  mutate(strain = "Strain  1") %>% 
  group_by(spentMedia, strain, time) %>% 
  summarize(avgOD = mean(blankSubtractedOD),
            maxOD = max(blankSubtractedOD),
            minOD = min(blankSubtractedOD)) %>% 
  # Smooth curves over 1/2 hr intervals.
  group_by(spentMedia) %>% 
  arrange(spentMedia, time) %>% 
  mutate(lagOD5 = lag(avgOD, n=5),
         lagOD4 = lag(avgOD, n=4),
         lagOD3 = lag(avgOD, n=3),
         lagOD2 = lag(avgOD, n=2),
         lagOD = lag(avgOD),
         leadOD = lead(avgOD),
         leadOD2 = lead(avgOD, n=2),
         leadOD3 = lead(avgOD, n=3),
         leadOD4 = lead(avgOD, n=4),
         leadOD5 = lead(avgOD, n=5),
         lagMax5 = lag(maxOD, n=5),
         lagMax4 = lag(maxOD, n=4),
         lagMax3 = lag(maxOD, n=3),
         lagMax2 = lag(maxOD, n=2),
         lagMax = lag(maxOD),
         leadMax = lead(maxOD),
         leadMax2 = lead(maxOD, n=2),
         leadMax3 = lead(maxOD, n=3),
         leadMax4 = lead(maxOD, n=4),
         leadMax5 = lead(maxOD, n=5),
         lagMin5 = lag(minOD, n=5),
         lagMin4 = lag(minOD, n=4),
         lagMin3 = lag(minOD, n=3),
         lagMin2 = lag(minOD, n=2),
         lagMin = lag(minOD),
         leadMin = lead(minOD),
         leadMin2 = lead(minOD, n=2),
         leadMin3 = lead(minOD, n=3),
         leadMin4 = lead(minOD, n=4),
         leadMin5 = lead(minOD, n=5)) %>% 
  group_by(spentMedia, time) %>% 
  mutate(avgSmoothedOD = mean(c(lagOD5,lagOD4,lagOD3,lagOD2,lagOD,avgOD,leadOD,leadOD2,leadOD3,leadOD4,leadOD5), na.rm=TRUE),
         avgSmoothedMax = mean(c(lagMax5,lagMax4,lagMax3,lagMax2,lagMax,maxOD,leadMax,leadMax2,leadMax3,leadMax4,leadMax5), na.rm=TRUE),
         avgSmoothedMin = mean(c(lagMin5,lagMin4,lagMin3,lagMin2,lagMin,minOD,leadMin,leadMin2,leadMin3,leadMin4,leadMin5), na.rm=TRUE)) %>%
  ungroup() %>% 
  group_by(spentMedia) %>% 
  mutate(maxInMedia = max(avgSmoothedOD))
p_faecalisLactis <- faecalisLactisExample %>% 
  ggplot() +
  geom_line(aes(x=time, y=log10(avgSmoothedOD), group=spentMedia, color=spentMedia), size=0.4) +
  geom_ribbon(aes(x=time, ymin=log10(avgSmoothedMin), ymax=log10(avgSmoothedMax), group=spentMedia, fill=spentMedia), alpha=0.1) +
  geom_hline(aes(yintercept=log10(unique(maxInMedia)[1])), linetype="dotted", color="black", size=0.4) +
  geom_hline(aes(yintercept=log10(unique(maxInMedia)[2])), linetype="dotted", color="#BE1E2D", size=0.4) +
  scale_color_manual(name="Medium", values = c("black", "#BE1E2D")) +
  scale_fill_manual(name="Medium", values = c("black", "#BE1E2D"), guide="none") +
  xlab("Time (min)") +
  scale_y_continuous(name="OD600",
                     limits=c(-2.5,0), breaks=c(-2,-1,0),
                     labels=label_math(10^.x)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  theme(legend.key.size=unit(0.5,"lines")) +
  #theme(legend.position="none") +
  facet_wrap(~strain) +
  DEFAULTS.THEME_PRINT
p_faecalisLactis
save_plot("out/faecalisLactisExample.pdf", p_faecalisLactis, base_width=1.2, base_height=1.4)

faecalisLactisLegend <- get_legend(p_faecalisLactis)
save_plot("out/faecalisLactisExampleLegend.pdf", faecalisLactisLegend, base_width=0.75, base_height=0.5)

# Plot growth curves combined on same plots.
p_allCurvesCombined <- dataCombined %>% 
  filter(time <= 1440) %>% 
  ungroup() %>% 
  filter(strain %in% c("B. fragilis","P. goldsteinii","E. faecalis",
                       "E. casseliflavus","L. garvieae","L. lactis") &
           spentMedia %in% c("B. fragilis","P. goldsteinii","E. faecalis",
                             "E. casseliflavus","L. garvieae","L. lactis","Community D1","Community D2","fresh")) %>%
  mutate(spentMedia = ifelse(spentMedia=="fresh","Fresh mBHI", spentMedia)) %>% 
  group_by(spentMedia, strain, time) %>% 
  summarize(avgOD = mean(blankSubtractedOD),
            maxOD = max(blankSubtractedOD),
            minOD = min(blankSubtractedOD)) %>% 
  # Smooth curves over 1/2 hr intervals.
  group_by(strain, spentMedia) %>% 
  arrange(strain, spentMedia, time) %>% 
  mutate(lagOD5 = lag(avgOD, n=5),
         lagOD4 = lag(avgOD, n=4),
         lagOD3 = lag(avgOD, n=3),
         lagOD2 = lag(avgOD, n=2),
         lagOD = lag(avgOD),
         leadOD = lead(avgOD),
         leadOD2 = lead(avgOD, n=2),
         leadOD3 = lead(avgOD, n=3),
         leadOD4 = lead(avgOD, n=4),
         leadOD5 = lead(avgOD, n=5),
         lagMax5 = lag(maxOD, n=5),
         lagMax4 = lag(maxOD, n=4),
         lagMax3 = lag(maxOD, n=3),
         lagMax2 = lag(maxOD, n=2),
         lagMax = lag(maxOD),
         leadMax = lead(maxOD),
         leadMax2 = lead(maxOD, n=2),
         leadMax3 = lead(maxOD, n=3),
         leadMax4 = lead(maxOD, n=4),
         leadMax5 = lead(maxOD, n=5),
         lagMin5 = lag(minOD, n=5),
         lagMin4 = lag(minOD, n=4),
         lagMin3 = lag(minOD, n=3),
         lagMin2 = lag(minOD, n=2),
         lagMin = lag(minOD),
         leadMin = lead(minOD),
         leadMin2 = lead(minOD, n=2),
         leadMin3 = lead(minOD, n=3),
         leadMin4 = lead(minOD, n=4),
         leadMin5 = lead(minOD, n=5)) %>% 
  group_by(strain, spentMedia, time) %>% 
  mutate(avgSmoothedOD = mean(c(lagOD5,lagOD4,lagOD3,lagOD2,lagOD,avgOD,leadOD,leadOD2,leadOD3,leadOD4,leadOD5), na.rm=TRUE),
         avgSmoothedMax = mean(c(lagMax5,lagMax4,lagMax3,lagMax2,lagMax,maxOD,leadMax,leadMax2,leadMax3,leadMax4,leadMax5), na.rm=TRUE),
         avgSmoothedMin = mean(c(lagMin5,lagMin4,lagMin3,lagMin2,lagMin,minOD,leadMin,leadMin2,leadMin3,leadMin4,leadMin5), na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(strain = fct_relevel(strain, rev(c("P. goldsteinii","B. fragilis","L. lactis","L. garvieae",
                                            "E. casseliflavus","E. faecalis"))),
         spentMedia = fct_relevel(spentMedia, c("Fresh mBHI","E. faecalis","E. casseliflavus","L. garvieae","L. lactis",
                                                "B. fragilis","P. goldsteinii","Community D1","Community D2"))) %>% 
  ggplot() +
  geom_line(aes(x=time, y=log10(avgSmoothedOD), group=spentMedia, color=spentMedia), size=0.4, alpha=0.6) +
  xlab("Time (min)") +
  scale_y_continuous(name="OD600",
                     limits=c(-2.5,0), breaks=c(-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name="Medium",
                     values=c("Fresh mBHI" = "black", 
                              "B. fragilis" = "#762D5D",
                              "P. goldsteinii" = "#B62025",
                              "E. faecalis"="#C1A92F",
                              "E. casseliflavus"="#BF732E",
                              "L. garvieae"="#6398A4",
                              "L. lactis"="#343795",
                              "Community D1"="#295135",
                              "Community D2"="#8CB369")) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        strip.text = element_text(face="bold.italic")) +
  theme(legend.key.size=unit(0.5,"lines")) +
  facet_grid(~strain) +
  DEFAULTS.THEME_PRINT
p_allCurvesCombined
save_plot("out/allGrowthCurvesCombined.pdf", p_allCurvesCombined, base_width=7, base_height=2)

# Plot niche overlaps.
p_nicheOverlapTiles <- nicheOverlapsSeparateReps %>% 
  ungroup() %>% 
  mutate(strain = fct_relevel(strain, c("P. goldsteinii","B. fragilis","L. lactis","L. garvieae",
                                        "E. casseliflavus","E. faecalis")),
         spentMedia = fct_relevel(spentMedia, c("E. faecalis","E. casseliflavus","L. garvieae","L. lactis",
                                                "B. fragilis","P. goldsteinii","Community D1","Community D2")),
         textColor = ifelse(nicheOverlap>0.5, TRUE, FALSE)) %>% 
  ggplot() +
  geom_tile(aes(x=spentMedia, y=strain, fill=nicheOverlap), color="grey30", size=0.01) +
  geom_text(aes(x=spentMedia, y=strain, label=round(nicheOverlap,2), color=textColor), size=2) +
  scale_color_manual(values=c("white","black"), guide="none") +
  xlab("\n\nSpent media") +
  ylab("Strain\n\n") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  theme(legend.key.size=unit(0.65,"lines")) +
  #theme(legend.position="bottom") +
  theme(legend.position="none") +
  scale_fill_viridis(name = "Niche overlap", option="viridis",
                     limits = c(0,1), breaks=c(0,0.25,0.5,0.75,1),
                     labels = c(0,0.25,0.5,0.75,1)) +
  DEFAULTS.THEME_PRINT
p_nicheOverlapTiles
save_plot("out/nicheOverlapTiles.png", p_nicheOverlapTiles, base_width=2.8, base_height=2.4)
save_plot("out/nicheOverlapTiles.pdf", p_nicheOverlapTiles, base_width=3.2, base_height=2.6)

nicheOverlapLegend <- get_legend(p_nicheOverlapTiles)
save_plot("out/nicheOverlapLegendHorizontal.pdf", nicheOverlapLegend, base_width=1.5, base_height=0.5)

# Plot distribution of niche overlaps within vs across orders.
p_orderOverlapComparison <- nicheOverlapsSeparateReps %>% 
  filter(spentMedia %in% c("B. fragilis","P. goldsteinii","E. faecalis","E. casseliflavus","L. garvieae","L. lactis")) %>% 
  mutate(strainOrder = ifelse(strain %in% c("B. fragilis","P. goldsteinii"),"Bacteroidales","Lactobacillales"),
         mediaOrder = ifelse(spentMedia %in% c("B. fragilis","P. goldsteinii"),"Bacteroidales","Lactobacillales"),
         diffOrder = ifelse(strainOrder==mediaOrder, "Same order", "Different order"),
         diffOrder = fct_relevel(diffOrder, c("Same order", "Different order"))) %>% 
  ggplot() +
  geom_jitter(aes(x=diffOrder, y=nicheOverlap), width=0.15, size=0.25) +
  geom_boxplot(aes(x=diffOrder, y=nicheOverlap), alpha=0, size=0.25) +
  xlab("Comparison") +
  ylab("Niche overlap") +
  #ylim(0,1) +
  DEFAULTS.THEME_PRINT
p_orderOverlapComparison
save_plot("out/orderOverlapComparison.png", p_orderOverlapComparison)
save_plot("out/orderOverlapComparison.pdf", p_orderOverlapComparison, base_width=1.75, base_height=1.5)

orderOverlapCorTest <- wilcox.test(nicheOverlapsSeparateReps %>% 
                                     filter(spentMedia %in% c("B. fragilis","P. goldsteinii","E. faecalis","E. casseliflavus","L. garvieae","L. lactis")) %>% 
                                     mutate(strainOrder = ifelse(strain %in% c("B. fragilis","P. goldsteinii"),"Bacteroidales","Lactobacillales"),
                                            mediaOrder = ifelse(spentMedia %in% c("B. fragilis","P. goldsteinii"),"Bacteroidales","Lactobacillales"),
                                            diffOrder = ifelse(strainOrder==mediaOrder, "Same order", "Different order")) %>%
                                     filter(diffOrder=="Different order") %>%
                                     pull(nicheOverlap),
                                   nicheOverlapsSeparateReps %>% 
                                     filter(spentMedia %in% c("B. fragilis","P. goldsteinii","E. faecalis","E. casseliflavus","L. garvieae","L. lactis")) %>% 
                                     mutate(strainOrder = ifelse(strain %in% c("B. fragilis","P. goldsteinii"),"Bacteroidales","Lactobacillales"),
                                            mediaOrder = ifelse(spentMedia %in% c("B. fragilis","P. goldsteinii"),"Bacteroidales","Lactobacillales"),
                                            diffOrder = ifelse(strainOrder==mediaOrder, "Same order", "Different order")) %>%
                                     filter(diffOrder=="Same order") %>%
                                     pull(nicheOverlap),
                                   alternative="less")

# Compare correlation of spent media vs metabolomics niche overlaps. -------------------------------------------------------------------------


# Import strain overlaps from metabolomics.
metabolomicNicheOverlap <- fread("../../analysis-DG-clean/strainOverlaps.txt") %>% 
  ungroup() %>% 
  select(comparison, overlapAsymmetric) %>% 
  unique() %>% 
  mutate(strain1 = sub("-.*","",comparison),
         strain2 = sub(".*-","",comparison)) %>% 
  select(!comparison)

metabolomicNicheOverlapReplicates <- fread("../../analysis-DG-clean/replicateOverlaps.txt") %>% 
  select(strain, overlapAsymmetric) %>% 
  rename(strain1 = strain) %>% 
  mutate(strain2 = strain1)

nicheOverlapCombined <- rbind(
  metabolomicNicheOverlap, metabolomicNicheOverlapReplicates) %>%
  filter(!(strain1 %in% c("XFA","XFB"))) %>% 
  mutate(strain1 = case_when(strain1=="Bacte0126" ~ "B. fragilis",
                             strain1=="Tanne0007" ~ "P. goldsteinii",
                             strain1=="EnteC0002" ~ "E. faecalis",
                             strain1=="EnteC0003" ~ "E. casseliflavus",
                             strain1=="Strep0007" ~ "L. garvieae",
                             strain1=="Strep0017" ~ "L. lactis"),
         strain2 = case_when(strain2=="Bacte0126" ~ "B. fragilis",
                             strain2=="Tanne0007" ~ "P. goldsteinii",
                             strain2=="EnteC0002" ~ "E. faecalis",
                             strain2=="EnteC0003" ~ "E. casseliflavus",
                             strain2=="Strep0007" ~ "L. garvieae",
                             strain2=="Strep0017" ~ "L. lactis",
                             strain2=="XFA" ~ "Community D1",
                             strain2=="XFB" ~ "Community D2")) %>% 
  rename(metabolomicOverlap = overlapAsymmetric) %>% 
  left_join(nicheOverlapsSeparateReps %>% 
              select(!foldChange) %>% 
              rename(strain1 = strain, strain2 = spentMedia, spentMediaOverlap = nicheOverlap)) %>% 
  mutate(strain1 = fct_relevel(strain1, c("B. fragilis","P. goldsteinii","E. faecalis","E. casseliflavus",
                                          "L. garvieae","L. lactis")),
         strain2 = fct_relevel(strain2, c("B. fragilis","P. goldsteinii","E. faecalis","E. casseliflavus",
                                          "L. garvieae","L. lactis","Community D1","Community D2")))

p_nicheOverlapMethodsComparison <- nicheOverlapCombined %>% 
  ggplot() +
  #geom_point(aes(x=spentMediaOverlap, y=metabolomicOverlap, color=strain2), size=1, pch=1) +
  geom_abline(aes(slope=1, intercept=0), linetype="dashed", color="black", size=0.5) +
  geom_point(aes(x=spentMediaOverlap, y=metabolomicOverlap, fill=strain1, color=strain2), size=1, pch=21) +
  #geom_smooth(aes(x=spentMediaOverlap, y=metabolomicOverlap), method="lm") +
  xlim(0,1) +
  ylim(0,1) +
  xlab("Niche overlap, spent media") +
  ylab("Niche overlap, metabolomics") +
  scale_fill_manual(name="Focal strain",
                    values=c("B. fragilis" = "#762D5D",
                             "P. goldsteinii" = "#B62025",
                             "E. faecalis"="#C1A92F",
                             "E. casseliflavus"="#BF732E",
                             "L. garvieae"="#6398A4",
                             "L. lactis"="#343795"),
                    guide="none") +
  scale_color_manual(name="Medium/comparison strain",
                     values=c("B. fragilis" = "#762D5D",
                              "P. goldsteinii" = "#B62025",
                              "E. faecalis"="#C1A92F",
                              "E. casseliflavus"="#BF732E",
                              "L. garvieae"="#6398A4",
                              "L. lactis"="#343795",
                              "Community D1"="#295135",
                              "Community D2"="#8CB369")) +
  theme(legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  DEFAULTS.THEME_PRINT
p_nicheOverlapMethodsComparison
save_plot("out/nicheOverlapMethodsComparison.pdf", p_nicheOverlapMethodsComparison, base_width=2, base_height=2)

focalStrainLegend <- get_legend(p_nicheOverlapMethodsComparison)
comparisonStrainLegend <- get_legend(p_nicheOverlapMethodsComparison)

save_plot("out/focalStrainLegendScatter.pdf", focalStrainLegend, base_width=0.75, base_height=0.75)
save_plot("out/comparisonStrainLegendScatter.pdf", comparisonStrainLegend, base_width=1.2, base_height=1)

# Check correlation of methods.
methodsComparisonCor <- cor.test(nicheOverlapCombined$metabolomicOverlap, 
                                 nicheOverlapCombined$spentMediaOverlap, 
                                 method="pearson")
