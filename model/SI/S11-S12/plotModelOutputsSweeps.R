library(R.matlab)
library(tidyverse)
library(cowplot)
library(scales)
library(viridis)

#setwd("C:/Users/kxue9/Dropbox/Dose dependence/model-KXtest")

# Import plot defaults.
source("../palettes/plotDefaults.R")
theme_set(theme_cowplot())

ratios = c("1000:1","100:1","10:1","1:1","1:10","1:100","1:1000")
doseList <- c("1:1000","1000:1")

VERSION <- "m6-normalized-longer"

# Analyze model versions with varying consumption rates. -------------------------------------------------

data <- readMat(paste0("data/",VERSION,"/relAbundance-",VERSION,".mat"))
dfData <- reshape2::melt(data)
colnames(dfData) <- c("r11","r22","rs1","rs2","idose","passage","species","relAbundance","df")

params <- readMat(paste0("data/",VERSION,"/params-",VERSION,".mat"))
dfParams <- reshape2::melt(params) %>% 
  select(!c(Var2, L1)) %>% 
  mutate(Var1 = case_when(Var1==1 ~ "r11",
                          Var1==2 ~ "r22",
                          Var1==3 ~ "rs1",
                          Var1==4 ~ "rs2"))
r11List <- unique(dfParams %>% filter(Var1=="r11") %>% pull(value))
r22List <- unique(dfParams %>% filter(Var1=="r22") %>% pull(value))
rs1List <- unique(dfParams %>% filter(Var1=="rs1") %>% pull(value))
rs2List <- unique(dfParams %>% filter(Var1=="rs2") %>% pull(value))

dfData <- dfData %>% 
  mutate(r11 = r11List[r11],
         r22 = r22List[r22],
         rs1 = rs1List[rs1],
         rs2 = rs2List[rs2],
         dose = doseList[idose])

# Filter to only species that are mutually invasible - do not always decrease in abundance.
dfDataFiltered <- dfData %>% 
  group_by(r11, r22, rs1, rs2, species, dose) %>% 
  mutate(decreasing = relAbundance[passage==5] < relAbundance[passage==1]) %>% 
  ungroup() %>% 
  group_by(r11, r22, rs1, rs2, species) %>% 
  mutate(increasingOnce = FALSE %in% decreasing) %>% 
  ungroup() %>% 
  group_by(r11, r22, rs1, rs2) %>%
  mutate(mutuallyInvasible = !(FALSE %in% increasingOnce)) %>% 
  #filter(!(FALSE %in% increasingOnce)) %>% 
  select(-c(decreasing, increasingOnce))

# Filter out species that go below 10^-3 relative abundance at any point in 5 passages.
dfDataFiltered <- dfDataFiltered %>% 
  filter(species!=3) %>% 
  ungroup() %>% 
  mutate(lowAbundance = relAbundance < 10^-3) %>% 
  group_by(r11, r22, rs1, rs2) %>% 
  mutate(lowAbundance = TRUE %in% lowAbundance)

dfDoseDependence <- dfDataFiltered %>% 
  select(-idose) %>% 
  pivot_wider(names_from=dose, values_from=relAbundance) %>%
  mutate(doseDependence=ifelse(species==1, `1000:1`/`1:1000`,`1:1000`/`1000:1`),
         uniqueRatio = r22/r11,
         sharedRatio = rs2/rs1)

write.table(dfDoseDependence, paste0("data/",VERSION,"/dfDoseDependence-",VERSION,".txt"), 
            quote=FALSE, row.names=FALSE, sep="\t")

# Plot distribution of dose dependence at each passage.
p_ddDistribution <- dfDoseDependence %>% 
  #filter(mutuallyInvasible & !lowAbundance) %>% 
  ggplot() +
  geom_histogram(aes(x=log10(doseDependence)), bins=40) +
  facet_grid(~passage~species) +
  DEFAULTS.THEME_PRINT
p_ddDistribution
save_plot(paste0("out/",VERSION,"/ddDistribution-",VERSION,".png"),p_ddDistribution)

# Plot dose dependence over passages.
p_ddLines <- dfDoseDependence %>% 
  filter(mutuallyInvasible & !lowAbundance) %>% 
  #filter(r11=="0.05" & r22=="0.95" & rs1=="0.05" & rs2=="0.95") %>%
  #filter(r11=="0.05" & r22=="0.35" & rs1=="0.15" & rs2=="0.65") %>% 
  ggplot() +
  geom_line(aes(x=passage, y=log10(doseDependence), group=interaction(r11,r22,rs1,rs2))) +
  facet_wrap(~species) +
  DEFAULTS.THEME_PRINT
p_ddLines
save_plot(paste0("out/",VERSION,"/ddLinesLowAbundance-",VERSION,".png"),p_ddLines)

# Plot all four params as heatmap of DD.
p_ddHeatmapSp1 <- dfDoseDependence %>%
  #filter(as.character(r11) %in% c("0.2","0.5","0.8") & as.character(r22) %in% c("0.2","0.5","0.8")) %>% 
  filter(as.character(r11) %in% c("0.05","0.35","0.65","0.95") & as.character(r22) %in% c("0.05","0.35","0.65","0.95")) %>% 
  #filter(as.character(r11) %in% c("0.75","2.25","3.75","5.25") & as.character(r22) %in% c("0.75","2.25","3.75","5.25")) %>% 
  filter(species==1 & passage==5) %>% 
  mutate(highAbundance = mutuallyInvasible & !lowAbundance) %>% 
  rename(R1ab = rs1, R2ab = rs2) %>% 
  ggplot() +
  #geom_tile(aes(x=as.character(r11), y=as.character(r22), 
  #              fill=log10(doseDependence), alpha=highAbundance), color="grey30", size=0.01) +
  geom_tile(aes(x=as.character(r11), y=as.character(r22), 
                fill=log10(doseDependence)), color="grey30", size=0.01) +
  scale_fill_viridis(name = "Dose dependence\n(passage 5)", option="inferno",
                      limits = c(0,2.5), breaks=c(0,1,2),
                      labels = label_math(10^.x)) +
  #scale_fill_viridis(name = "Dose dependence\n(passage 5)", option="inferno",
  #                   limits = c(0,5), breaks=c(0,1,2,3,4,5),
  #                   labels = label_math(10^.x)) +
  #scale_alpha_manual(values = c(0,1), guide="none") +
  xlab("R1,a") +
  ylab("R2,b") +
  theme(legend.key.size=unit(0.5,"lines"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_grid(~R2ab~R1ab, labeller=label_both) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_ddHeatmapSp1
save_plot(paste0("out/",VERSION,"/ddHeatmapSp1Filtered-",VERSION,".png"),p_ddHeatmapSp1, nrow=3, ncol=3)
save_plot(paste0("out/",VERSION,"/ddHeatmapSp1Filtered-",VERSION,".pdf"),p_ddHeatmapSp1, 
          base_width=4.25, base_height=4.25)

# Get heatmap legend.
heatmapLegend <- get_legend(p_ddHeatmapSp1)
save_plot("out/ddHeatmapLegend-5.pdf", heatmapLegend, base_width=0.9, base_height=1)

p_ddHeatmapSp2 <- dfDoseDependence %>% 
  filter(species==2 & passage==5) %>% 
  mutate(highAbundance = mutuallyInvasible & !lowAbundance) %>% 
  rename(R1ab = rs1, R2ab = rs2) %>% 
  ggplot() +
  #geom_tile(aes(x=r11, y=r22, fill=log10(doseDependence))) +
  geom_tile(aes(x=as.character(r11), y=as.character(r22), 
                fill=log10(doseDependence), alpha=highAbundance), color="grey30", size=0.01) +
  #geom_tile(aes(x=as.character(r11), y=as.character(r22), 
  #              fill=log10(doseDependence)), color="grey30", size=0.01) +
  # scale_fill_viridis(name = "Dose dependence\n(passage 5)", option="inferno",
  #                    limits = c(0,2.5), breaks=c(0,1,2),
  #                    labels = label_math(10^.x)) +
  scale_fill_viridis(name = "Dose dependence\n(passage 5)", option="inferno",
                     limits = c(0,5), breaks=c(0,1,2,3,4,5),
                     labels = label_math(10^.x)) +
  scale_alpha_manual(values = c(0,1), guide="none") +
  xlab("R1,a") +
  ylab("R2,b") +
  theme(legend.key.size=unit(0.5,"lines"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_grid(~R2ab~R1ab, labeller=label_both) +
  #theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_ddHeatmapSp2
save_plot(paste0("out/",VERSION,"/ddHeatmapSp2-",VERSION,".png"),p_ddHeatmapSp2, nrow=3, ncol=3)

p_ddHeatmapRatiosSp1 <- dfDoseDependence %>% 
  filter(species==1 & passage==5) %>% 
  mutate(highAbundance = mutuallyInvasible & !lowAbundance) %>% 
  ggplot() +
  geom_point(aes(x=uniqueRatio, y=sharedRatio, color=log10(doseDependence), alpha=highAbundance)) +
  scale_color_viridis(option="inferno") +
  xlab("Unique ratio r22/r11") +
  ylab("Shared ratio rs2/rs1") +
  DEFAULTS.THEME_PRINT
p_ddHeatmapRatiosSp1
save_plot(paste0("out/",VERSION,"/ddHeatmapRatiosSp1-",VERSION,".png"),p_ddHeatmapRatiosSp1, nrow=3, ncol=3)

p_ddHeatmapRatiosSp2 <- dfDoseDependence %>% 
  filter(species==2 & passage==5) %>% 
  mutate(highAbundance = mutuallyInvasible & !lowAbundance) %>% 
  ggplot() +
  geom_point(aes(x=uniqueRatio, y=sharedRatio, color=log10(doseDependence), alpha=highAbundance)) +
  scale_color_viridis(option="inferno") +
  xlab("Unique ratio r22/r11") +
  ylab("Shared ratio rs2/rs1") +
  DEFAULTS.THEME_PRINT
p_ddHeatmapRatiosSp2
save_plot(paste0("out/",VERSION,"/ddHeatmapRatiosSp2-",VERSION,".png"),p_ddHeatmapRatiosSp2, nrow=3, ncol=3)

# Plot a specific relative abundance example.
p_relAb <- dfDoseDependence %>% 
  ungroup() %>% 
  #filter(mutuallyInvasible) %>% 
  #filter(uniqueRatio==19 & sharedRatio==19)
  #filter(uniqueRatio==19 & sharedRatio==15) %>% 
  #filter(uniqueRatio==7 & sharedRatio==9) %>% 
  #filter(uniqueRatio==7 & round(sharedRatio,2)==4.33) %>% 
  #filter(round(uniqueRatio,2)==6.33 & round(sharedRatio,2)==4.33) %>% 
  #filter(round(uniqueRatio,2)==0.05 & round(sharedRatio,2)==0.06) %>% 
  #filter(uniqueRatio==5 & sharedRatio==17) %>% 
  #filter(uniqueRatio==5 & round(sharedRatio,2)==1.12) %>% 
  #filter(r11=="0.8" & r22=="0.8" & rs1=="0.95" & rs2=="0.95") %>% # Example 1.
  #filter(r11=="0.8" & r22=="0.8" & rs1=="0.8" & rs2=="0.95") %>%   # Example 2.
  #filter(r11=="0.8" & r22=="0.8" & rs1=="0.35" & rs2=="0.5") %>% # Example 3.
  #filter(r11=="0.95" & r22=="0.95" & rs1=="0.95" & rs2=="0.95") %>%  # Example 1.
  #filter(r11=="0.95" & r22=="0.95" & rs1=="0.8" & rs2=="0.95") %>%  # Example 2.
  filter(r11=="0.95" & r22=="0.95" & rs1=="0.35" & rs2=="0.5") %>%  # Example 3.
  #filter(r11=="0.05" & r22=="0.95" & rs1=="0.95" & rs2=="0.95") %>% 
  select(passage, species, `1:1000`, `1000:1`) %>% 
  pivot_longer(c(`1000:1`,`1:1000`),names_to = "ratio") %>% 
  mutate(species = ifelse(species==1, "Species 1", "Species 2")) %>% 
  mutate(ratio = fct_relevel(ratio, c("1000:1","1:1000"))) %>% 
  ggplot() +
  geom_point(aes(x=ratio, y=log10(value), color=as.character(passage)), alpha=0.8, size=0.5) +
  geom_line(aes(x=ratio, y=log10(value), color=as.character(passage), group=passage), alpha=0.6, size=0.4) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-3,0), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name = "Passage", 
                     values = c("#3BBFF7","#2B84AE","#1B4965","#222847","#290628")) +
  theme(legend.key.size=unit(0.5,"lines")) +
  facet_wrap(~species) +
  DEFAULTS.THEME_PRINT
p_relAb
save_plot(paste0("out/",VERSION,"/relAbExampleLowAbundance-",VERSION,".png"),p_relAb)
save_plot(paste0("out/",VERSION,"/relAbExample3",VERSION,".pdf"), p_relAb, base_width=2.5, base_height=1.5)


# Analyze model versions with varying resource levels. --------------------


data <- readMat(paste0("data/",VERSION,"/relAbundance-",VERSION,".mat"))
dfData <- reshape2::melt(data)
colnames(dfData) <- c("Y1","Y2","Y3","idose","passage","species","relAbundance","df")

params <- readMat(paste0("data/",VERSION,"/params-",VERSION,".mat"))
dfParams <- reshape2::melt(params) %>% 
  select(!c(Var2, L1)) %>% 
  mutate(Var1 = case_when(Var1==1 ~ "Y1",
                          Var1==2 ~ "Y2",
                          Var1==3 ~ "Y3"))
Y1List <- unique(dfParams %>% filter(Var1=="Y1") %>% pull(value))
Y2List <- unique(dfParams %>% filter(Var1=="Y2") %>% pull(value))
Y3List <- unique(dfParams %>% filter(Var1=="Y3") %>% pull(value))

dfData <- dfData %>% 
  mutate(Y1 = Y1List[Y1],
         Y2 = Y2List[Y2],
         Y3 = Y3List[Y3],
         dose = doseList[idose])

dfDoseDependence <- dfData %>% 
  ## FIX THIS LATER
  filter(Y2==1 & Y3==1) %>% 
  select(-idose) %>% 
  pivot_wider(names_from=dose, values_from=relAbundance) %>%
  mutate(doseDependence=ifelse(species==1, `1000:1`/`1:1000`,`1:1000`/`1000:1`),
         gamma=Y2/(Y1+Y2))

# Plot dose dependence across Y1 values.
p_DDLines <- dfDoseDependence %>% 
  filter(species!=3) %>% 
  ggplot() +
  geom_line(aes(x=gamma, y=log10(doseDependence), group=passage, color=as.character(passage))) +
  #ylim(0,1) +
  xlab("y_a,b") +
  ylab("Log10 dose dependence") +
  ggtitle("ab=1") +
  facet_wrap(~species)
p_DDLines
save_plot(paste0("out/",VERSION,"/DDLines-",VERSION,".png"), p_DDLines)

p_relAb <- dfDoseDependence %>% 
  ungroup() %>% 
  filter(Y1==1.1 & species!=3) %>%  
  select(passage, species, `1:1000`, `1000:1`) %>% 
  pivot_longer(c(`1000:1`,`1:1000`),names_to = "ratio") %>% 
  mutate(species = ifelse(species==1, "Species 1", "Species 2")) %>% 
  mutate(ratio = fct_relevel(ratio, c("1000:1","1:1000"))) %>% 
  ggplot() +
  geom_point(aes(x=ratio, y=log10(value), color=as.character(passage)), alpha=0.8, size=0.5) +
  geom_line(aes(x=ratio, y=log10(value), color=as.character(passage), group=passage), alpha=0.6, size=0.4) +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-2,0), breaks=c(-2,-1,0),
                     labels=label_math(10^.x)) +
  scale_color_manual(name = "Passage", 
                     values = c("#3BBFF7","#2B84AE","#1B4965","#222847","#290628")) +
  theme(legend.key.size=unit(0.5,"lines")) +
  ggtitle("Y1 = 1.1, y_a,b = .48")+
  facet_wrap(~species) +
  DEFAULTS.THEME_PRINT
p_relAb
save_plot(paste0("out/",VERSION,"/relAbExample-equalY1",VERSION,".png"),p_relAb)
save_plot(paste0("out/",VERSION,"/relAbExample",VERSION,".pdf"), p_relAb, base_width=2.5, base_height=1.5)

# P5 relative abundance at max gamma (.909), when Y1=0.1: 
# P5 relative abundance at med gamma (0.18), when Y1=4.6: 
# P5 relative abundance at min gamma (0.09), when Y1=9.6: 


# Analyze model version to recapitulate main text model. ------------------

data <- readMat(paste0("data/",VERSION,"/relAbundance-",VERSION,".mat"))
dfData <- reshape2::melt(data)
colnames(dfData) <- c("gamma","idose","passage","species","relAbundance","df")

params <- readMat(paste0("data/",VERSION,"/params-",VERSION,".mat"))
dfParams <- reshape2::melt(params) %>% 
  select(!c(Var2, L1)) %>% 
  mutate(Var1 = case_when(Var1==1 ~ "gamma"))
gammaList <- unique(dfParams %>% filter(Var1=="gamma") %>% pull(value))

dfData <- dfData %>% 
  mutate(gamma = gammaList[gamma],
         dose = doseList[idose])

# Filter to only species that are mutually invasible - do not always decrease in abundance.
# dfDataFiltered <- dfData %>% 
#   group_by(gamma, species, dose) %>% 
#   mutate(decreasing = relAbundance[passage==5] < relAbundance[passage==1]) %>% 
#   ungroup() %>% 
#   group_by(gamma, species) %>% 
#   mutate(increasingOnce = FALSE %in% decreasing) %>% 
#   ungroup() %>% 
#   group_by(gamma) %>%
#   mutate(mutuallyInvasible = !(FALSE %in% increasingOnce)) %>% 
#   #filter(!(FALSE %in% increasingOnce)) %>% 
#   select(-c(decreasing, increasingOnce))

# Filter out species that go below 10^-3 relative abundance at any point in 5 passages.
# dfDataFiltered <- dfDataFiltered %>% 
#   filter(species!=3) %>% 
#   ungroup() %>% 
#   mutate(lowAbundance = relAbundance < 10^-3) %>% 
#   group_by(r11, r22, rs1, rs2) %>% 
#   mutate(lowAbundance = TRUE %in% lowAbundance)

dfDoseDependence <- dfData %>% 
  select(-idose) %>% 
  pivot_wider(names_from=dose, values_from=relAbundance) %>%
  mutate(doseDependence=ifelse(species==1, `1000:1`/`1:1000`,`1:1000`/`1000:1`))

write.table(dfDoseDependence, paste0("data/",VERSION,"/dfDoseDependence-",VERSION,".txt"), 
            quote=FALSE, row.names=FALSE, sep="\t")

p_gammaDD <- dfDoseDependence %>%
  filter(species==1, passage %in% c(3,5)) %>%
  #filter(passage==3) %>% 
  ggplot() +
  #geom_hline(yintercept=10, linetype="dashed", size=0.4) +
  #geom_vline(xintercept=0.975, linetype="dashed", size=0.4) +
  geom_line(aes(x=gamma, y=log10(doseDependence), group=passage,
                linetype=factor(passage)), size=0.4) +
  geom_point(data=dfDoseDependence %>%
               filter(gamma %in% gammaList[c(0.75/0.025, 0.9/0.025, 0.975/0.025)],
                      #filter(gamma %in% gammaList[c(0.75/0.025)],    
                      #filter(gamma %in% gammaList[c(0.75/0.025, 0.9/0.025)],    
                      species==1, passage %in% c(3,5)),
             #species==1, passage==3),
             aes(x=gamma, y=log10(doseDependence), color=factor(gamma)), size=1, alpha=0.75) +
  scale_y_continuous(name="Dose dependence", limits=c(0,3), breaks = c(0,1,2,3), labels = label_math(10^.x)) +
  scale_x_continuous(name="Niche overlap (.)", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1), 
                     labels=c("0","0.25","0.5","0.75","1")) +
  scale_linetype_manual(values=c("dotted","solid"), name="Passage") +
  scale_color_manual(values = c("#77A440","#707FBE","#C33D54")) +
  theme(legend.position="none") +
  xlab("Niche overlap (.)") + 
  ylab("Dose dependence") +
  DEFAULTS.THEME_PRINT
p_gammaDD
