library(R.matlab)
library(tidyverse)
library(cowplot)
library(scales)
library(viridis)

# Import plot defaults.
source("../palettes/plotDefaults.R")
theme_set(theme_cowplot())

ratios = c("1000:1","100:1","10:1","1:1","1:10","1:100","1:1000")
doseList <- c("1:1000","1000:1")


# Vary rt1, rt3. ----------------------------------------------------------


VERSION <- "s1-m22"

data <- readMat(paste0("data/",VERSION,"/relAbundance-",VERSION,".mat"))
dfData <- reshape2::melt(data)
colnames(dfData) <- c("rt1","rt3","idose","passage","species","relAbundance","df")

params <- readMat(paste0("data/",VERSION,"/params-",VERSION,".mat"))
dfParams <- reshape2::melt(params) %>% 
  select(!c(Var2, L1)) %>% 
  mutate(Var1 = case_when(Var1==1 ~ "rt1",
                          Var1==2 ~ "rt3"))
rt1List <- unique(dfParams %>% filter(Var1=="rt1") %>% pull(value))
rt3List <- unique(dfParams %>% filter(Var1=="rt3") %>% pull(value))

dfData <- dfData %>% 
  mutate(rt1 = rt1List[rt1],
         rt3 = rt3List[rt3],
         dose = doseList[idose])

dfDoseDependence <- dfData %>% 
  select(-idose) %>% 
  pivot_wider(names_from=dose, values_from=relAbundance) %>%
  mutate(doseDependence=ifelse(species==1, `1000:1`/`1:1000`,`1:1000`/`1000:1`))

write.table(dfDoseDependence, paste0("data/",VERSION,"/dfDoseDependence-",VERSION,".txt"), 
            quote=FALSE, row.names=FALSE, sep="\t")

# Plot distribution of dose dependence at each passage.
p_ddDistribution <- dfDoseDependence %>% 
  ggplot() +
  geom_histogram(aes(x=log10(doseDependence)), bins=40) +
  facet_grid(~passage~species) +
  DEFAULTS.THEME_PRINT
p_ddDistribution
save_plot(paste0("out/",VERSION,"/ddDistribution-",VERSION,".png"),p_ddDistribution)

# Plot dose dependence over passages.
p_ddLines <- dfDoseDependence %>% 
  ggplot() +
  geom_line(aes(x=passage, y=log10(doseDependence), group=interaction(rt1,rt3))) +
  facet_wrap(~species) +
  DEFAULTS.THEME_PRINT
p_ddLines
save_plot(paste0("out/",VERSION,"/ddLines-",VERSION,".png"),p_ddLines)

# Plot rt1 vs rt3 as heatmap of DD.
p_ddHeatmap <- dfDoseDependence %>% 
  ggplot() +
  geom_tile(aes(x=rt3, y=rt1, fill=log10(doseDependence))) +
  scale_fill_viridis(option="inferno") +
  xlab("rt3") +
  ylab("rt1") +
  facet_grid(~passage~species) +
  DEFAULTS.THEME_PRINT
p_ddHeatmap
save_plot(paste0("out/",VERSION,"/ddHeatmap-",VERSION,".png"),p_ddHeatmap)

# Plot a specific relative abundance example.
p_relAb <- dfDoseDependence %>% 
  ungroup() %>% 
  select(passage, species, `1:1000`, `1000:1`) %>% 
  pivot_longer(c(`1000:1`,`1:1000`),names_to = "ratio") %>% 
  mutate(ratio = fct_relevel(ratio, c("1000:1","1:1000"))) %>% 
  ggplot() +
  geom_point(aes(x=ratio, y=log10(value), color=as.character(passage))) +
  geom_line(aes(x=ratio, y=log10(value), color=as.character(passage), group=passage)) +
  facet_wrap(~species) +
  DEFAULTS.THEME_PRINT
p_relAb
save_plot(paste0("out/",VERSION,"/relAbExample-",VERSION,".png"),p_relAb)


# Vary rs1, rs2, rt1, rt3. ------------------------------------------------


VERSION <- "s2-m42-normalized-longer"

data <- readMat(paste0("data/",VERSION,"/relAbundance-",VERSION,".mat"))
dfData <- reshape2::melt(data)
colnames(dfData) <- c("rs1", "rs2", "rt1","rt3","idose","passage","species","relAbundance","df")

params <- readMat(paste0("data/",VERSION,"/params-",VERSION,".mat"))
dfParams <- reshape2::melt(params) %>% 
  select(!c(Var2, L1)) %>% 
  mutate(Var1 = case_when(Var1==1 ~ "rs1",
                          Var1==2 ~ "rs2",
                          Var1==3 ~ "rt1",
                          Var1==4 ~ "rt3"))
rs1List <- unique(dfParams %>% filter(Var1=="rs1") %>% pull(value))
rs2List <- unique(dfParams %>% filter(Var1=="rs2") %>% pull(value))
rt1List <- unique(dfParams %>% filter(Var1=="rt1") %>% pull(value))
rt3List <- unique(dfParams %>% filter(Var1=="rt3") %>% pull(value))

dfData <- dfData %>% 
  mutate(rs1 = rs1List[rs1],
         rs2 = rs2List[rs2],
         rt1 = rt1List[rt1],
         rt3 = rt3List[rt3],
         dose = doseList[idose])

# Filter to only species that are mutually invasible - do not always decrease in abundance.
dfDataFiltered <- dfData %>% 
  group_by(rs1, rs2, rt1, rt3, species, dose) %>% 
  mutate(decreasing = relAbundance[passage==5] < relAbundance[passage==1]) %>% 
  ungroup() %>% 
  group_by(rs1, rs2, rt1, rt3, species) %>% 
  mutate(increasingOnce = FALSE %in% decreasing) %>% 
  ungroup() %>% 
  group_by(rs1, rs2, rt1, rt3) %>%
  mutate(mutuallyInvasible = !(FALSE %in% increasingOnce)) %>% 
  #filter(!(FALSE %in% increasingOnce)) %>% 
  select(-c(decreasing, increasingOnce))

# Filter out species that go below 10^-3 relative abundance at any point in 5 passages.
dfDataFiltered <- dfDataFiltered %>% 
  ungroup() %>% 
  mutate(lowAbundance = relAbundance < 10^-3) %>% 
  group_by(rs1, rs2, rt1, rt3) %>% 
  mutate(lowAbundance = TRUE %in% lowAbundance)

dfDoseDependence <- dfDataFiltered %>% 
  select(-idose) %>% 
  pivot_wider(names_from=dose, values_from=relAbundance) %>%
  mutate(doseDependence=ifelse(species==1, `1000:1`/`1:1000`,`1:1000`/`1000:1`)) %>% 
  mutate(sharedSratio = rs2/rs1,
         sharedTratio = rt3/rt1)

write.table(dfDoseDependence, paste0("data/",VERSION,"/dfDoseDependence-",VERSION,".txt"), 
            quote=FALSE, row.names=FALSE, sep="\t")

# Plot distribution of dose dependence at each passage.
p_ddDistribution <- dfDoseDependence %>% 
  ggplot() +
  geom_histogram(aes(x=log10(doseDependence)), bins=40) +
  facet_grid(~passage~species) +
  DEFAULTS.THEME_PRINT
p_ddDistribution
save_plot(paste0("out/",VERSION,"/ddDistribution-",VERSION,".png"),p_ddDistribution)

# Plot dose dependence over passages.
p_ddLines <- dfDoseDependence %>% 
  ggplot() +
  geom_line(aes(x=passage, y=log10(doseDependence), group=interaction(rs1, rs2, rt1,rt3))) +
  facet_wrap(~species) +
  DEFAULTS.THEME_PRINT
p_ddLines
save_plot(paste0("out/",VERSION,"/ddLines-",VERSION,".png"),p_ddLines)

# Plot rt1 vs rt3 as heatmap of DD, species 1.
p_ddHeatmapSp1 <- dfDoseDependence %>%
  filter(as.character(rt1) %in% c("0.05","0.35","0.65","0.95") & 
           as.character(rt3) %in% c("0.05","0.35","0.65","0.95")) %>% 
  #filter(as.character(rt1) %in% c("0.2","0.5","0.8") & 
  #         as.character(rt3) %in% c("0.2","0.5","0.8")) %>% 
  #filter(rs1 %in% c("0.8","0.95") & rs2 %in% c("0.8","0.95") &
  #         rt1 %in% c("0.8","0.95") & rt3 %in% c("0.8","0.95")) %>% 
  filter(passage==5 & species==1) %>% 
  mutate(highAbundance = mutuallyInvasible & !lowAbundance) %>% 
  rename(R1ab = rs1, R2ab = rs2) %>% 
  ggplot() +
  #geom_tile(aes(x=as.character(rt1), y=as.character(rt3), fill=log10(doseDependence), alpha=highAbundance), 
  #          color="grey30", size=0.01) +
  geom_tile(aes(x=as.character(rt1), y=as.character(rt3), fill=log10(doseDependence)), 
            color="grey30", size=0.01) +
  scale_fill_viridis(name = "Dose dependence\n(passage 5)", option="inferno",
                     limits = c(0,2.5), breaks=c(0,0.5,1,1.5,2,2.5)) +
  #scale_fill_viridis(name = "Dose dependence\n(passage 5)", option="inferno",
  #                   limits = c(-0.1,5), breaks=c(0,1,2,3,4,5),
  #                   labels = label_math(10^.x)) +
  #scale_alpha_manual(values = c(0,1), guide="none") +
  xlab("R1ac") +
  ylab("R3ac") +
  theme(legend.key.size=unit(0.5,"lines"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  theme(legend.position = "none") +
  facet_grid(~R2ab~R1ab, labeller=label_both) +
  DEFAULTS.THEME_PRINT
p_ddHeatmapSp1
save_plot(paste0("out/",VERSION,"/ddHeatmapSp1-",VERSION,".png"),p_ddHeatmapSp1, nrow=3, ncol=3)
save_plot(paste0("out/",VERSION,"/ddHeatmapSp1-",VERSION,".pdf"), p_ddHeatmapSp1, 
          base_width=4.25, base_height=4.25)

heatmapLegend <- get_legend(p_ddHeatmapSp1)
save_plot("out/ddHeatmapLegend-3speciesHighGamma.pdf", heatmapLegend, base_width=0.9, base_height=1)

# Plot rt1 vs rt3 as heatmap of DD, species 2.
p_ddHeatmapSp2 <- dfDoseDependence %>% 
  filter(passage==5 & species==2) %>% 
  ggplot() +
  geom_tile(aes(x=rt3, y=rt1, fill=log10(doseDependence))) +
  scale_fill_viridis(option="inferno") +
  xlab("rt3") +
  ylab("rt1") +
  facet_grid(~rs2~rs1) +
  DEFAULTS.THEME_PRINT
p_ddHeatmapSp2
save_plot(paste0("out/",VERSION,"/ddHeatmapSp2-",VERSION,".png"),p_ddHeatmapSp2, nrow=3, ncol=3)

# Plot rt1 vs rt3 as heatmap of DD, species 3.
p_ddHeatmapSp3 <- dfDoseDependence %>% 
  filter(passage==5 & species==3) %>% 
  ggplot() +
  geom_tile(aes(x=rt3, y=rt1, fill=log10(doseDependence))) +
  scale_fill_viridis(option="inferno") +
  xlab("rt3") +
  ylab("rt1") +
  facet_grid(~rs2~rs1) +
  DEFAULTS.THEME_PRINT
p_ddHeatmapSp3
save_plot(paste0("out/",VERSION,"/ddHeatmapSp3-",VERSION,".png"),p_ddHeatmapSp3, nrow=3, ncol=3)

# Plot a specific relative abundance example.
p_relAb <- dfDoseDependence %>% 
  #filter(rs1==0.5 & rs2==0.5 & rt1==0.5 & rt3==0.5) %>% 
  #filter(rt1=="0.2" & rt3=="0.8" & rs1=="0.95" & rs2=="0.95") %>%
  #filter(rt1=="0.35" & rt3=="0.95" & rs1=="0.05" & rs2=="0.05") %>% 
  #filter(rt1=="0.5" & rt3=="0.8" & rs1=="0.8" & rs2=="0.8") %>% 
  #filter(rt1=="0.05" & rt3=="0.95" & rs1=="0.95" & rs2=="0.95") %>% # Example 1, symmetric model
  #filter(rt1=="0.95" & rt3=="0.95" & rs1=="0.05" & rs2=="0.95") %>% # Example 2, symmetric model
  filter(rt1=="0.35" & rt3=="0.65" & rs1=="0.2" & rs2=="0.5") %>% # Example 3, symmetric model
  #filter(rt1=="0.05" & rt3=="0.95" & rs1=="0.95" & rs2=="0.95") %>% # Example 1, low/high gamma.
  #filter(rt1=="0.05" & rt3=="0.95" & rs1=="0.8" & rs2=="0.95") %>% # Example 2, low/high gamma.
  #filter(rt1=="0.05" & rt3=="0.95" & rs1=="0.5" & rs2=="0.65") %>% # Example 3, low/high gamma.
  #filter(rt1=="0.95" & rt3=="0.95" & rs1=="0.8" & rs2=="0.95") %>% # New example 2, low/high gamma.
  #filter(rt1=="0.05" & rt3=="0.95" & rs1=="0.8" & rs2=="0.95") %>% # New example 3, low/high gamma.
  #filter(rt1=="7.5" & rt3=="7.5" & rs1=="7.5" & rs2=="7.5") %> 
  ungroup() %>% 
  select(passage, species, `1:1000`, `1000:1`) %>% 
  pivot_longer(c(`1000:1`,`1:1000`),names_to = "ratio") %>% 
  mutate(species = case_when(species==1 ~ "Species 1", 
                             species==2 ~ "Species 2",
                             species==3 ~ "Species 3")) %>%
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
save_plot(paste0("out/",VERSION,"/relAbExample-",VERSION,".png"),p_relAb)
save_plot(paste0("out/",VERSION,"/relAbExample2",VERSION,".pdf"), p_relAb, base_width=3.75, base_height=1.5)

# Plot DD against rs ratio and rt ratio.
p_ddHeatmap2D <- dfDoseDependence %>% 
  filter(passage==5 & species==1) %>% 
  ggplot() +
  geom_tile(aes(x=rt3, y=rt1, fill=log10(doseDependence))) +
  scale_fill_viridis(option="inferno") +
  facet_wrap(~factor(sharedSratio)) +
  DEFAULTS.THEME_PRINT
p_ddHeatmap2D
