library(R.matlab)
library(tidyverse)
library(cowplot)
library(scales)

#setwd("C:/Users/kxue9/Dropbox/Dose dependence/model-KXtest")

# Import plot defaults.
source("palettes/plotDefaults.R")
theme_set(theme_cowplot())

ratios = c("1000:1","100:1","10:1","1:1","1:10","1:100","1:1000")

# Dose dependence by gamma, with insets -------------------------------------------------------

# I converted the original model parameter alpha (unique resource / shared resource)
# to gamma (shared resource / (unique resource + shared resource))

# Import the data that relates dose dependence over passages to gamma.
relAbundanceByGamma <- readMat("relAbundanceByGamma-5resources.mat")
dfGamma <- reshape2::melt(relAbundanceByGamma)
colnames(dfGamma) <- c("igamma", "idose", "passage", "species", "relAbundance", "df")

# Import the list of values of gamma and dose.
gammaList <- seq(0.025,1,0.025)
#gammaList[40] <- 0.999
doseList <- c("1:1000","1:100","1:10","1:1","10:1","100:1","1000:1")
# Bind the values of gamma and dose to the matrix.
# Also calculate the original model parameter alpha based on gamma.
# Alpha is the ratio of the unique niche to the shared niche.
dfGamma <- dfGamma %>%
  mutate(gamma=gammaList[igamma], dose=doseList[idose],
         alpha=(1-gamma)/(2*gamma))

# Plots the species trajectories across dose 
# at individual values of gamma.
p_trajectoriesAtGamma <- dfGamma %>%
  filter(gamma %in% gammaList[c(0.75/0.025, 0.9/0.025, 0.975/0.025)], 
         species!=3, passage %in% c(3,5)) %>%
  #filter(gamma %in% gammaList[c(0.75/0.025, 0.9/0.025, 0.975/0.025)], 
  #       species!=3, passage==3) %>%
  mutate(plot = case_when(gamma == "0.75" ~ ".1 = 0.75",
                          gamma == "0.9" ~ ".1 = 0.9",
                          gamma == "0.975" ~ ".1 = 0.975")) %>% 
  ggplot() +
  geom_line(aes(x=fct_relevel(dose, ratios), y=log10(relAbundance), group=interaction(species, passage),
                linetype=factor(passage), color=factor(species)), alpha=0.6, size=0.4) +
  geom_point(aes(x=fct_relevel(dose, ratios), y=log10(relAbundance), shape=factor(passage), color=factor(species)), 
             alpha=0.8, size=0.5) +
  scale_color_manual(name="Species",
                     labels=c("1"="Species 1", "2"="Species 2"), 
                     values=c("1"="#943F92", "2"="#148989")) +
  scale_linetype_manual(labels=c("Passage 3","Passage 5"),
                        values=c("dotted","solid"), guide="none") +
  scale_shape_manual(values=c(1, 19), guide="none") +
  xlab("Mixture ratio") +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-1.5,0), breaks=c(-1,0),
                     labels=label_math(10^.x)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  facet_wrap(~plot) +
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_trajectoriesAtGamma
save_plot("figs2/trajectoriesAtGamma.pdf", p_trajectoriesAtGamma, base_width=2.8, base_height=1.5)
save_plot("figs2/trajectoriesAtGamma.png", p_trajectoriesAtGamma, base_width=2.8, base_height=1.5)

speciesLegend1and2 <- get_legend(p_trajectoriesAtGamma)
save_plot("figs/speciesLegend1and2.pdf", speciesLegend1and2, base_width=1.25, base_height=0.25)

dfGammaDoseDependence <- dfGamma %>%
  dplyr::select(-igamma, -idose) %>%
  pivot_wider(names_from=dose, values_from=relAbundance) %>%
  mutate(doseDependence=ifelse(species==1, `1000:1`/`1:1000`,`1:1000`/`1000:1`))
p_gammaDD <- dfGammaDoseDependence %>%
  filter(species==1, passage %in% c(3,5)) %>%
  #filter(passage==3) %>% 
  ggplot() +
  #geom_hline(yintercept=10, linetype="dashed", size=0.4) +
  #geom_vline(xintercept=0.975, linetype="dashed", size=0.4) +
  geom_line(aes(x=gamma, y=log10(doseDependence), group=passage,
                linetype=factor(passage)), size=0.4) +
  geom_point(data=dfGammaDoseDependence %>%
               filter(gamma %in% gammaList[c(0.75/0.025, 0.9/0.025, 0.975/0.025)],
               #filter(gamma %in% gammaList[c(0.75/0.025)],    
               #filter(gamma %in% gammaList[c(0.75/0.025, 0.9/0.025)],    
                      species==1, passage %in% c(3,5)),
                      #species==1, passage==3),
             aes(x=gamma, y=log10(doseDependence), color=factor(gamma)), size=1, alpha=0.75) +
  scale_y_continuous(name="Dose dependence", limits=c(0,3.01), breaks = c(0,1,2,3), labels = label_math(10^.x)) +
  scale_x_continuous(name="Niche overlap (.)", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1), 
                     labels=c("0","0.25","0.5","0.75","1")) +
  scale_linetype_manual(values=c("dotted","solid"), name="Passage") +
  scale_color_manual(values = c("#77A440","#707FBE","#C33D54")) +
  theme(legend.position="none") +
  xlab("Niche overlap (.)") + 
  ylab("Dose dependence") +
  DEFAULTS.THEME_PRINT
p_gammaDD
save_plot("figs2/gammaDD.pdf", p_gammaDD, base_width=1.4, base_height=1.2)
save_plot("figs2/gammaDD.png", p_gammaDD, base_width=1.4, base_height=1.2)

# Dose dependence by gamma, across beta -----------------------------------

# Import the data that relates dose dependence over passages to gamma.
relAbundanceByBetaGamma <- readMat("relAbundanceByBetaGammaAll-5resources.mat")
dfBetaGamma <- reshape2::melt(relAbundanceByBetaGamma)
colnames(dfBetaGamma) <- c("igamma","ibeta", "idose", "passage", "species", "relAbundance", "df")

# Import the list of values of gamma, beta, and dose.
gammaList <- seq(0.025,1,0.025)
#gammaList[40] <- 0.999
betaList <- c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.8, 0.85, 0.9, 0.975)
doseList <- c("1:1000","1000:1")

# Bind the values of beta, gamma, and dose to the matrix.
dfBetaGamma <- dfBetaGamma %>%
  mutate(beta=betaList[ibeta], gamma=gammaList[igamma], dose=doseList[idose])

dfBetaGammaDoseDependence <- dfBetaGamma %>%
  dplyr::select(-ibeta, -igamma, -idose) %>%
  pivot_wider(names_from=dose, values_from=relAbundance) %>% 
  mutate(doseDependence=ifelse(species==1, `1000:1`/`1:1000`,`1:1000`/`1000:1`))
  # Calculate dose dependence as logit.
  #mutate(doseDependence = ifelse(`1000:1` > `1:1000`, 
  #                               (`1000:1`/(1-`1:1000`))/(`1:1000`/(1-`1:1000`)), 
  #                               (`1:1000`/(1-`1000:1`))/(`1000:1`/(1-`1000:1`))))

p_gammaDDAtBeta <- dfBetaGammaDoseDependence %>%
  #filter(beta==0) %>% 
  #filter(beta %in% c(0,0.1)) %>%
  #filter(beta %in% c(0, 0.1, 0.5)) %>% 
  #filter(gamma!=0.999) %>% 
  filter(passage==5) %>% 
  filter(beta %in% c(0, 0.1, 0.5, 0.9)) %>% 
  #filter(passage %in% c(1,3,5)) %>% 
  #filter(beta!=0.975) %>% 
  #filter(passage==5, species %in% c(1,2,3),
  #       beta %in% c(0, 0.1, 0.5, 0.75, 0.85)) %>%
  mutate(plot = case_when(species==1 ~ "Species 1", 
                          species==2 ~ "Species 2",
                          species==3 ~ "Species 3")) %>% 
  ggplot() +
  geom_line(aes(x=gamma, y=log10(doseDependence), group=interaction(species, beta),
                color=factor(species), linetype=factor(beta)), size=0.4) +
  #scale_y_log10() +
  scale_x_continuous(name="Niche overlap, species 1 and 2 (.)", limits=c(0,1), 
                     breaks=c(0,0.25,0.5,0.75,1), labels=c("0","0.25","0.5","0.75","1")) +
  scale_y_continuous(name="Dose dependence", limits=c(-0.05,3.05), breaks = c(0,1,2,3), labels = label_math(10^.x)) +
  scale_linetype_manual(values=c("solid","dashed","dotted", "dotdash"), 
                        name="Niche overlap\nbetween species\n1 and 3 (B)") +
  scale_color_manual(name="Species",
                     labels=c("1"="Species 1", "2"="Species 2"), 
                     values=c("1"="#943F92", "2"="#148989"), guide="none") +
  xlab("Niche overlap, species 1 and 2 (.)") + 
  ylab("Dose dependence") +
  theme(legend.key.size=unit(0.5,"lines")) +
  facet_grid(~plot) +
  DEFAULTS.THEME_PRINT
p_gammaDDAtBeta
save_plot("figs/gammaDDAtBeta.pdf", p_gammaDDAtBeta, base_width=3.6, base_height=1.6)
save_plot("figs/gammaDDAtBeta.png", p_gammaDDAtBeta, base_width=3.6, base_height=1.6)

save_plot("figs2/gammaDDAtBeta.pdf", p_gammaDDAtBeta, base_width=4.5, base_height=1.6)
save_plot("figs2/gammaDDAtBeta.png", p_gammaDDAtBeta, base_width=4.5, base_height=1.6)
save_plot("figs2/gammaDDAtBeta.png", p_gammaDDAtBeta, base_width=6, base_height=4.5)

# Plot relative abundance trajectories.
p_relAb <- dfBetaGamma %>% 
  filter(passage %in% c(1,3,5)) %>% 
  #filter(gamma %in% c("0.25", "0.375", "0.5", "0.625", "0.75")) %>%
  filter(beta=="0.9") %>% 
  mutate(dose = fct_relevel(dose, c("1000:1","1:1000"))) %>% 
  ggplot() +
  geom_point(aes(x=dose, y=log10(relAbundance), color=gamma), alpha=0.8, size=0.5) +
  geom_line(aes(x=dose, y=log10(relAbundance), color=gamma, group=gamma), alpha=0.6, size=0.4) +
  xlab("Mixture ratio") +
  #scale_y_continuous(name="Relative abundance",
  #                   limits=c(-3,0), breaks=c(-3,-2,-1,0),
  #                   labels=label_math(10^.x)) +
  #scale_color_manual(name = "Passage", 
  #                   values = c("#3BBFF7","#2B84AE","#1B4965","#222847","#290628")) +
  theme(legend.key.size=unit(0.5,"lines")) +
  facet_grid(~passage~species) +
  #facet_grid(~gamma~species) +
  DEFAULTS.THEME_PRINT
p_relAb
save_plot("figs2/relAbBeta0.9.png", p_relAb, base_width=6, base_height=4.5)

p_relAbSp1 <- dfBetaGamma %>% 
  filter(passage==5 & species==1) %>% 
  mutate(dose = fct_relevel(dose, c("1000:1","1:1000"))) %>% 
  ggplot() + 
  geom_tile(aes(x=gamma, y=as.character(beta), fill=log10(relAbundance))) +
  #geom_text(aes(x=gamma, y=as.character(beta), label=round(log10(relAbundance), 1))) +
  facet_wrap(~dose)
p_relAbSp1

# Annotate mutual invasibility of each species.
dfBetaGammaAnnotated <- dfBetaGamma %>% 
  group_by(beta, gamma, species, dose) %>% 
  mutate(decreasing = relAbundance[passage==5] < relAbundance[passage==1]) %>% 
  ungroup() %>% 
  group_by(beta, gamma, species) %>% 
  mutate(increasingOnce = FALSE %in% decreasing) %>% 
  ungroup() %>% 
  group_by(beta, gamma) %>%
  mutate(mutuallyInvasible = !(FALSE %in% increasingOnce)) %>% 
  #filter(!(FALSE %in% increasingOnce)) %>% 
  select(-c(decreasing, increasingOnce))

# Dose dependence by gamma, across r21 (relative consumption of shared resource) -----------------------------------

# Import the data that relates dose dependence over passages to gamma.
relAbundanceByR21Gamma <- readMat("relAbundanceByR21GammaAll-5resources.mat")
dfR21Gamma <- reshape2::melt(relAbundanceByR21Gamma)
colnames(dfR21Gamma) <- c("igamma","ir21", "idose", "passage", "species", "relAbundance", "df")

# Import the list of values of gamma, beta, and dose.
gammaList <- seq(0.025,1,0.025)
r21List <- c(0.9, 1, 1.1)
doseList <- c("1:1000","1000:1")

# Bind the values of r21, gamma, and dose to the matrix.
dfR21Gamma <- dfR21Gamma %>%
  mutate(r21=r21List[ir21], gamma=gammaList[igamma], dose=doseList[idose])

dfR21GammaDoseDependence <- dfR21Gamma %>%
  dplyr::select(-ir21, -igamma, -idose) %>%
  pivot_wider(names_from=dose, values_from=relAbundance) %>%
  mutate(doseDependence=ifelse(species==1, `1000:1`/`1:1000`,`1:1000`/`1000:1`))

p_gammaDDatR21 <- dfR21GammaDoseDependence %>%
  #filter(passage==5, species %in% c(1)) %>%
  filter(passage==5 & species %in% c(1,2)) %>%
  mutate(plot = ifelse(species==1, "Species 1", "Species 2")) %>% 
  ggplot() +
  geom_line(aes(x=gamma, y=log10(doseDependence), group=interaction(species, r21),
                color=factor(r21)), size=0.4) +
  scale_color_manual(name = "Relative consumption\nrate of resource ab,\nspecies 1 and 2",
                     values = c("#9BC3D4","#8C70C2","#BB62BD")) +
  scale_y_continuous(name="Dose dependence", 
                     limits=c(0,4.5), breaks=c(0,1,2,3,4), 
                     labels=label_math(10^.x)) +
  scale_x_continuous(name="Niche overlap, species 1 and 2 (.)",
                     limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1),
                     labels=c("0","0.25","0.5","0.75","1")) +
  theme(legend.key.size=unit(0.5,"lines")) +
  #theme(legend.position="none") +
  facet_wrap(~plot) +
  DEFAULTS.THEME_PRINT
p_gammaDDatR21
save_plot("figs2/gammaDDAtR21BothSpecies.pdf", p_gammaDDatR21, base_width=3.6, base_height=1.6)

consumptionRateLegend <- get_legend(p_gammaDDatR21)
save_plot("figs/consumptionRateLegend.pdf", consumptionRateLegend, base_width=1, base_height=0.75)


# Dose dependence with asymmetric r21. ------------------------------------

# Import data.
relAbundanceByR21Gamma <- readMat("relAbundanceByR21GammaAsymmetric-5resources.mat")
dfR21Gamma <- reshape2::melt(relAbundanceByR21Gamma)
colnames(dfR21Gamma) <- c("igamma","ir21", "idose", "passage", "species", "relAbundance", "df")

# Import the list of values of gamma, beta, and dose.
gammaList <- seq(0.025,1,0.025)
r21List <- c(1/1.2, 1/1.15, 1/1.1)
doseList <- c("1:1000","1:100","1:10","1:1","10:1","100:1","1000:1")

# Bind the values of r21, gamma, and dose to the matrix.
dfR21Gamma <- dfR21Gamma %>%
  mutate(r21=r21List[ir21], gamma=gammaList[igamma], dose=doseList[idose])

dfR21GammaDoseDependence <- dfR21Gamma %>%
  dplyr::select(-ir21, -igamma, -idose) %>%
  pivot_wider(names_from=dose, values_from=relAbundance) %>%
  mutate(doseDependence=ifelse(species==1, `1000:1`/`1:1000`,`1:1000`/`1000:1`))

p_trajectoriesAsymmetricR21 <- dfR21Gamma %>% 
  filter(gamma=="0.975" & species!=3 & passage %in% c(3,5)) %>% 
  mutate(plot = case_when(r21 == 1/1.2 ~ "R2ab/R1ab = 1.2",
                          r21 == 1/1.15 ~ "R2ab/R1ab = 1.15",
                          r21 == 1/1.1 ~ "R2ab/R1ab = 1.1")) %>% 
  ggplot() +
  geom_line(aes(x=fct_relevel(dose, ratios), y=log10(relAbundance), group=interaction(species, passage),
                linetype=factor(passage), color=factor(species)), alpha=0.6, size=0.4) +
  geom_point(aes(x=fct_relevel(dose, ratios), y=log10(relAbundance), shape=factor(passage), color=factor(species)), 
             alpha=0.8, size=0.5) +
  #geom_line(aes(x=fct_relevel(dose, ratios), y=relAbundance, group=interaction(species, passage),
  #              linetype=factor(passage), color=factor(species)), alpha=0.6, size=0.4) +
  #geom_point(aes(x=fct_relevel(dose, ratios), y=relAbundance, shape=factor(passage), color=factor(species)), 
  #           alpha=0.8, size=0.5) +
  scale_color_manual(name="Species",
                     labels=c("1"="Species 1", "2"="Species 2"), 
                     values=c("1"="#943F92", "2"="#148989")) +
  scale_linetype_manual(labels=c("Passage 3","Passage 5"),
                        values=c("dotted","solid"), guide="none") +
  scale_shape_manual(values=c(1, 19), guide="none") +
  xlab("Mixture ratio") +
  #scale_y_continuous(name="Relative abundance",
  #                   limits=c(-2,0), breaks=c(-2,-1,0),
  #                   labels=label_math(10^.x)) +
  scale_y_continuous(name="Relative abundance",
                     limits=c(-1.5,0), breaks=c(-1,0),
                     labels=label_math(10^.x)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  facet_wrap(~plot) +
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_trajectoriesAsymmetricR21
save_plot("figs2/trajectoriesAsymmetricR21.pdf", p_trajectoriesAsymmetricR21, base_width=2.8, base_height=1.5)
save_plot("figs/trajectoriesAsymmetricR21.png", p_trajectoriesAsymmetricR21, base_width=2.8, base_height=1.5)

# Dose dependence by gamma, across rho (relative consumption of shared resource u13) -----------------------------------

# Import the data that relates dose dependence over passages to gamma.
relAbundanceByR3Gamma <- readMat("relAbundanceByR3GammaAll-5resources.mat")
dfR3Gamma <- reshape2::melt(relAbundanceByR3Gamma)
colnames(dfR3Gamma) <- c("igamma","ir3", "idose", "passage", "species", "relAbundance", "df")

# Import the list of values of gamma, beta, and dose.
gammaList <- seq(0.025,1,0.025)
r3List <- 10^seq(-2,2,0.5)
doseList <- c("1:1000","1000:1")

# Bind the values of r3, gamma, and dose to the matrix.
dfR3Gamma <- dfR3Gamma %>%
  mutate(r3=r3List[ir3], gamma=gammaList[igamma], dose=doseList[idose])

dfR3GammaDoseDependence <- dfR3Gamma %>%
  dplyr::select(-ir3, -igamma, -idose) %>%
  pivot_wider(names_from=dose, values_from=relAbundance) %>%
  mutate(doseDependence=ifelse(species==1, `1000:1`/`1:1000`,`1:1000`/`1000:1`))

p_gammaDDatR3 <- dfR3GammaDoseDependence %>%
  filter(passage==5, species %in% c(1,2,3), r3 %in% 10^seq(-2,2,2)) %>% 
  mutate(plot = case_when(species==1 ~ "Species 1", 
                          species==2 ~ "Species 2",
                          species==3 ~ "Species 3")) %>%
  ggplot() +
  geom_line(aes(x=gamma, y=log10(doseDependence), group=interaction(species, r3),
                color=factor(species), linetype=factor(r3)), size=0.4) +
  scale_linetype_manual(values=c("dashed","solid","dotted", "dotdash"), 
                        name="Relative consumption\nrate of resource ac,\nspecies 3 and 1 (r)") +
  scale_color_manual(name="Species",
                     labels=c("1"="Species 1", "2"="Species 2"), 
                     values=c("1"="#943F92", "2"="#148989"), guide="none") +
  scale_y_continuous(name="Dose dependence", 
                     limits=c(-0.05,3.05), breaks=c(0,1,2,3), 
                     labels=label_math(10^.x)) +
  scale_x_continuous(name="Niche overlap, species 1 and 2 (.)",
                     limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1),
                     labels=c("0","0.25","0.5","0.75","1")) +
  theme(legend.key.size=unit(0.5,"lines")) +
  #theme(legend.position="none") +
  facet_wrap(~plot) +
  DEFAULTS.THEME_PRINT
p_gammaDDatR3
save_plot("figs2/gammaDDAtR3.pdf", p_gammaDDatR3, base_width=4.5, base_height=1.6)
#save_plot("figs/gammaDDAtR3.pdf", p_gammaDDatR3, base_width=3.6, base_height=1.6)

consumptionRateLegend <- get_legend(p_gammaDDatR3)
save_plot("figs/consumptionRateRhoLegend.pdf", consumptionRateLegend, base_width=1, base_height=0.75)

# Dose dependence by gamma, more passages -------------------------------------------------------

# I converted the original model parameter alpha (unique resource / shared resource)
# to gamma (shared resource / (unique resource + shared resource))

# Import the data that relates dose dependence over passages to gamma.
relAbundanceByGamma <- readMat("relAbundanceByGammaMorePassages-5resources.mat")
dfGamma <- reshape2::melt(relAbundanceByGamma)
colnames(dfGamma) <- c("igamma", "idose", "passage", "species", "relAbundance", "df")

# Import the list of values of gamma and dose.
gammaList <- seq(0.025,1,0.025)
#gammaList[40] <- 0.999
doseList <- c("1:1000","1000:1")
# Bind the values of gamma and dose to the matrix.
# Also calculate the original model parameter alpha based on gamma.
# Alpha is the ratio of the unique niche to the shared niche.
dfGamma <- dfGamma %>%
  mutate(gamma=gammaList[igamma], dose=doseList[idose],
         alpha=(1-gamma)/(2*gamma))

dfGammaDoseDependence <- dfGamma %>%
  dplyr::select(-igamma, -idose) %>%
  pivot_wider(names_from=dose, values_from=relAbundance) %>%
  mutate(doseDependence=ifelse(species==1, `1000:1`/`1:1000`,`1:1000`/`1000:1`))

p_gammaDDpassages <- dfGammaDoseDependence %>%
  filter(species==1, gamma %in% gammaList[c(0.75/0.025, 0.9/0.025, 0.975/0.025)]) %>%  
  ggplot() +
  geom_line(aes(x=passage, y=log10(doseDependence), group=gamma, color=factor(gamma))) +
  #scale_y_log10() + 
  xlim(0,75) +
  scale_color_manual(values = c("#77A440","#707FBE","#C33D54"), 
                     name="Niche overlap\n(.)") +
  xlab("# passages") + 
  ylab("Dose dependence") +
  scale_y_continuous(name="Dose dependence", limits=c(0,2.1), breaks = c(0,1,2), labels = label_math(10^.x)) +
  theme(legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_gammaDDpassages
save_plot("figs2/gammaDDPassages.pdf", p_gammaDDpassages, base_width=1.4, base_height=1.2)
save_plot("figs2/gammaDDPassages.png", p_gammaDDpassages, base_width=1.4, base_height=1.2)

gammaLegend <- get_legend(p_gammaDDpassages)
save_plot("figs/gammaLegend.pdf", gammaLegend, base_width=0.6, base_height=0.75)