library(tidyverse)
library(cowplot)
library(foreach)
library(RColorBrewer)
library(gcplyr)
library(data.table)
library(scales)

theme_set(theme_cowplot())
source("../../config/palettes/plotDefaults.R")

# Import e0069 plate reader data.
e0069 <- fread("data/e0069.C1-cleaned.txt") %>% 
  select(!`T° 600`) %>% 
  separate(Time, into=c("hour", "minute", "second"), sep=":") %>% 
  mutate(time = as.numeric(as.integer(hour)*60 + as.integer(minute) + as.integer(second)/60)) %>% 
  select(-c(hour, minute, second)) %>% 
  pivot_longer(!time, names_to="well", values_to="OD600") %>% 
  mutate(row = substr(well, 1, 1), col = as.numeric(substr(well, 2, length(well)))) %>% 
  mutate(strain = case_when(row %in% c("A","H") | col %in% c(1,11,12) ~ "blank",
                            row=="B" & col %in% c(2:10) ~ "L. lactis",
                            row=="C" & col %in% c(2:10) ~ "L. garvieae",
                            row=="D" & col %in% c(2:10) ~ "E. casseliflavus",
                            row=="E" & col %in% c(2:10) ~ "E. faecalis",
                            row=="F" & col %in% c(2:10) ~ "P. goldsteinii",
                            row=="G" & col %in% c(2:10) ~ "B. fragilis"),
         rep = case_when(row %in% c("B","C","D","E","F","G") & col==2 ~ 1,
                         row %in% c("B","C","D","E","F","G") & col==3 ~ 2,
                         row %in% c("B","C","D","E","F","G") & col==4 ~ 3,
                         row %in% c("B","C","D","E","F","G") & col==5 ~ 4,
                         row %in% c("B","C","D","E","F","G") & col==6 ~ 5,
                         row %in% c("B","C","D","E","F","G") & col==7 ~ 6,
                         row %in% c("B","C","D","E","F","G") & col==8 ~ 7,
                         row %in% c("B","C","D","E","F","G") & col==9 ~ 8,
                         row %in% c("B","C","D","E","F","G") & col==10 ~ 9,
                         .default = NA))

# Subtract blanks from e0069 based on first three timepoints of each well.
e0069 <- e0069 %>% 
  group_by(well) %>% 
  mutate(blank = mean(OD600[time %in% unique(time)[1:3]]),
         blankSubtractedOD = OD600 - blank,
         blankSubtractedOD = ifelse(blankSubtractedOD<0.005, 0.005, blankSubtractedOD),
         blankSubtractedDerivCalc = calc_deriv(x=time, y=log10(blankSubtractedOD)))

# Plots. ------------------------------------------------------------------

# Plot blank-subtracted growth curves from e0069.
p_gcAllStrainse0069 <- e0069 %>% 
  filter(strain!="blank") %>% 
  #filter(col %in% c(2,3)) %>% 
  #filter(strain!="blank"& time<=600) %>% 
  group_by(strain, time) %>% 
  summarize(avgOD = mean(log10(blankSubtractedOD)),
           minOD = mean(log10(blankSubtractedOD))-(sd(log10(blankSubtractedOD))/3),
           maxOD = mean(log10(blankSubtractedOD))+(sd(log10(blankSubtractedOD))/3)) %>%
  group_by(strain) %>% 
  arrange(strain, time) %>% 
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
  group_by(strain, time) %>% 
  mutate(avgSmoothedOD = mean(c(lagOD5,lagOD4,lagOD3,lagOD2,lagOD,avgOD,leadOD,leadOD2,leadOD3,leadOD4,leadOD5), na.rm=TRUE),
         avgSmoothedMax = mean(c(lagMax5,lagMax4,lagMax3,lagMax2,lagMax,maxOD,leadMax,leadMax2,leadMax3,leadMax4,leadMax5), na.rm=TRUE),
         avgSmoothedMin = mean(c(lagMin5,lagMin4,lagMin3,lagMin2,lagMin,minOD,leadMin,leadMin2,leadMin3,leadMin4,leadMin5), na.rm=TRUE)) %>%
  ungroup() %>% 
  mutate(order = ifelse(strain %in% c("B. fragilis","P. goldsteinii"), "Bacteroidales", "Lactobacillales")) %>% 
  ggplot() +
  geom_line(aes(x=time, y=avgSmoothedOD, color=strain, linetype=order), alpha=0.9) +
  geom_ribbon(aes(x=time, ymin=avgSmoothedMin, ymax=avgSmoothedMax, fill=strain), alpha=0.1) +
  scale_color_manual(values=c("B. fragilis" = "#762D5D",
                              "P. goldsteinii" = "#B62025", 
                              "E. faecalis"="#C1A92F",
                              "E. casseliflavus"="#BF732E",
                              "L. garvieae"="#6398A4",
                              "L. lactis"="#343795")) +
  scale_fill_manual(values=c("B. fragilis" = "#762D5D",
                             "P. goldsteinii" = "#B62025",
                             "E. faecalis"="#C1A92F",
                             "E. casseliflavus"="#BF732E",
                             "L. garvieae"="#6398A4",
                             "L. lactis"="#343795")) +
  scale_linetype_manual(values = c(2,1)) +
  theme(legend.position="none") +
  xlab("Time (min)") +
  scale_y_continuous(name="OD600",
                     limits=c(-2.5,0.1), breaks=c(-2,-1,0),
                     labels=label_math(10^.x)) +
  xlim(0, 1200) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  guides(color=guide_legend(ncol=2)) +
  DEFAULTS.THEME_PRINT
p_gcAllStrainse0069
save_plot("out/gcAllStrainse0069.png", p_gcAllStrainse0069)
save_plot("out/gcAllStrainse0069.pdf", p_gcAllStrainse0069, base_width=1.7, base_height=1.3)

p_derivse0069 <- e0069 %>% 
  filter(strain!="blank") %>% 
  #filter(col %in% c(2,3)) %>%
  group_by(strain, time) %>% 
  summarize(avgDeriv = mean(blankSubtractedDerivCalc, na.rm=TRUE),
            minDeriv = mean(blankSubtractedDerivCalc, na.rm=TRUE)-(sd(blankSubtractedDerivCalc, na.rm=TRUE)/3),
            maxDeriv = mean(blankSubtractedDerivCalc, na.rm=TRUE)+(sd(blankSubtractedDerivCalc, na.rm=TRUE)/3)) %>% 
  ungroup() %>% 
  group_by(strain) %>% 
  arrange(strain, time) %>% 
  mutate(lag5 = lag(avgDeriv, n=5),
         lag4 = lag(avgDeriv, n=4),
         lag3 = lag(avgDeriv, n=3),
         lag2 = lag(avgDeriv, n=2),
         lag = lag(avgDeriv),
         lead = lead(avgDeriv),
         lead2 = lead(avgDeriv, n=2),
         lead3 = lead(avgDeriv, n=3),
         lead4 = lead(avgDeriv, n=4),
         lead5 = lead(avgDeriv, n=5)) %>% 
  group_by(strain, time) %>% 
  mutate(avgSmoothedDeriv = mean(c(lag5,lag4,lag3,lag2,lag,avgDeriv,lead,lead2,lead3,lead4,lead5), na.rm=TRUE)) %>%
  mutate(order = ifelse(strain %in% c("B. fragilis","P. goldsteinii"), "Bacteroidales", "Lactobacillales")) %>% 
  ggplot() +
  geom_line(aes(x=time, y=avgSmoothedDeriv, group=strain, color=strain, linetype=order), alpha=0.9) +
  #geom_ribbon(aes(x=time, ymin=minDeriv, ymax=maxDeriv, group=strain, fill=strain), alpha=0.1) +
  #geom_ribbon(aes(x=time, ymin=avgDeriv-sdDeriv, ymax=avgDeriv+sdDeriv, group=strain, fill=strain), alpha=0.1) +
  xlim(0, 1200) +
  xlab("Time (min)") +
  ylim(-0.001,0.02) +
  #ylim(-0.005, 0.03) +
  ylab("Growth rate (OD600/min)") +
  scale_color_manual(values=c("B. fragilis" = "#762D5D",
                              "P. goldsteinii" = "#B62025", 
                              "E. faecalis"="#C1A92F",
                              "E. casseliflavus"="#BF732E",
                              "L. garvieae"="#6398A4",
                              "L. lactis"="#343795")) +
  scale_linetype_manual(values = c(2,1)) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_derivse0069
save_plot("out/derivse0069.pdf", p_derivse0069, base_width=1.7, base_height=1.3)
