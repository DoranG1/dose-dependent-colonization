library(tidyverse)
library(cowplot)
library(foreach)
library(RColorBrewer)
library(gcplyr)
library(data.table)
library(scales)

theme_set(theme_cowplot())
source("../../config/palettes/plotDefaults.R")

# Import plate reader data for e0049 passage 2.
e0049p2 <- fread("data/e0049.2-cleaned.txt") %>% 
  select(!`T° 600`) %>% 
  separate(Time, into=c("hour", "minute", "second"), sep=":") %>% 
  mutate(time = as.numeric(as.integer(hour)*60 + as.integer(minute) + as.integer(second)/60)) %>% 
  select(-c(hour, minute, second)) %>% 
  pivot_longer(!time, names_to="well", values_to="OD600") %>% 
  mutate(row = substr(well, 1, 1), col = as.numeric(substr(well, 2, length(well))))

# Import plate reader data for e0043 passage 2.
e0043p2 <- fread("data/e0043.2-cleaned.txt") %>% 
  select(!`T° 600`) %>% 
  separate(Time, into=c("hour", "minute", "second"), sep=":") %>% 
  mutate(time = as.numeric(as.integer(hour)*60 + as.integer(minute) + as.integer(second)/60)) %>% 
  select(-c(hour, minute, second)) %>% 
  pivot_longer(!time, names_to="well", values_to="OD600") %>% 
  mutate(row = substr(well, 1, 1), col = as.numeric(substr(well, 2, length(well))))

# Import plate reader data for e0049 passage 2.
e0043p5 <- fread("data/e0043.5-cleaned.txt") %>% 
  select(!`T° 600`) %>% 
  separate(Time, into=c("hour", "minute", "second"), sep=":") %>% 
  mutate(time = as.numeric(as.integer(hour)*60 + as.integer(minute) + as.integer(second)/60)) %>% 
  select(-c(hour, minute, second)) %>% 
  pivot_longer(!time, names_to="well", values_to="OD600") %>% 
  mutate(row = substr(well, 1, 1), col = as.numeric(substr(well, 2, length(well))))

# Import e0051 plate reader data.
e0051 <- fread("data/e0051.1-cleaned.txt") %>% 
  select(!`T° 600`) %>% 
  separate(Time, into=c("hour", "minute", "second"), sep=":") %>% 
  mutate(time = as.numeric(as.integer(hour)*60 + as.integer(minute) + as.integer(second)/60)) %>% 
  select(-c(hour, minute, second)) %>% 
  pivot_longer(!time, names_to="well", values_to="OD600") %>% 
  mutate(row = substr(well, 1, 1), col = as.numeric(substr(well, 2, length(well)))) %>% 
  mutate(strain = ifelse(row %in% c("A","B","G","H") | col %in% c(1,12), "blank", "other"),
         strain = case_when(strain=="blank" ~ "blank",
                            strain=="other" & (row %in% c("C","E") & col%%4==2) |
                              (row %in% c("D","F") & col%%4==0) ~ "EnteC-2",
                            strain=="other" & (row %in% c("C","E") & col%%4==3) |
                              (row %in% c("D","F") & col%%4==1) ~ "EnteC-3",
                            strain=="other" & (row %in% c("C","E") & col%%4==0) |
                              (row %in% c("D","F") & col%%4==2) ~ "Strep-7",
                            strain=="other" & (row %in% c("C","E") & col%%4==1) |
                              (row %in% c("D","F") & col%%4==3) ~ "Strep-17"))

# Import e0052 plate reader data.
e0052 <- fread("data/e0052.1-cleaned.txt") %>% 
  select(!`T° 600`) %>% 
  separate(Time, into=c("hour", "minute", "second"), sep=":") %>% 
  mutate(time = as.numeric(as.integer(hour)*60 + as.integer(minute) + as.integer(second)/60)) %>% 
  select(-c(hour, minute, second)) %>% 
  pivot_longer(!time, names_to="well", values_to="OD600") %>% 
  mutate(row = substr(well, 1, 1), col = as.numeric(substr(well, 2, length(well)))) %>% 
  mutate(dilution = ifelse(row %in% c("A","B","G","H") | col %in% c(1,12), "blank", "other"),
         dilution = case_when(dilution=="blank" ~ "blank",
                              dilution=="other" & (row %in% c("C","E") & col%%4==2) |
                                (row %in% c("D","F") & col%%4==0) ~ "1",
                              dilution=="other" & (row %in% c("C","E") & col%%4==3) |
                                (row %in% c("D","F") & col%%4==1) ~ "1:10",
                              dilution=="other" & (row %in% c("C","E") & col%%4==0) |
                                (row %in% c("D","F") & col%%4==2) ~ "1:100",
                              dilution=="other" & (row %in% c("C","E") & col%%4==1) |
                                (row %in% c("D","F") & col%%4==3) ~ "1:1000"))

# Import e0054 plate reader data.
e0054 <- fread("data/e0054.1-redo-cleaned.txt") %>% 
  select(!`T° 600`) %>% 
  separate(Time, into=c("hour", "minute", "second"), sep=":") %>% 
  mutate(time = as.numeric(as.integer(hour)*60 + as.integer(minute) + as.integer(second)/60)) %>% 
  select(-c(hour, minute, second)) %>% 
  pivot_longer(!time, names_to="well", values_to="OD600") %>% 
  mutate(row = substr(well, 1, 1), col = as.numeric(substr(well, 2, length(well)))) %>% 
  mutate(strain = case_when(row=="C" & col %in% c(2:11) ~ "EnteC-2",
                            row=="D" & col %in% c(2:11) ~ "EnteC-3",
                            row=="E" & col %in% c(2:11) ~ "Strep-7",
                            row=="F" & col %in% c(2:11) ~ "Strep-17",
                            row %in% c("A","B","G","H") | col %in% c(1,12) ~ "blank"))

# Subtract background from e0054.
e0054 <- e0054 %>% 
  group_by(time) %>% 
  mutate(blankSubtractedOD600 = OD600 - min(OD600))

blankAvg <- mean(e0054 %>% filter(strain=="blank") %>% pull(OD600))
blankMin <- min(e0054 %>% pull(OD600))

# e0054 derivatives.
e0054Deriv <- e0054 %>%
  filter(strain!="blank") %>% 
  group_by(strain, col) %>% 
  mutate(deriv = calc_deriv(x=time, y=log10(OD600)),
         blankSubtractedDeriv = calc_deriv(x=time, y=log10(blankSubtractedOD600)),
         blankSubtractedDerivCalc = calc_deriv(x=time, y=OD600, percapita=TRUE, trans_y="log", blank=blankMin),
         blankSubtractedDerivCalcSmoothed = calc_deriv(x=time, y=OD600, percapita=TRUE, trans_y="log", blank=blankMin))

# Import e0055.A plate reader data (EnteC-3)
e0055A <- fread("data/e0055.A1-redo-cleaned.txt") %>% 
  select(!`T° 600`) %>% 
  separate(Time, into=c("hour", "minute", "second"), sep=":") %>% 
  mutate(time = as.numeric(as.integer(hour)*60 + as.integer(minute) + as.integer(second)/60)) %>% 
  select(-c(hour, minute, second)) %>% 
  pivot_longer(!time, names_to="well", values_to="OD600") %>% 
  mutate(row = substr(well, 1, 1), col = as.numeric(substr(well, 2, length(well)))) %>% 
  mutate(dilution = case_when(row=="C" & col %in% c(2:11) ~ "1",
                              row=="D" & col %in% c(2:11) ~ "1:10",
                              row=="E" & col %in% c(2:11) ~ "1:100",
                              row=="F" & col %in% c(2:11) ~ "1:1000",
                              row %in% c("A","B","G","H") | col %in% c(1,12) ~ "blank"))

# Plots. ------------------------------------------------------------------

# Plot growth curves for entire plate, e0049 p2.
p_e0049p2 <- e0049p2 %>% 
  ggplot() +
  geom_line(aes(x=time, y=OD600)) +
  facet_grid(~row~col) +
  DEFAULTS.THEME_PRINT
p_e0049p2
save_plot("out/e0049p2.png", p_e0049p2, base_width=12, base_height=8)

# Plot growth curves for entire plate, e0043 p2.
p_e0043p2 <- e0043p2 %>% 
  ggplot() +
  geom_line(aes(x=time, y=OD600)) +
  facet_grid(~row~col) +
  DEFAULTS.THEME_PRINT
p_e0043p2
save_plot("out/e0043p2.png", p_e0043p2, base_width=12, base_height=8)

# Plot growth curves for entire plate, e0043 p5.
p_e0043p5 <- e0043p5 %>% 
  ggplot() +
  geom_line(aes(x=time, y=OD600)) +
  facet_grid(~row~col) +
  DEFAULTS.THEME_PRINT
p_e0043p5
save_plot("out/e0043p5.png", p_e0043p5, base_width=12, base_height=8)

# Plot Strep17 and EnteC2 growth curves overlaid.
p_EnteC2Strep17gc <- e0049p2 %>%
  filter(row=="A", col %in% c(1,2,3,10,11,12),
         time <= 1440) %>%
  mutate(strain = ifelse(col %in% c(1,2,3), "EnteC-2", "Strep-17")) %>% 
  ggplot() +
  geom_line(aes(x=time, y=OD600, group=col, color=strain)) +
  scale_color_manual(values=c("EnteC-2"="#8D61A9",
                              "Strep-17"="#BBD585")) +
  theme(legend.position="none") +
  xlab("Time (min)") +
  ylim(0,1.75) +
  DEFAULTS.THEME_PRINT
p_EnteC2Strep17gc
save_plot("out/EnteC2Strep17gc.png", p_EnteC2Strep17gc, base_width=4, base_height=2.5)

# Plot Strep17 and EnteC3 groth curves overlaid.
p_EnteC3Strep17gc <- e0049p2 %>%
  filter(row=="A", col %in% c(5,6,7,10,11,12), 
         time <= 1440) %>%
  mutate(strain = ifelse(col %in% c(5,6,7), "EnteC-3", "Strep-17")) %>% 
  ggplot() +
  geom_line(aes(x=time, y=OD600, group=col, color=strain)) +
  scale_color_manual(values=c("EnteC-3"="#6091AB",
                              "Strep-17"="#BBD585")) +
  theme(legend.position="none") +
  xlab("Time (min)") +
  ylim(0,1.75) +
  DEFAULTS.THEME_PRINT
p_EnteC3Strep17gc
save_plot("out/EnteC3Strep17gc.png", p_EnteC3Strep17gc, base_width=4, base_height=2.5)

# Plot Strep7 and Strep17 growth curves overlaid, passage 2.
p_Strep7Strep17gc2 <- e0043p2 %>% 
  filter(row=="A" & col %in% c(1,2,3,4,5,6),
         time <= 1440) %>% 
  mutate(strain = ifelse(col %in% c(1,2,3), "Strep-7", "Strep-17")) %>% 
  ggplot() +
  geom_line(aes(x=time, y=OD600, group=col, color=strain)) +
  scale_color_manual(values=c("Strep-7"="#6D2E46",
                              "Strep-17"="#BBD585")) +
  theme(legend.position="none") +
  xlab("Time (min)") +
  ylim(0,1.75) +
  DEFAULTS.THEME_PRINT
p_Strep7Strep17gc2
save_plot("out/Strep7Strep17gc2.png", p_Strep7Strep17gc2, base_width=4, base_height=2.5)

# Plot Strep7 and Strep17 growth curves overlaid, passage 5.
p_Strep7Strep17gc5 <- e0043p5 %>% 
  filter(row=="A" & col %in% c(1,2,3,4,5,6),
         time <= 1440) %>% 
  mutate(strain = ifelse(col %in% c(1,2,3), "Strep-7", "Strep-17")) %>% 
  ggplot() +
  geom_line(aes(x=time, y=OD600, group=col, color=strain)) +
  scale_color_manual(values=c("Strep-7"="#6D2E46",
                              "Strep-17"="#BBD585")) +
  theme(legend.position="none") +
  xlab("Time (min)") +
  ylim(0,1.75) +
  DEFAULTS.THEME_PRINT
p_Strep7Strep17gc5
save_plot("out/Strep7Strep17gc5.png", p_Strep7Strep17gc5, base_width=4, base_height=2.5)

# Plot growth curves together.
p_allStrainsgc <- e0049p2 %>%
  filter(row=="A" & col %in% c(1,2,3,5,6,7,10,11,12) & time <= 1440) %>% 
  mutate(strain = case_when(col %in% c(1,2,3) ~ "EnteC-2",
                            col %in% c(5,6,7) ~ "EnteC-3",
                            col %in% c(10,11,12) ~ "Strep-17")) %>% 
  rbind(e0043p5 %>%
          filter(row=="A" & col %in% c(1,2,3,4,5,6) & time <= 1440) %>% 
          mutate(strain = case_when(col %in% c(1,2,3) ~ "Strep-7",
                                    col %in% c(4,5,6) ~ "Strep-17"))) %>% 
  ggplot() + 
  geom_line(aes(x=time, y=OD600, group=interaction(col, strain), color=strain)) +
  scale_color_manual(labels = c("EnteC-2","EnteC-3","Strep-7","Strep-17"),
                     values=c("EnteC-2"="#8D61A9",
                              "EnteC-3"="#6091AB",
                              "Strep-7"="#6D2E46",
                              "Strep-17"="#BBD585")) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  xlab("Time (min)") +
  ylim(0,1.75) +
  DEFAULTS.THEME_PRINT
p_allStrainsgc
#save_plot("out/allStrainsgc.png", p_allStrainsgc, base_width=4, base_height=2.5)
save_plot("out/allStrainsgc.pdf", p_allStrainsgc, base_width=2.5, base_height=1.75)

gcStrainLegend <- get_legend(p_allStrainsgc)
save_plot("out/gcStrainLegend.pdf", gcStrainLegend, base_width=0.5, base_height=0.5)

# Plot Strep-17 growth curves at different dilutions.
p_Strep17Dilutiongc <- e0049p2 %>% 
  filter(time <= 1440 & ((row=="A" & col %in% c(10,11,12)) |
           (row %in% c("B","C","D","E","F","G") & col==9))) %>% 
  mutate(dilution = case_when(row=="A" ~ "1",
                              row %in% c("B","C") ~ "1:10",
                              row %in% c("D","E") ~ "1:100",
                              row %in% c("F","G") ~ "1:1000"),
         rep = case_when(col==10 | row %in% c("B","D","F") ~ 1,
                         col==11 | row %in% c("C","E","G") ~ 2,
                         col==12 ~ 3)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=OD600, group=interaction(dilution, rep), color=dilution)) +
  scale_color_brewer(palette="Dark2") +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  #theme(legend.position="none") +
  xlab("Time (min)") +
  ylim(0,1.75) +
  DEFAULTS.THEME_PRINT
p_Strep17Dilutiongc
save_plot("out/Strep17Dilutiongc.pdf", p_Strep17Dilutiongc, base_width=2.5, base_height=1.75)

gcStrep17DilutionLegend <- get_legend(p_Strep17Dilutiongc)
save_plot("out/gcStrep17DilutionLegend.pdf", gcStrep17DilutionLegend, base_width=0.5, base_height=0.5)

# Plot growth curves for entire plate, e0051.
p_e0051 <- e0051 %>% 
  ggplot() +
  geom_line(aes(x=time, y=OD600)) +
  facet_grid(~row~col) +
  DEFAULTS.THEME_PRINT
p_e0051
save_plot("out/e0051.png", p_e0051, base_width=12, base_height=8)

# Plot 10 growth curves for each strain from e0051.
p_gcAllStrains <- e0051 %>% 
  filter(strain!="blank"& time<=1440) %>% 
  group_by(strain, time) %>% 
  summarize(avgOD = mean(OD600), minOD = mean(OD600)-sd(OD600), maxOD = mean(OD600)+sd(OD600)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=avgOD, color=strain)) +
  geom_ribbon(aes(x=time, ymin=minOD, ymax=maxOD, fill=strain), alpha=0.1) +
  scale_color_manual(labels = c("EnteC-2","EnteC-3","Strep-7","Strep-17"),
                     values=c("EnteC-2"="#8D61A9",
                              "EnteC-3"="#6091AB",
                              "Strep-7"="#6D2E46",
                              "Strep-17"="#BBD585")) +
  scale_fill_manual(labels = c("EnteC-2","EnteC-3","Strep-7","Strep-17"),
                     values=c("EnteC-2"="#8D61A9",
                              "EnteC-3"="#6091AB",
                              "Strep-7"="#6D2E46",
                              "Strep-17"="#BBD585")) +
  theme(legend.position="none") +
  xlab("Time (min)") +
  ylab("OD600") +
  DEFAULTS.THEME_PRINT
p_gcAllStrains
#save_plot("out/gcAllStrains.png", p_gcAllStrains)
save_plot("out/gcAllStrains.pdf", p_gcAllStrains, base_width=2.5, base_height=1.75)

# Plot growth curves for entire plate, e0052.
p_e0052 <- e0052 %>% 
  ggplot() +
  geom_line(aes(x=time, y=OD600)) +
  facet_grid(~row~col) +
  DEFAULTS.THEME_PRINT
p_e0052
save_plot("out/e0052.png", p_e0052, base_width=12, base_height=8)

# Plot 10 growth curves for each strain from e0051.
p_gcAllStrep17Dilutions <- e0052 %>% 
  filter(dilution!="blank"& time<=1440) %>% 
  group_by(dilution, time) %>% 
  summarize(avgOD = mean(OD600), minOD = mean(OD600)-sd(OD600), maxOD = mean(OD600)+sd(OD600)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=avgOD, color=dilution)) +
  geom_ribbon(aes(x=time, ymin=minOD, ymax=maxOD, fill=dilution), alpha=0.1) +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="none") +
  xlab("Time (min)") +
  ylab("OD600") +
  DEFAULTS.THEME_PRINT
p_gcAllStrep17Dilutions
# save_plot("out/gcAllStrep17Dilutions.png", p_gcAllStrep17Dilutions)
save_plot("out/gcAllStrep17Dilutions.pdf", p_gcAllStrep17Dilutions, base_width=2.5, base_height=1.6)

# Plot growth curves for entire plate, e0054.
p_e0054 <- e0054 %>% 
  ggplot() +
  geom_line(aes(x=time, y=OD600)) +
  facet_grid(~row~col) +
  DEFAULTS.THEME_PRINT
p_e0054
save_plot("out/e0054.png", p_e0054, base_width=12, base_height=8)

# Plot 10 growth curves for each strain from e0054.
p_gcAllStrainse0054 <- e0054 %>% 
  filter(strain!="blank"& time<=1440) %>% 
  group_by(strain, time) %>% 
  summarize(avgOD = mean(OD600), minOD = mean(OD600)-sd(OD600), maxOD = mean(OD600)+sd(OD600)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=avgOD, color=strain)) +
  geom_ribbon(aes(x=time, ymin=minOD, ymax=maxOD, fill=strain), alpha=0.1) +
  scale_color_manual(labels = c("EnteC-2","EnteC-3","Strep-7","Strep-17"),
                     values=c("EnteC-2"="#8D61A9",
                              "EnteC-3"="#6091AB",
                              "Strep-7"="#EF4338",
                              "Strep-17"="#FFC43D")) +
  scale_fill_manual(labels = c("EnteC-2","EnteC-3","Strep-7","Strep-17"),
                    values=c("EnteC-2"="#8D61A9",
                             "EnteC-3"="#6091AB",
                             "Strep-7"="#EF4338",
                             "Strep-17"="#FFC43D")) +
  theme(legend.position="none") +
  xlab("Time (min)") +
  ylab("OD600") +
  ylim(0,1.3) +
  DEFAULTS.THEME_PRINT
p_gcAllStrainse0054
save_plot("out/gcAllStrainse0054.png", p_gcAllStrainse0054)
save_plot("out/gcAllStrainse0054.pdf", p_gcAllStrainse0054, base_width=2.5, base_height=1.6)

# Plot blank-subtracted growth curves for all strains from e0054.
p_gcAllStrainse0054blankSubtracted <- e0054 %>% 
  filter(strain!="blank"& time<=600) %>% 
  group_by(strain, time) %>% 
  summarize(avgOD = mean(log10(blankSubtractedOD600)), 
            minOD = mean(log10(blankSubtractedOD600))-sd(log10(blankSubtractedOD600)), 
            maxOD = mean(log10(blankSubtractedOD600))+sd(log10(blankSubtractedOD600))) %>% 
  ggplot() +
  geom_line(aes(x=time, y=avgOD, color=strain)) +
  geom_ribbon(aes(x=time, ymin=minOD, ymax=maxOD, fill=strain), alpha=0.1) +
  scale_color_manual(labels = c("EnteC-2"="E. faecalis",
                                "EnteC-3"="E. casseliflavus",
                                "Strep-7"="L. garvieae",
                                "Strep-17"="L. lactis"),
                     values=c("EnteC-2"="#C1A92F",
                              "EnteC-3"="#BF732E",
                              "Strep-7"="#6398A4",
                              "Strep-17"="#343795")) +
  scale_fill_manual(labels = c("EnteC-2"="E. faecalis",
                               "EnteC-3"="E. casseliflavus",
                               "Strep-7"="L. garvieae",
                               "Strep-17"="L. lactis"),
                    values=c("EnteC-2"="#C1A92F",
                             "EnteC-3"="#BF732E",
                             "Strep-7"="#6398A4",
                             "Strep-17"="#343795")) +
  theme(legend.position="none") +
  xlab("Time (min)") +
  scale_y_continuous(name="OD600",
                     limits=c(-3,0.1), breaks=c(-3,-2,-1,0),
                     labels=label_math(10^.x)) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  guides(color=guide_legend(ncol=2)) +
  #ylim(0,1.3) +
  DEFAULTS.THEME_PRINT
p_gcAllStrainse0054blankSubtracted
save_plot("out/gcAllStrainse0054blankSubtracted.png", p_gcAllStrainse0054blankSubtracted)
save_plot("out/gcAllStrainse0054blankSubtracted.pdf", p_gcAllStrainse0054blankSubtracted, base_width=1.5, base_height=1.3)

gcLegende0054 <- get_legend(p_gcAllStrainse0054blankSubtracted)
save_plot("out/gcLegende0054-2col.pdf", gcLegende0054, base_width=1.5, base_height=0.5)

# Plot derivatives of growth curves for all strains from e0054.
p_allStrainsDerivativese0054 <- e0054Deriv %>% 
  filter(time<=1440) %>% 
  group_by(strain, time) %>% 
  summarize(avgDeriv = mean(deriv), minDeriv=mean(deriv)-sd(deriv), maxDeriv=mean(deriv)+sd(deriv)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=avgDeriv, color=strain)) +
  geom_ribbon(aes(x=time, ymin=minDeriv, ymax=maxDeriv, fill=strain), alpha=0.1) +
  scale_color_manual(labels = c("EnteC-2","EnteC-3","Strep-7","Strep-17"),
                     values=c("EnteC-2"="#8D61A9",
                              "EnteC-3"="#6091AB",
                              "Strep-7"="#EF4338",
                              "Strep-17"="#FFC43D")) +
  scale_fill_manual(labels = c("EnteC-2","EnteC-3","Strep-7","Strep-17"),
                    values=c("EnteC-2"="#8D61A9",
                             "EnteC-3"="#6091AB",
                             "Strep-7"="#EF4338",
                             "Strep-17"="#FFC43D")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT
p_allStrainsDerivativese0054
save_plot("out/allStrainsDerivativese0054.png", p_allStrainsDerivativese0054)

# Plot derivatives of blank-subtracted growth curves for all strains from e0054.
p_allStrainsDerivativese0054blankSubtracted <- e0054Deriv %>% 
  filter(time<=600) %>% 
  group_by(strain, time) %>% 
  summarize(avgDeriv = mean(blankSubtractedDerivCalc)) %>%  
            #minDeriv=mean(blankSubtractedDerivCalc)-sd(blankSubtractedDerivCalc), 
            #maxDeriv=mean(blankSubtractedDerivCalc)+sd(blankSubtractedDerivCalc)) %>% 
  # Smooth curve over 5 min intervals.
  ungroup() %>% 
  group_by(strain) %>% 
  arrange(strain, time) %>% 
  mutate(lag2 = lag(avgDeriv, n=2),
         lag = lag(avgDeriv),
         lead = lead(avgDeriv),
         lead2 = lead(avgDeriv, n=2)) %>% 
  group_by(strain, time) %>% 
  mutate(avgSmoothedDeriv = mean(c(lag2,lag,avgDeriv,lead,lead2), na.rm=TRUE)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=avgSmoothedDeriv, color=strain)) +
  #geom_ribbon(aes(x=time, ymin=minDeriv, ymax=maxDeriv, fill=strain), alpha=0.1) +
  scale_color_manual(labels = c("EnteC-2","EnteC-3","Strep-7","Strep-17"),
                     values=c("EnteC-2"="#C1A92F",
                              "EnteC-3"="#BF732E",
                              "Strep-7"="#6398A4",
                              "Strep-17"="#343795")) +
  scale_fill_manual(labels = c("EnteC-2","EnteC-3","Strep-7","Strep-17"),
                    values=c("EnteC-2"="#C1A92F",
                             "EnteC-3"="#BF732E",
                             "Strep-7"="#6398A4",
                             "Strep-17"="#343795")) +
  theme(legend.position="none") +
  xlab("Time (min)") +
  ylab("Growth rate (OD/min)") +
  ylim(-.009,.030) +
  DEFAULTS.THEME_PRINT
p_allStrainsDerivativese0054blankSubtracted
save_plot("out/allStrainsDerivativese0054blankSubtracted.png", p_allStrainsDerivativese0054blankSubtracted)
save_plot("out/allStrainsDerivativese0054blankSubtracted.pdf", p_allStrainsDerivativese0054blankSubtracted, base_width=1.5, base_height=1.3)
