library(dada2)
library(phyloseq)
library(tidyverse)
library(ggfittext)
library(seqinr)
library(cowplot)
library(exactRankTests)
library(forcats)
library(rlist)
library(foreach)
library(reshape2)

setwd("../../")
source("analysis/generateMixturesDataframe.R")
theme_set(theme_cowplot())
plotdir <- "analysis/betaDiversityPlots"

# Calculate JSD. ----------------------------------------------------------

# Calculate JSD between all samples.
#JSD <- distance(ps_theoretical, "jsd")
#JSD_df <- melt(as.matrix(JSD), varnames=c("sample1","sample2")) %>% rename(JSD=value)

# Parse sample metadata from JSD dataframe.
#JSD_df <- JSD_df %>% 
  #filter(as.character(sample1) >= as.character(sample2)) %>% 
#  separate(sample1, into=c("line1", "label1", "passage1"), sep="_") %>% 
#  separate(sample2, into=c("line2", "label2", "passage2"), sep="_") %>%
#  separate(line1, into=c("mixtureType1", "inoculationReplicate1"), sep="-", remove=FALSE) %>%
#  separate(line2, into=c("mixtureType2", "inoculationReplicate2"), sep="-", remove=FALSE) %>%
#  separate(label1, into=c("subject1a", "ratio1", "subject1b"), sep="-", remove=FALSE) %>% 
#  separate(label2, into=c("subject2a", "ratio2", "subject2b"), sep="-", remove=FALSE) %>% 
#  mutate(combo1=paste(subject1a,subject1b,sep="-"),
#         combo2=paste(subject2a,subject2b,sep="-"),
#         ratio1=fct_relevel(ratio1, ratios),
#         ratio2=fct_relevel(ratio2, ratios),
#         combo1=fct_relevel(combo1, combos),
#         combo2=fct_relevel(combo2, combos))

#write.table(JSD_df, paste0(plotdir, "/JSD.txt"), row.names=FALSE, quote=FALSE, sep="\t")

# Read in JSD dataframe for all samples.
JSD_df <- read.table(paste0(plotdir, "/JSD.txt"), header=TRUE, sep="\t")

# Plot JSD between actual and theoretical mixtures. -----------------------

# Generate list of each unique combination, passage, and ratio.
comboPassageRatio <- expand.grid(comboPassages, ratios) %>% mutate(comboPassageRatio=paste(Var1,Var2,sep="_"))
comboPassageRatio <- comboPassageRatio$comboPassageRatio

# Calculate JSD between different lines within each combo, passage, ratio.
theoreticalActual_JSD_df <- foreach(x=comboPassageRatio, .combine=rbind) %do% {
  combo=gsub("_.*","",x)
  passage=strsplit(x,"_")[[1]][2]
  ratio=gsub(".*_","",x)
  ps_comboPassageRatio <- prune_samples(sample_data(ps_theoretical)$combo==combo &
                                          sample_data(ps_theoretical)$passage==passage & 
                                          sample_data(ps_theoretical)$ratio==ratio,
                                        ps_theoretical)
  comboPassageRatioJSD <- distance(ps_comboPassageRatio, "jsd")
  comboPassageRatioJSD_df <- melt(as.matrix(comboPassageRatioJSD), varnames=c("sample1", "sample2")) %>% rename(JSD=value)
  return(comboPassageRatioJSD_df)
}

# Parse sample metadata from JSD dataframe.
theoreticalActual_JSD_df <- theoreticalActual_JSD_df %>% 
  filter(as.character(sample1) > as.character(sample2)) %>% 
  separate(sample1, into=c("mixtureType1", "label", "passage"), sep="_") %>%
  separate(sample2, into=c("mixtureType2", NA, NA), sep="_") %>% 
  separate(label, into=c("subject1", "ratio", "subject2"), sep="-", remove=FALSE) %>% 
  mutate(mixtureType1=sub("-.*","",mixtureType1),
         mixtureType2=sub("-.*","",mixtureType2),
         combo=paste0(subject1, "-", subject2),
         ratio=fct_relevel(ratio, ratios),
         combo=fct_relevel(combo, combos)) %>% 
  filter(mixtureType1!=mixtureType2)

# Plot theoretical-actual JSD, passage 3.
p_theoreticalActual_JSD_p3 <- theoreticalActual_JSD_df %>% 
  filter(passage==3) %>% 
  ggplot() +
  geom_boxplot(aes(x=ratio, y=JSD, color=combo)) +
  geom_point(aes(x=ratio, y=JSD)) +
  ylim(c(0,0.4)) +
  xlab("ratio") +
  ylab("JSD") +
  ggtitle("Theoretical-actual JSD, p3") +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10)) +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
p_theoreticalActual_JSD_p3
save_plot(paste0(plotdir, "/theoreticalActual_JSD_p3.png"), p_theoreticalActual_JSD_p3, nrow=1.5, ncol=1.75)

# Plot theoretical-actual JSD, passage 5.
p_theoreticalActual_JSD_p5 <- theoreticalActual_JSD_df %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_boxplot(aes(x=ratio, y=JSD, color=combo)) +
  geom_point(aes(x=ratio, y=JSD)) +
  ylim(c(0,0.4)) +
  xlab("ratio") +
  ylab("JSD") +
  ggtitle("Theoretical-actual JSD, p5") +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10)) +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
p_theoreticalActual_JSD_p5
save_plot(paste0(plotdir, "/theoreticalActual_JSD_p5.png"), p_theoreticalActual_JSD_p5, nrow=1.5, ncol=1.75)

# Plot JSD between different mixture ratios. -------------------------------

# Generate list of each unique combination of combo, passage, and line, excluding missing ones.
comboPassageLine <- expand.grid(comboPassages, c("theoretical-1","theoretical-2","actual-1","actual-2","actual-3")) %>% 
  mutate(comboPassageLine=paste(Var1,Var2,sep="_")) %>% 
  filter(!(comboPassageLine %in% c("XBA-XBB_3_theoretical-1",
                                   "XBA-XCA_3_theoretical-1",
                                   "XFA-XBA_3_theoretical-1")))
comboPassageLine <- comboPassageLine$comboPassageLine

# Calculate JSD between different mixture ratios within each combo, passage, line.
mixtureRatio_JSD_df <- foreach(x=comboPassageLine, .combine=rbind) %do% {
  combo=gsub("_.*","",x)
  passage=strsplit(x,"_")[[1]][2]
  line=gsub(".*_","",x)
  ps_comboPassageLine <- prune_samples(sample_data(ps_theoretical)$combo==combo &
                                  sample_data(ps_theoretical)$passage==passage &
                                  sample_data(ps_theoretical)$line==line,
                                ps_theoretical)
  comboPassageLineJSD <- distance(ps_comboPassageLine, "jsd")
  comboPassageLineJSD_df <- melt(as.matrix(comboPassageLineJSD), varnames=c("sample1", "sample2")) %>% rename(JSD=value)
  return(comboPassageLineJSD_df)
}

# Parse sample metadata from JSD dataframe.
mixtureRatio_JSD_df <- mixtureRatio_JSD_df %>% 
  separate(sample1, into=c("line", "label1", "passage"), sep="_") %>%
  separate(sample2, into=c(NA, "label2", NA), sep="_") %>% 
  separate(label1, into=c("subject1", "ratio1", "subject2"), sep="-") %>%
  separate(label2, into=c(NA, "ratio2", NA), sep="-") %>% 
  mutate(combo=paste0(subject1, "-", subject2),
         combo=fct_relevel(combo, combos),
         ratioComparison=paste0(ratio1, "-", ratio2)) %>% 
  filter(ratioComparison %in% c("1:0-1000:1",
                                "1000:1-100:1",
                                "100:1-10:1",
                                "10:1-1:1",
                                "1:1-1:10",
                                "1:10-1:100",
                                "1:100-1:1000",
                                "1:1000-0:1",
                                "1:0-0:1")) %>% 
  mutate(ratioComparison=fct_relevel(ratioComparison, c("1:0-1000:1",
                                                        "1000:1-100:1",
                                                        "100:1-10:1",
                                                        "10:1-1:1",
                                                        "1:1-1:10",
                                                        "1:10-1:100",
                                                        "1:100-1:1000",
                                                        "1:1000-0:1",
                                                        "1:0-0:1")))

# Plot JSD between adjacent mixture ratios at passage 3.
p_mixtureRatio_JSD_p3 <- mixtureRatio_JSD_df %>% 
  filter(passage==3 & ratioComparison!="1:0-0:1") %>% 
  ggplot() +
  geom_line(aes(x=ratioComparison, y=JSD, group=line, color=line)) +
  geom_point(aes(x=ratioComparison, y=JSD)) +
  ylim(c(0,0.32)) +
  xlab("ratio comparison") +
  ylab("JSD") +
  ggtitle("JSD between mixture ratios, p3") +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("theoretical-1"="#313695", 
                              "theoretical-2"="#4575B4", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10)) +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
p_mixtureRatio_JSD_p3
save_plot(paste0(plotdir, "/mixtureRatio_JSD_p3.png"), p_mixtureRatio_JSD_p3, nrow=1.5, ncol=1.75)

# Plot JSD between adjacent mixture ratios at passage 5.
p_mixtureRatio_JSD_p5 <- mixtureRatio_JSD_df %>% 
  filter(passage==5 & ratioComparison!="1:0-0:1") %>% 
  ggplot() +
  geom_line(aes(x=ratioComparison, y=JSD, group=line, color=line)) +
  geom_point(aes(x=ratioComparison, y=JSD)) +
  ylim(c(0,0.32)) +
  xlab("ratio comparison") +
  ylab("JSD") +
  ggtitle("JSD between mixture ratios, p5") +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2",
                              "actual-1"="Actual 1",
                              "actual-2"="Actual 2",
                              "actual-3"="Actual 3"), 
                     values=c("theoretical-1"="#313695", 
                              "theoretical-2"="#4575B4", 
                              "actual-1"="#D73027", 
                              "actual-2"="#F46D43",
                              "actual-3"="#FDAE61")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10)) +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
p_mixtureRatio_JSD_p5
save_plot(paste0(plotdir, "/mixtureRatio_JSD_p5.png"), p_mixtureRatio_JSD_p5, nrow=1.5, ncol=1.75)

# Plot JSD between different parents and passages. ------------------------

p_parents_JSD <- mixtureRatio_JSD_df %>% 
  filter(ratioComparison=="1:0-0:1") %>% 
  ggplot() +
  geom_boxplot(aes(x=passage, y=JSD)) +
  geom_jitter(aes(x=passage, y=JSD, color=line), width=0.1) +
  xlab("passage") +
  ylab("JSD") +
  ggtitle("JSD between parent communities") +
  scale_color_manual(name="Mixture type",
                     labels=c("theoretical-1"="Theoretical 1",
                              "theoretical-2"="Theoretical 2"), 
                     values=c("theoretical-1"="#313695", 
                              "theoretical-2"="#4575B4")) +
  DEFAULTS.THEME_PRES
p_parents_JSD
save_plot(paste0(plotdir, "/parentCommunities_JSD.png"), p_parents_JSD)

# Plot JSD from mixtures to each parent community. ------------------------

# Generate dataframe with JSD to first parent community in one column and JSD to second in another.
dominance_JSD_df <- left_join(
  JSD_df %>% 
    filter(combo1==combo2 & passage1==passage2 & ratio2=="1:0" &
             (line1==line2 | sub("-.*","",line1)=="actual")) %>% 
    select(combo1, passage1, line1, ratio1, ratio2, JSD),
  JSD_df %>% 
    filter(combo1==combo2 & passage1==passage2 & ratio2=="0:1" &
             (line1==line2 | sub("-.*","",line1)=="actual")) %>% 
    select(combo1, passage1, line1, ratio1, ratio2, JSD),
  by=c("combo1","passage1", "line1", "ratio1")) %>% 
  rename(combo=combo1, passage=passage1,
         ratio=ratio1, line=line1,
         parent_1=ratio2.x, parent_2=ratio2.y,
         JSD_1=JSD.x, JSD_2=JSD.y,
         )
dominance_JSD_df <- dominance_JSD_df %>% 
  rbind(dominance_JSD_df %>% 
          filter(ratio %in% c("1:0","0:1")) %>% 
          mutate(line=case_when(line=="theoretical-1" ~ "actual-1",
                                line=="theoretical-2" ~ "actual-2"))) %>% 
  mutate(ratio=fct_relevel(ratio, ratiosFull),
         combo=fct_relevel(combo, combos))

p_dominance_JSD_p3_theoretical <- dominance_JSD_df %>%
  filter(passage==3 & line %in% c("theoretical-1","theoretical-2")) %>% 
  ggplot() +
  geom_point(aes(x=JSD_1, y=JSD_2, color=ratio)) +
  scale_color_brewer(palette="Paired", name="Ratio") +
  xlab("JSD to parent 1:0") +
  ylab("JSD to parent 0:1") +
  xlim(c(0,0.45)) +
  ylim(c(0,0.45)) +
  ggtitle("Theoretical JSD to each parent, p3") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
save_plot(paste0(plotdir, "/dominance_JSD_p3_theoretical.png"), 
          p_dominance_JSD_p3_theoretical, nrow=1.5, ncol=1.75)

p_dominance_JSD_p3_actual <- dominance_JSD_df %>%
  filter(passage==3 & line %in% c("actual-1","actual-2","actual-3")) %>% 
  ggplot() +
  geom_point(aes(x=JSD_1, y=JSD_2, color=ratio)) +
  scale_color_brewer(palette="Paired", name="Ratio") +
  xlab("JSD to parent 1:0") +
  ylab("JSD to parent 0:1") +
  xlim(0,0.45) +
  ylim(0,0.45) +
  ggtitle("JSD to each parent, p3") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
save_plot(paste0(plotdir, "/dominance_JSD_p3_actual.png"), 
          p_dominance_JSD_p3_actual, nrow=1.5, ncol=1.75)

p_dominance_JSD_p3_actual_avg <- dominance_JSD_df %>%
  filter(passage==3 & line %in% c("actual-1","actual-2","actual-3")) %>% 
  group_by(combo,ratio) %>% 
  mutate(JSD_1_avg=mean(JSD_1), JSD_2_avg=mean(JSD_2),
         JSD_1_range=max(JSD_1)-min(JSD_1), JSD_2_range=max(JSD_2)-min(JSD_2)) %>% 
  ggplot() +
  geom_point(aes(x=JSD_1_avg, y=JSD_2_avg, color=ratio)) +
  geom_errorbar(aes(x=JSD_1_avg, color=ratio, ymin=JSD_2_avg-JSD_2_range/2, 
                    ymax=JSD_2_avg+JSD_2_range/2), alpha=0.2) +
  geom_errorbarh(aes(y=JSD_2_avg, color=ratio, xmin=JSD_1_avg-JSD_1_range/2,
                     xmax=JSD_1_avg+JSD_1_range/2), alpha=0.2) +
  scale_color_brewer(palette="Paired", name="Ratio") +
  xlab("JSD to parent 1:0") +
  ylab("JSD to parent 0:1") +
  xlim(0,0.45) +
  ylim(0,0.45) +
  ggtitle("Avg JSD to each parent, p3") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
save_plot(paste0(plotdir, "/dominance_JSD_p3_actual_avg.png"), 
          p_dominance_JSD_p3_actual_avg, nrow=1.5, ncol=1.75)

p_dominance_JSD_p3_all <- dominance_JSD_df %>%
  mutate(lineType=sub("-.*","",line),
         color=ifelse(lineType=="actual",as.character(ratio),"theoretical")) %>% 
  filter(passage==3 & (lineType=="actual" | ratio %in% ratios)) %>% 
  group_by(combo,ratio,lineType) %>% 
  mutate(JSD_1_avg=mean(JSD_1), JSD_2_avg=mean(JSD_2),
         JSD_1_range=max(JSD_1)-min(JSD_1), JSD_2_range=max(JSD_2)-min(JSD_2)) %>% 
  ggplot() +
  geom_point(aes(x=JSD_1_avg, y=JSD_2_avg, color=color)) +
  geom_errorbar(aes(x=JSD_1_avg, color=color, ymin=JSD_2_avg-JSD_2_range/2, 
                    ymax=JSD_2_avg+JSD_2_range/2), alpha=0.2) +
  geom_errorbarh(aes(y=JSD_2_avg, color=color, xmin=JSD_1_avg-JSD_1_range/2,
                     xmax=JSD_1_avg+JSD_1_range/2), alpha=0.2) +
  geom_line(aes(x=JSD_1_avg, y=JSD_2_avg, group=color), linetype="dotted", alpha=0.5) +
  scale_color_manual(values=c("1:0"="#A6CEE3",
                              "1000:1"="#1F78B4",
                              "100:1"="#B2DF8A",
                              "10:1"="#33A02C",
                              "1:1"="#FB9A99",
                              "1:10"="#E31A1C",
                              "1:100"="#FDBF6F",
                              "1:1000"="#FF7F00",
                              "0:1"="#CAB2D6",
                              "theoretical"="#808080"),
                     name="Ratio") +
  xlab("JSD to parent 1:0") +
  ylab("JSD to parent 0:1") +
  xlim(0,0.45) +
  ylim(0,0.45) +
  ggtitle("Avg JSD to each parent, p3") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
save_plot(paste0(plotdir, "/dominance_JSD_p3_all.png"), 
          p_dominance_JSD_p3_all, nrow=1.5, ncol=1.75)

p_dominance_JSD_p5_theoretical <- dominance_JSD_df %>%
  filter(passage==5 & line %in% c("theoretical-1","theoretical-2")) %>% 
  ggplot() +
  geom_point(aes(x=JSD_1, y=JSD_2, color=ratio)) +
  scale_color_brewer(palette="Paired", name="Ratio") +
  xlab("JSD to parent 1:0") +
  ylab("JSD to parent 0:1") +
  xlim(c(0,0.45)) +
  ylim(c(0,0.45)) +
  ggtitle("Theoretical JSD to each parent, p5") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
save_plot(paste0(plotdir, "/dominance_JSD_p5_theoretical.png"), 
          p_dominance_JSD_p5_theoretical, nrow=1.5, ncol=1.75)

p_dominance_JSD_p5_theoretical_avg <- dominance_JSD_df %>%
  filter(passage==5 & line %in% c("theoretical-1","theoretical-2")) %>% 
  group_by(combo,ratio) %>% 
  mutate(JSD_1_avg=mean(JSD_1), JSD_2_avg=mean(JSD_2),
         JSD_1_range=max(JSD_1)-min(JSD_1), JSD_2_range=max(JSD_2)-min(JSD_2)) %>%
  ggplot() +
  geom_point(aes(x=JSD_1_avg, y=JSD_2_avg, color=ratio)) +
  geom_errorbar(aes(x=JSD_1_avg, color=ratio, ymin=JSD_2_avg-JSD_2_range/2, 
                    ymax=JSD_2_avg+JSD_2_range/2), alpha=0.25) +
  geom_errorbarh(aes(y=JSD_2_avg, color=ratio, xmin=JSD_1_avg-JSD_1_range/2,
                     xmax=JSD_1_avg+JSD_1_range/2), alpha=0.25) +
  geom_line(aes(x=JSD_1_avg, y=JSD_2_avg), linetype="dotted", alpha=0.5) +
  scale_color_brewer(palette="Paired", name="Ratio") +
  xlab("JSD to parent 1:0") +
  ylab("JSD to parent 0:1") +
  xlim(c(0,0.45)) +
  ylim(c(0,0.45)) +
  ggtitle("JSD to each parent, p5") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
save_plot(paste0(plotdir, "/dominance_JSD_p5_theoretical_avg.png"), 
          p_dominance_JSD_p5_theoretical_avg, nrow=1.5, ncol=1.75)

p_dominance_JSD_p5_actual <- dominance_JSD_df %>%
  filter(passage==5 & line %in% c("actual-1","actual-2","actual-3")) %>% 
  ggplot() +
  geom_point(aes(x=JSD_1, y=JSD_2, color=ratio)) +
  scale_color_brewer(palette="Paired", name="Ratio") +
  xlab("JSD to parent 1:0") +
  ylab("JSD to parent 0:1") +
  xlim(0,0.45) +
  ylim(0,0.45) +
  ggtitle("JSD to each parent, p5") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
save_plot(paste0(plotdir, "/dominance_JSD_p5_actual.png"), 
          p_dominance_JSD_p5_actual, nrow=1.5, ncol=1.75)

p_dominance_JSD_p5_actual_avg <- dominance_JSD_df %>%
  filter(passage==5 & line %in% c("actual-1","actual-2","actual-3")) %>% 
  group_by(combo,ratio) %>% 
  mutate(JSD_1_avg=mean(JSD_1), JSD_2_avg=mean(JSD_2),
         JSD_1_sd=sd(JSD_1), JSD_2_sd=sd(JSD_2)) %>% 
  ggplot() +
  geom_point(aes(x=JSD_1_avg, y=JSD_2_avg, color=ratio)) +
  geom_errorbar(aes(x=JSD_1_avg, color=ratio, ymin=JSD_2_avg-JSD_2_sd/2, 
                    ymax=JSD_2_avg+JSD_2_sd/2), alpha=0.1) +
  geom_errorbarh(aes(y=JSD_2_avg, color=ratio, xmin=JSD_1_avg-JSD_1_sd/2,
                     xmax=JSD_1_avg+JSD_1_sd/2), alpha=0.1) +
  scale_color_brewer(palette="Paired", name="Ratio") +
  xlab("JSD to parent 1:0") +
  ylab("JSD to parent 0:1") +
  xlim(0,0.45) +
  ylim(0,0.45) +
  ggtitle("Avg JSD to each parent, p5") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
save_plot(paste0(plotdir, "/dominance_JSD_p5_actual_avg.png"), 
          p_dominance_JSD_p5_actual_avg, nrow=1.5, ncol=1.75)

p_dominance_JSD_p5_all <- dominance_JSD_df %>%
  filter(passage==5) %>%
  ungroup() %>% 
  mutate(lineType=sub("-.*","",line),
         color=ifelse(lineType=="actual",as.character(ratio),"theoretical"),
         color=fct_relevel(color, c("1:0","1000:1","100:1","10:1","1:1","1:10","1:100","1:1000","0:1","theoretical"))) %>% 
  mutate(combo = case_when(combo=="XBA-XBB" ~ "A1/A2",
                           combo=="XCA-XCB" ~ "B1/B2",
                           combo=="XDA-XDB" ~ "C1/C2",
                           combo=="XFA-XFB" ~ "D1/D2",
                           combo=="XBA-XCA" ~ "A1/B1",
                           combo=="XCA-XDA" ~ "B1/C1",
                           combo=="XDA-XFA" ~ "C1/D1",
                           combo=="XFA-XBA" ~ "D1/A1"),
         combo = fct_relevel(combo, c("A1/A2","B1/B2","C1/C2","D1/D2","A1/B1","B1/C1","C1/D1","D1/A1"))) %>% 
  group_by(combo,ratio,lineType) %>% 
  summarize(color, 
            JSD_1_avg=mean(JSD_1), JSD_2_avg=mean(JSD_2),
            JSD_1_max=max(JSD_1), JSD_1_min=min(JSD_1), 
            JSD_2_max=max(JSD_2), JSD_2_min=min(JSD_2)) %>% 
  unique() %>% 
  arrange(desc(lineType)) %>%
  #filter(combo %in% c("B1/C1","D1/D2")) %>% 
  #mutate(combo = fct_relevel(combo, "B1/C1","D1/D2")) %>% 
  ggplot() +
  geom_line(aes(x=JSD_1_avg, y=JSD_2_avg, group=color), linetype="dotted", alpha=0.25, size=0.5) +
  geom_point(aes(x=JSD_1_avg, y=JSD_2_avg, color=color), alpha=0.6, size=0.75) +
  geom_errorbar(aes(x=JSD_1_avg, color=color, ymin=JSD_2_min, ymax=JSD_2_max), alpha=0.6, size=0.5) +
  geom_errorbarh(aes(y=JSD_2_avg, color=color, xmin=JSD_1_min, xmax=JSD_1_max), alpha=0.6, size=0.5) +
  # scale_color_manual(values=c("1:0"="#A6CEE3",
  #                             "1000:1"="#1F78B4",
  #                             "100:1"="#B2DF8A",
  #                             "10:1"="#33A02C",
  #                             "1:1"="#FB9A99",
  #                             "1:10"="#E31A1C",
  #                             "1:100"="#FDBF6F",
  #                             "1:1000"="#FF7F00",
  #                             "0:1"="#CAB2D6",
  #                             "theoretical"="#808080"),
  #                    labels=c("1:0"="Parent 1 (1:0)",
  #                             "1000:1"="1000:1",
  #                             "100:1"="100:1",
  #                             "10:1"="10:1",
  #                             "1:1"="1:1",
  #                             "1:10"="1:10",
  #                             "1:100"="1:100",
  #                             "1:1000"="1:1000",
  #                             "0:1"="Parent 2 (0:1)",
  #                             "theoretical"="Theoretical"),
  #                    name="Ratio") +
  scale_color_manual(values=c("1:0"="#260503",
                              "1000:1"="#5F0D07",
                              "100:1"="#99140B",
                              "10:1"="#BE370E",
                              "1:1"="#EF5B2E",
                              "1:10"="#FF7950",
                              "1:100"="#FFA970",
                              "1:1000"="#FFD770",
                              "0:1"="#F8C134",
                              "theoretical"="#808080"),
                     labels=c("1:0"="Parent 1 (1:0)",
                              "1000:1"="1000:1",
                              "100:1"="100:1",
                              "10:1"="10:1",
                              "1:1"="1:1",
                              "1:10"="1:10",
                              "1:100"="1:100",
                              "1:1000"="1:1000",
                              "0:1"="Parent 2 (0:1)",
                              "theoretical"="Theoretical"),
                     name="Ratio") +
  xlab("JSD to parent 1") +
  ylab("JSD to parent 2") +
  xlim(0,0.45) +
  ylim(0,0.45) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) +
  theme(legend.title=element_blank(),
        #axis.text.x=element_text(angle=90),
        legend.key.size=unit(0.5,"lines")) +
  theme(legend.position="none") +
  DEFAULTS.THEME_PRINT +
  facet_wrap(~combo, nrow=2, scales="free") +
  #facet_wrap(~combo, scales="free") +
  theme(strip.text = element_text(face="bold"))
p_dominance_JSD_p5_all
save_plot(paste0(plotdir, "/dominance_JSD_p5_all.pdf"), 
          p_dominance_JSD_p5_all, base_width=6, base_height=3)
#save_plot(paste0(plotdir, "/dominance_JSD_p5_2ex.pdf"),
#                 p_dominance_JSD_p5_2ex, base_width=3, base_height=1.5)


dominanceJSDLegend <- get_legend(p_dominance_JSD_p5_all)
save_plot(paste0(plotdir, "/dominanceJSDLegend.pdf"), dominanceJSDLegend, base_width=0.75, base_height=1.2)

p_dominance_JSD_p5_example_theoretical <- dominance_JSD_df %>%
  filter(passage==5 & line %in% c("theoretical-1","theoretical-2") & combo=="XCA-XDA") %>% 
  group_by(combo,ratio) %>% 
  mutate(JSD_1_avg=mean(JSD_1), JSD_2_avg=mean(JSD_2),
         JSD_1_range=max(JSD_1)-min(JSD_1), JSD_2_range=max(JSD_2)-min(JSD_2)) %>%
  ggplot() +
  geom_point(aes(x=JSD_1_avg, y=JSD_2_avg, color=ratio)) +
  geom_errorbar(aes(x=JSD_1_avg, color=ratio, ymin=JSD_2_avg-JSD_2_range/2, 
                    ymax=JSD_2_avg+JSD_2_range/2), alpha=0.25) +
  geom_errorbarh(aes(y=JSD_2_avg, color=ratio, xmin=JSD_1_avg-JSD_1_range/2,
                     xmax=JSD_1_avg+JSD_1_range/2), alpha=0.25) +
  geom_line(aes(x=JSD_1_avg, y=JSD_2_avg), linetype="dotted", alpha=0.5) +
  scale_color_brewer(palette="Paired", name="Ratio") +
  xlab("JSD to parent 1:0") +
  ylab("JSD to parent 0:1") +
  xlim(c(0,0.45)) +
  ylim(c(0,0.45)) +
  ggtitle("JSD to each parent, p5") +
  DEFAULTS.THEME_PRES +
  facet_wrap(~combo, nrow=2)
p_dominance_JSD_p5_example_theoretical
save_plot(paste0(plotdir, "/dominance_JSD_p5_example_theoretical.png"), 
          p_dominance_JSD_p5_example_theoretical, nrow=1, ncol=.75)

# Plot example with one mixture.
p_dominance_JSD_p5_exampleXDAXFA1ratioActualTheoretical <- dominance_JSD_df %>%
  #filter(passage==5 & combo=="XFA-XFB") %>% 
  #filter(passage==5 & combo=="XCA-XDA") %>%
  #filter(passage==5 & combo=="XFA-XFB") %>% 
  filter(passage==5  & combo=="XDA-XFA") %>% 
  #filter(ratio %in% c("1:0","0:1")) %>% 
  #filter(line %in% c("theoretical-1","theoretical-2") | ratio %in% c("1:0","0:1")) %>% 
  #filter((line %in% c("theoretical-1", "theoretical-2") & ratio=="1:1") |
  #         ratio %in% c("1:0","0:1") & line %in% c("actual-1","actual-2")) %>%
  #filter(ratio=="1:1" | ratio %in% c("1:0","0:1") & line %in% c("actual-1","actual-2")) %>% 
  filter(ratio=="1:1" | (ratio %in% c("1:0","0:1") & line %in% c("actual-1","actual-2")) | 
           line %in% c("theoretical-1","theoretical-2")) %>% 
  mutate(lineType=sub("-.*","",line),
         color=ifelse(lineType=="actual",as.character(ratio),"theoretical")) %>% 
  filter(!(ratio %in% c("1:0","0:1") & lineType=="theoretical")) %>% 
  group_by(combo,ratio,lineType) %>% 
  mutate(JSD_1_avg=mean(JSD_1), JSD_2_avg=mean(JSD_2),
         JSD_1_range=max(JSD_1)-min(JSD_1), JSD_2_range=max(JSD_2)-min(JSD_2)) %>% 
  # Remove errorbars for theoreticals and parents.
  mutate(JSD_1_range = ifelse(line %in% c("theoretical-1","theoretical-2") | ratio %in% c("1:0","0:1"), 0, JSD_1_range),
         JSD_2_range = ifelse(line %in% c("theoretical-1","theoretical-2") | ratio %in% c("1:0","0:1"), 0, JSD_2_range)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(aes(x=JSD_1_avg, y=JSD_2_avg, color=color), size=0.75) +
  geom_errorbar(aes(x=JSD_1_avg, color=color, ymin=JSD_2_avg-JSD_2_range/2, 
                    ymax=JSD_2_avg+JSD_2_range/2), alpha=0.2, size=0.5, width=0, show.legend=FALSE) +
  geom_errorbarh(aes(y=JSD_2_avg, color=color, xmin=JSD_1_avg-JSD_1_range/2,
                     xmax=JSD_1_avg+JSD_1_range/2), alpha=0.2, size=0.5, height=0, show.legend=FALSE) +
  geom_line(aes(x=JSD_1_avg, y=JSD_2_avg, group=color), linetype="dotted", alpha=0.5, size=0.5) +
  scale_color_manual(values=c("1:0"="#A6CEE3",
                              "1000:1"="#1F78B4",
                              "100:1"="#B2DF8A",
                              "10:1"="#33A02C",
                              "1:1"="#FB9A99",
                              "1:10"="#E31A1C",
                              "1:100"="#FDBF6F",
                              "1:1000"="#FF7F00",
                              "0:1"="#CAB2D6",
                              "theoretical"="#808080"),
                     labels=c("1:0"="Parent 1",
                              "1000:1"="1000:1",
                              "100:1"="100:1",
                              "10:1"="10:1",
                              "1:1"="1:1",
                              "1:10"="1:10",
                              "1:100"="1:100",
                              "1:1000"="1:1000",
                              "0:1"="Parent 2",
                              "theoretical"="Theoretical"),
                     name="Ratio") +
  xlab("JSD to parent 1") +
  ylab("JSD to parent 2") +
  xlim(0,0.45) +
  ylim(0,0.45) +
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.5,"lines")) +
  #theme(legend.position="none") +
  #ggtitle("JSD to each parent") +
  #DEFAULTS.THEME_PRES +
  DEFAULTS.THEME_PRINT
  #facet_wrap(~combo, nrow=2)
p_dominance_JSD_p5_exampleXDAXFA1ratioActualTheoretical
save_plot(paste0(plotdir, "/dominance_JSD_p5_exampleXDAXFA1ratioActualTheoretical.png"), 
          p_dominance_JSD_p5_exampleXDAXFA1ratioActualTheoretical, base_width=2.75, base_height=2)
#save_plot(paste0(plotdir, "/dominance_JSD_p5_exampleXCAXDA.png"), 
#          p_dominance_JSD_p5_exampleXCAXDA, base_width=2.025, base_height=1.5)

# Summary JSD. ------------------------------------------------------------

# Summary JSD of all parents, replicates, and theory vs. actual.
summary_JSD_df <- rbind(
  # Different parents from the same passage and replicate.
  JSD_df %>%
    filter(ratio1 %in% c("1:0","0:1") & ratio2 %in% c("1:0","0:1") &
             passage1==passage2 & inoculationReplicate1==inoculationReplicate2 & JSD!=0) %>% 
    mutate(combo1=case_when(ratio1=="1:0"~subject1a, ratio1=="0:1"~subject1b),
           combo2=case_when(ratio2=="1:0"~subject2a, ratio2=="0:1"~subject2b),
           comparison="parents", type="parents") %>% 
    group_by(JSD) %>% 
    slice(1),
  # Replicates from all ratios of the actual mixtures.
  JSD_df %>% 
    filter(mixtureType1==mixtureType2 & ratio1==ratio2 & passage1==passage2 & combo1==combo2 &
             (mixtureType1=="actual" | ratio1 %in% c("1:0","0:1")) & JSD!=0) %>% 
    mutate(comparison="replicates", type="replicates") %>% 
    group_by(JSD) %>% 
    slice(1),
  # Theoretical vs actual for all mixture ratios.
  JSD_df %>% 
    filter(mixtureType1!=mixtureType2 & ratio1==ratio2 & passage1==passage2 & combo1==combo2 & 
             JSD!=0) %>%
    mutate(comparison=combo1, type="theoretical-actual") %>% 
    group_by(JSD) %>% 
    slice(1),
  # Passage 3 vs 5 for all mixture ratios.
  JSD_df %>% 
    filter(combo1==combo2 & ratio1==ratio2 & line1==line2 & passage1!=passage2) %>% 
    mutate(comparison="passage3-passage5", type="passage3-passage5") %>% 
    group_by(JSD) %>% 
    slice(1),
  # Mixtures to parents for all actual mixture ratios.
  JSD_df %>% 
    filter(combo1==combo2 & passage1==passage2 & ratio2=="1:0" & ratio1!="0:1" &
             sub("-.*","",line1)=="actual" & JSD!=0) %>% 
    mutate(comparison=paste0(ratio1,"-",ratio2), type="toParent"),
  JSD_df %>% 
    filter(combo1==combo2 & passage1==passage2 & ratio2=="0:1" & ratio1!="1:0" &
             sub("-.*","",line1)=="actual" & JSD!=0) %>% 
    mutate(comparison=paste0(ratio1,"-",ratio2), type="toParent")
  )

# Plot distribution of JSD between parents, theoretical/actual, and replicates.
p_summary_JSD <- summary_JSD_df %>% 
  filter(passage1==5) %>% 
  filter(type %in% c("parents","theoretical-actual","replicates")) %>%
  mutate(type = case_when(type=="parents" ~ "parents",
                          type=="theoretical-actual" ~ "theoretical\nvs actual",
                          type=="replicates" ~ "replicates")) %>% 
  ggplot(aes(x=factor(type, levels=c("parents","theoretical\nvs actual","replicates")), y=JSD)) +
  #geom_point(aes(color=type), position=position_jitter(width=0.15, seed=1), size=0.25, alpha=0.5) +
  geom_violin(aes(fill=type), size=0.25) +
  geom_boxplot(size=0.25, width=0.1, outlier.shape=NA) +
  scale_fill_manual(values=c("parents"="#4DAF4A",
                              "theoretical\nvs actual"="#377EB8",
                              "replicates"="#E41A1C")) +
  theme(legend.position="none") +
  xlab("Comparison") +
  ylab("JSD") +
  #ggtitle("JSD distributions") +
  ylim(c(0,0.45)) +
  DEFAULTS.THEME_PRINT
p_summary_JSD
save_plot(paste0(plotdir, "/summary_JSD.pdf"), p_summary_JSD, base_width=2, base_height=1.8)
save_plot(paste0(plotdir, "/summary_JSD.png"), p_summary_JSD, nrow=1.5,ncol=1.25)

# Test whether differences in distribution are significant.
aovdata <- summary_JSD_df %>% filter(passage1==5 &
                                       type %in% c("parents","theoretical-actual","replicates")) %>% select(type, JSD) %>% 
  mutate(JSD = as.numeric(JSD))
summaryAOV <- aov(JSD ~ type, data=aovdata)
summary(summaryAOV)
summaryTukey <- TukeyHSD(summaryAOV, conf.level=0.95)

t.test(aovdata %>% filter(type=="parents") %>% pull(JSD),
       aovdata %>% filter(type=="theoretical-actual") %>% pull(JSD))
t.test(aovdata %>% filter(type=="theoretical-actual") %>% pull(JSD),
       aovdata %>% filter(type=="replicates") %>% pull(JSD))
t.test(aovdata %>% filter(type=="parents") %>% pull(JSD),
       aovdata %>% filter(type=="replicates") %>% pull(JSD))

# Wilcoxon test.
JSDparentTheoreticalDiff <- wilcox.test(
  summary_JSD_df %>% filter(passage1==5 & type=="parents") %>% pull(JSD),
  summary_JSD_df %>% filter(passage1==5 & type=="theoretical-actual") %>% pull(JSD))

JSDTheoreticalReplicatesDiff <- wilcox.test(
  summary_JSD_df %>% filter(passage1==5 & type=="theoretical-actual") %>% pull(JSD),
  summary_JSD_df %>% filter(passage1==5 & type=="replicates") %>% pull(JSD))

JSDparentReplicatesDiff <- wilcox.test(
  summary_JSD_df %>% filter(passage1==5 & type=="parents") %>% pull(JSD),
  summary_JSD_df %>% filter(passage1==5 & type=="replicates") %>% pull(JSD))

# Plot distribution of JSD between mixture ratios and both parents for each combo separately.
p_summary_toParentsAllRatios <- summary_JSD_df %>% 
  filter(passage1==5 & type=="toParent") %>% 
  ggplot() +
  #geom_point(aes(x=ratio2, y=JSD), position=position_jitter(width=0.15, seed=1), size=0.25, alpha=0.5) +
  geom_violin(aes(x=ratio2, y=JSD), size=0.25) +
  geom_boxplot(aes(x=ratio2, y=JSD), size=0.25, width=0.1, outlier.shape=NA) +
  xlab("Comparison") +
  ylab("JSD") +
  facet_wrap(~combo1, nrow=2) +
  DEFAULTS.THEME_PRINT
p_summary_toParentsAllRatios

# Plot scatterplot distribution of JSD to both parents at 1:1 mixture ratio.
p_summary_toParentsScatter <- summary_JSD_df %>% 
  filter(passage1==5 & type=="toParent" & ratio1=="1:1" & ratio2=="1:0") %>% 
  select(combo1, line1, line2, JSD) %>% 
  rename(JSD1=JSD) %>% 
  left_join(summary_JSD_df %>% 
              filter(passage1==5 & type=="toParent" & ratio1=="1:1" & ratio2=="0:1") %>% 
              select(combo1, line1, line2, JSD) %>% 
              rename(JSD2=JSD), 
            by = c("combo1","line1","line2")) %>% 
  group_by(combo1) %>% 
  mutate(avgJSD1=mean(JSD1), avgJSD2=mean(JSD2)) %>% 
  ggplot() +
  geom_point(aes(x=avgJSD1, y=avgJSD2), size=0.25) +
  geom_abline(aes(slope=1, intercept=0), linetype="dashed", size=0.25) +
  xlab("JSD to parent 1") +
  ylab("JSD to parent 2") +
  xlim(0,0.4) +
  ylim(0,0.4) +
  DEFAULTS.THEME_PRINT
p_summary_toParentsScatter
save_plot(paste0(plotdir, "/summary_toParentsScatter.pdf"), p_summary_toParentsScatter, base_width=2, base_height=2)

p_summary_1ratio <- summary_JSD_df %>% 
  filter(passage1==5 & comparison %in% c("parents", "replicates", combos, "1:1-1:0", "1:1-0:1")) %>% 
  filter(ratio1=="1:1" | comparison=="parents") %>% 
  group_by(label1, label2, comparison) %>% 
  summarize(combo1, combo2, type, avgJSD=mean(JSD)) %>% 
  unique() %>% 
  group_by(type, combo1) %>% 
  mutate(maxJSD=max(avgJSD),
         comparison=replace(comparison, type=="toParent" & avgJSD==maxJSD, "1:1 vs\nparent 1"),
         comparison=replace(comparison, type=="toParent" & avgJSD<maxJSD, "1:1 vs\nparent 2"),
         comparison=replace(comparison, comparison %in% combos, "theoretical\nvs actual"),
         grouping=ifelse(type=="toParent", combo1, paste(label1, label2, type))) %>% 
  ggplot(aes(x=factor(comparison, levels=c(
    "parents", "theoretical\nvs actual", "1:1 vs\nparent 1", "1:1 vs\nparent 2", "replicates")),
    y=avgJSD)) +
  geom_point(aes(color=comparison), position=position_jitter(width=0.15, seed=1), size=0.25) +
  geom_line(aes(group=grouping, color=type), size=0.15) +
  geom_boxplot(alpha=0, size=0.25) +
  scale_color_manual(values=c("parents"="#4DAF4A",
                              "theoretical\nvs actual"="#377EB8",
                              "replicates"="#E41A1C",
                              "1:1 vs\nparent 1"="#FF7F00",
                              "1:1 vs\nparent 2"="#FF7F00",
                              "toParent"="#FF7F00")) +
  theme(legend.position="none") +
  xlab("Comparison") +
  ylab("JSD") +
  ylim(c(0,0.4)) +
  #ggtitle("JSD between communities")
  DEFAULTS.THEME_PRINT
p_summary_1ratio

p_summary_JSD_p3 <- summary_JSD_df %>% 
  filter(type %in% c("parents","theoretical-actual","replicates"), passage1==3) %>%
  ggplot(aes(x=factor(type, levels=c("parents","theoretical-actual","replicates")), y=JSD)) +
  geom_boxplot() +
  geom_jitter(aes(color=type), width=0.15, alpha=0.25) +
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A")) +
  theme(legend.position="none") +
  xlab("Comparison") +
  ylab("JSD") +
  ggtitle("JSD distributions, p3") +
  ylim(c(0,0.45))
p_summary_JSD_p3
save_plot(paste0(plotdir, "/summary_JSD_p3.png"), p_summary_JSD_p3, nrow=1.5,ncol=1.25)

p_summary_JSD_p5 <- summary_JSD_df %>% 
  filter(type %in% c("parents","theoretical-actual","replicates"), passage1==5) %>%
  ggplot(aes(x=factor(type, levels=c("parents","theoretical-actual","replicates")), y=JSD)) +
  geom_boxplot() +
  geom_jitter(aes(color=type), width=0.15, alpha=0.25) +
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A")) +
  theme(legend.position="none") +
  xlab("Comparison") +
  ylab("JSD") +
  ggtitle("JSD distributions, p5") +
  ylim(c(0,0.45))
p_summary_JSD_p5
save_plot(paste0(plotdir, "/summary_JSD_p5.png"), p_summary_JSD_p5, nrow=1.5,ncol=1.25)

p_summary_JSD_byCombo <- summary_JSD_df %>% 
  filter(type %in% c("parents","theoretical-actual","replicates")) %>% 
  ggplot(aes(x=factor(comparison, 
                      levels=c("parents","XBA-XBB","XCA-XCB","XDA-XDB","XFA-XFB",
                               "XBA-XCA","XCA-XDA","XDA-XFA","XFA-XBA","replicates")), y=JSD)) +
  geom_boxplot() +
  geom_jitter(aes(color=type), width=0.15, alpha=0.25) +
  scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A")) +
  theme(legend.position="none") +
  xlab("Comparison") +
  ylab("JSD") +
  ggtitle("JSD distributions") +
  ylim(c(0,0.45))
p_summary_JSD_byCombo
save_plot(paste0(plotdir, "/summary_JSD_byCombo.png"), p_summary_JSD_byCombo, nrow=1.5,ncol=1.75)

p_summary_JSD_byCombo_p3 <- summary_JSD_df %>% 
  filter(type %in% c("parents","theoretical-actual","replicates"), passage1==3) %>% 
  ggplot(aes(x=factor(comparison, 
                      levels=c("parents","XBA-XBB","XCA-XCB","XDA-XDB","XFA-XFB",
                               "XBA-XCA","XCA-XDA","XDA-XFA","XFA-XBA","replicates")), y=JSD)) +
  geom_boxplot() +
  geom_jitter(aes(color=type), width=0.15, alpha=0.25) +
  scale_color_manual(values=c("parents"="#4DAF4A",
                              "theoretical-actual"="#377EB8",
                              "replicates"="#E41A1C")) +
  theme(legend.position="none") +
  xlab("Comparison") +
  ylab("JSD") +
  ggtitle("JSD distributions, p3") +
  ylim(c(0,0.45))
p_summary_JSD_byCombo_p3
save_plot(paste0(plotdir, "/summary_JSD_byCombo_p3.png"), p_summary_JSD_byCombo_p3, nrow=1.5,ncol=1.75)

p_summary_JSD_byCombo_p5 <- summary_JSD_df %>% 
  filter(type %in% c("parents","theoretical-actual","replicates"), passage1==5) %>% 
  ggplot(aes(x=factor(comparison, 
                      levels=c("parents","XBA-XBB","XCA-XCB","XDA-XDB","XFA-XFB",
                               "XBA-XCA","XCA-XDA","XDA-XFA","XFA-XBA","replicates")), y=JSD)) +
  geom_boxplot() +
  geom_jitter(aes(color=type), width=0.15, alpha=0.25) +
  scale_color_manual(values=c("parents"="#4DAF4A",
                              "theoretical-actual"="#377EB8",
                              "replicates"="#E41A1C")) +
  theme(legend.position="none") +
  xlab("Comparison") +
  ylab("JSD") +
  ggtitle("JSD distributions, p5") +
  ylim(c(0,0.45))
p_summary_JSD_byCombo_p5
save_plot(paste0(plotdir, "/summary_JSD_byCombo_p5.png"), p_summary_JSD_byCombo_p5, nrow=1.5,ncol=1.75)

p_summary_JSD_passages <- summary_JSD_df %>% 
  filter(type %in% c("parents","theoretical-actual","replicates","passage3-passage5")) %>% 
  ggplot(aes(x=factor(type, levels=c("parents","theoretical-actual","replicates","passage3-passage5")), y=JSD)) +
  geom_boxplot() +
  geom_jitter(aes(color=type), width=0.15, alpha=0.25) +
  scale_color_manual(values=c("#E41A1C","#377EB8","#984EA3","#4DAF4A")) +
  theme(legend.position="none") +
  xlab("Comparison") +
  ylab("JSD") +
  ylim(c(0,0.45))
p_summary_JSD_passages
save_plot(paste0(plotdir, "/summary_JSD_passages.png"), p_summary_JSD_passages, nrow=1.5,ncol=1.25)

p_summary_JSD_toParent <- summary_JSD_df %>% 
  ggplot(aes(x=factor(type, 
                      levels=c("parents","theoretical-actual","replicates","passage3-passage5","toParent")),
             y=JSD)) +
  geom_boxplot() +
  geom_jitter(aes(color=type), width=0.15, alpha=0.25) +
  scale_color_manual(values=c("#E41A1C","#FF7F00","#377EB8","#984EA3","#4DAF4A")) +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10)) +
  xlab("Comparison") +
  ylab("JSD") +
  ylim(c(0,0.45))
p_summary_JSD_toParent
save_plot(paste0(plotdir, "/summary_JSD_toParent.png"), p_summary_JSD_toParent, nrow=1.5,ncol=2)

p_summary_JSD_all_p3 <- summary_JSD_df %>% 
  filter(passage1==3) %>% 
  mutate(comparison=replace(comparison,type=="theoretical-actual","theoretical-actual")) %>% 
  ggplot(aes(x=factor(comparison, 
                      levels=c("parents","theoretical-actual","replicates","passage3-passage5",
                               "1000:1-1:0","100:1-1:0","10:1-1:0","1:1-1:0",
                               "1:10-1:0","1:100-1:0","1:1000-1:0",
                               "1000:1-0:1","100:1-0:1","10:1-0:1","1:1-0:1",
                               "1:10-0:1","1:100-0:1","1:1000-0:1")), y=JSD)) +
  geom_boxplot() +
  geom_jitter(aes(color=type), width=0.15, alpha=0.25) +
  scale_color_manual(values=c("parents"="#4DAF4A",
                              "theoretical-actual"="#377EB8",
                              "replicates"="#E41A1C",
                              "passage3-passage5"="#984EA3",
                              "toParent"="#FF7F00")) +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10)) +
  xlab("Comparison") +
  ylab("JSD") +
  ylim(c(0,0.45)) +
  ggtitle("JSD distributions, p3")
p_summary_JSD_all_p3
save_plot(paste0(plotdir, "/summary_JSD_all_p3.png"), p_summary_JSD_all_p3, nrow=1.5,ncol=2)

p_summary_JSD_all_p5 <- summary_JSD_df %>% 
  filter(passage1==5) %>% 
  mutate(comparison=replace(comparison,type=="theoretical-actual","theoretical-actual")) %>% 
  ggplot(aes(x=factor(comparison, 
                      levels=c("parents","theoretical-actual","replicates","passage3-passage5",
                               "1000:1-1:0","100:1-1:0","10:1-1:0","1:1-1:0",
                               "1:10-1:0","1:100-1:0","1:1000-1:0",
                               "1000:1-0:1","100:1-0:1","10:1-0:1","1:1-0:1",
                               "1:10-0:1","1:100-0:1","1:1000-0:1")), y=JSD)) +
  geom_boxplot() +
  geom_jitter(aes(color=type), width=0.15, alpha=0.25) +
  scale_color_manual(values=c("parents"="#4DAF4A",
                              "theoretical-actual"="#377EB8",
                              "replicates"="#E41A1C",
                              "passage3-passage5"="#984EA3",
                              "toParent"="#FF7F00")) +
  theme(legend.position="none") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10)) +
  xlab("Comparison") +
  ylab("JSD") +
  ylim(c(0,0.45)) +
  ggtitle("JSD distributions, p5")
p_summary_JSD_all_p5
save_plot(paste0(plotdir, "/summary_JSD_all_p5.png"), p_summary_JSD_all_p5, nrow=1.5,ncol=2)

# Plot summary boxplots with toParent separation for all combos and passages separately.
foreach(x=comboPassages) %do% {
  combo <- sub("_.*","",x)
  passage <- sub(".*_","",x)
  p_summary_JSD_combo <- summary_JSD_df %>% 
    filter(combo1==combo | (combo1==sub(".*-","",combo) & combo2==sub("-.*","",combo))) %>%
    filter(passage1==passage) %>% 
    mutate(comparison=replace(comparison,type=="theoretical-actual","theoretical-actual")) %>%
    ggplot(aes(x=factor(comparison, 
                        levels=c("parents","theoretical-actual","replicates","passage3-passage5",
                                 "1000:1-1:0","100:1-1:0","10:1-1:0","1:1-1:0",
                                 "1:10-1:0","1:100-1:0","1:1000-1:0",
                                 "1000:1-0:1","100:1-0:1","10:1-0:1","1:1-0:1",
                                 "1:10-0:1","1:100-0:1","1:1000-0:1")), y=JSD)) +
    geom_boxplot() +
    geom_jitter(aes(color=type), width=0.15, alpha=0.25) +
    scale_color_manual(values=c("parents"="#4DAF4A",
                                "theoretical-actual"="#377EB8",
                                "replicates"="#E41A1C",
                                "passage3-passage5"="#984EA3",
                                "toParent"="#FF7F00")) +
    theme(legend.position="none") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=10)) +
    xlab("Comparison") +
    ylab("JSD") +
    ggtitle(paste0("JSD distributions, ",combo,", p",passage)) +
    ylim(c(0,0.45))
  save_plot(paste0(plotdir, "/summary_JSD_", combo, "_p", passage, ".png"), p_summary_JSD_combo, nrow=1.5,ncol=2)
}

# Plot summary boxplots for parents, theoretical-actual (1:1), distance to parents (1:1), and replicates.
p_summary_1ratio <- summary_JSD_df %>% 
  filter(passage1==5 & comparison %in% c("parents", "replicates", combos, "1:1-1:0", "1:1-0:1")) %>% 
  filter(ratio1=="1:1" | comparison=="parents") %>% 
  group_by(label1, label2, comparison) %>% 
  summarize(combo1, combo2, type, avgJSD=mean(JSD)) %>% 
  unique() %>% 
  group_by(type, combo1) %>% 
  mutate(maxJSD=max(avgJSD),
         comparison=replace(comparison, type=="toParent" & avgJSD==maxJSD, "1:1 vs\nparent 1"),
         comparison=replace(comparison, type=="toParent" & avgJSD<maxJSD, "1:1 vs\nparent 2"),
         comparison=replace(comparison, comparison %in% combos, "theoretical\nvs actual"),
         grouping=ifelse(type=="toParent", combo1, paste(label1, label2, type))) %>% 
  ggplot(aes(x=factor(comparison, levels=c(
    "parents", "theoretical\nvs actual", "1:1 vs\nparent 1", "1:1 vs\nparent 2", "replicates")),
             y=avgJSD)) +
  geom_point(aes(color=comparison), position=position_jitter(width=0.15, seed=1), size=0.25) +
  geom_line(aes(group=grouping, color=type), size=0.15) +
  geom_boxplot(alpha=0, size=0.25) +
  scale_color_manual(values=c("parents"="#4DAF4A",
                              "theoretical\nvs actual"="#377EB8",
                              "replicates"="#E41A1C",
                              "1:1 vs\nparent 1"="#FF7F00",
                              "1:1 vs\nparent 2"="#FF7F00",
                              "toParent"="#FF7F00")) +
  theme(legend.position="none") +
  xlab("Comparison") +
  ylab("JSD") +
  ylim(c(0,0.4)) +
  #ggtitle("JSD between communities")
  DEFAULTS.THEME_PRINT
  #DEFAULTS.THEME_PRES
p_summary_1ratio
save_plot(paste0(plotdir, "/summary_1ratio.pdf"), p_summary_1ratio, base_width=2.5, base_height=1.8)
#save_plot(paste0(plotdir, "/summary_1ratio.png"), p_summary_1ratio, nrow=1, ncol=1.2)

# Create table for supplementary figure.
summary_allRatios_JSD_df <- summary_JSD_df %>% 
  filter(passage1==5 & comparison %in% c("parents", "replicates", combos)) %>% 
  group_by(label1, label2, comparison) %>% 
  summarize(combo1, combo2, type, avgJSD=mean(JSD)) %>% 
  unique() %>% 
  filter(type!="parents" | paste(combo1, combo2, sep="-") %in% combos | paste(combo2, combo1, sep="-") %in% combos) %>% 
  group_by(combo1, combo2, type) %>% 
  summarize(avgJSD=mean(avgJSD)) %>% 
  unique() %>% 
  ungroup() %>% 
  mutate(combo=case_when(combo1 %in% combos ~ combo1,
                         paste(combo1, combo2, sep="-") %in% combos ~ paste(combo1, combo2, sep="-"),
                         paste(combo2, combo1, sep="-") %in% combos ~ paste(combo2, combo1, sep="-"))) %>% 
  select(-c(combo1, combo2)) %>% 
  pivot_wider(names_from=type, values_from=avgJSD) %>% 
  left_join(
    toParentJSDsSummary <- left_join(
      summary_JSD_df %>% 
        filter(passage1==5 & type=="toParent" & ratio2=="1:0") %>% 
        group_by(combo1, ratio1) %>%
        summarize(ratio2, avgJSDto1=mean(JSD)) %>% 
        rename(combo=combo1, parent1=ratio2) %>% 
        unique(),
      summary_JSD_df %>% 
        filter(passage1==5 & type=="toParent" & ratio2=="0:1") %>% 
        group_by(combo1, ratio1) %>% 
        summarize(ratio2, avgJSDto2=mean(JSD)) %>%
        rename(combo=combo1, parent2=ratio2) %>% 
        unique(),
      by=c("combo","ratio1")
    ) %>% 
      group_by(combo) %>% 
      summarize(avgJSDto1=mean(avgJSDto1),
                avgJSDto2=mean(avgJSDto2)) %>% 
      unique() %>% 
      ungroup(),
    by=c("combo")
  )

# Analyze correlation of JSD to each parent with proportion of OC/UC ASVs. --------

# Get 1:1 JSD to each parent of actual communities from dominance JSD df.
dominanceToParents <- dominance_JSD_df %>% 
  filter(ratio=="1:1" & sub("-.*","",line)=="actual") %>% 
  group_by(combo, passage) %>% 
  summarize(parent_1, parent_2, 
            avgJSD_1 = mean(JSD_1), minJSD_1 = min(JSD_1), maxJSD_1 = max(JSD_1),
            avgJSD_2 = mean(JSD_2), minJSD_2 = min(JSD_2), maxJSD_2 = max(JSD_2),
            avgJSDratio = mean(JSD_1/JSD_2), minJSDratio = min(JSD_1/JSD_2), maxJSDratio = max(JSD_1/JSD_2)) %>% 
  unique()

OCUCRelAbAll <- theoretical_mixtures_df_categorized %>% 
  filter(category %in% c("overcolonizing","undercolonizing") & ratio %in% c("1:0","0:1")) %>% 
  ungroup() %>% 
  select(combo, passage, line, category, ratio, rel_abundance_old) %>% 
  group_by(combo, passage, ratio, line, category) %>%
  summarize(totalAbundance = sum(rel_abundance_old, na.rm=TRUE)) %>% 
  unique() %>% 
  pivot_wider(names_from=ratio, values_from=totalAbundance) %>% 
  mutate(diff = `0:1`-`1:0`) %>% 
  group_by(combo, passage, category) %>% 
  summarize(avgDiff = mean(diff), minDiff = min(diff), maxDiff = max(diff))
  
dominanceOCUC <- left_join(OCUCRelAbAll, dominanceToParents)

# Plot JSD ratio and total relative abundance of OC ASVs of parent 1.
p_JSDOC <- dominanceOCUC %>%
  filter(passage==5 & category=="overcolonizing") %>% 
  ggplot() +
  geom_point(aes(x=avgDiff, y=avgJSDratio), alpha=0.6, size=0.75) +
  geom_errorbar(aes(x=avgDiff, ymin=minJSDratio, ymax=maxJSDratio), width=0, alpha=0.6, size=0.5) +
  geom_errorbarh(aes(xmin=minDiff, xmax=maxDiff, y=avgJSDratio), height=0, alpha=0.6, size=0.5) +
  #geom_smooth(aes(x=OCSubtracted, y=JSDratio), method="lm", se=FALSE) +
  xlim(-0.2,0.2) +
  ylim(0.5, 5.5) +
  xlab("OC relative abundance diff") +
  ylab("JSD to parents ratio") +
  DEFAULTS.THEME_PRINT
p_JSDOC
save_plot(paste0(plotdir,"/dominanceOC.pdf"), p_JSDOC, base_width=2, base_height=2)

JSDOC <- lm(dominanceOCUC %>% filter(passage==5 & category=="overcolonizing") %>% pull(avgJSDratio) ~ 
              dominanceOCUC %>% filter(passage==5 & category=="overcolonizing") %>% pull(avgDiff))

# Plot JSD ratio and total relative abundance of OC ASVs of parent 1.
p_JSDUC <- dominanceOCUC %>% 
  filter(passage==5 & category=="undercolonizing") %>% 
  ggplot() +
  geom_point(aes(x=avgDiff, y=avgJSDratio), alpha=0.6, size=0.75) +
  geom_errorbar(aes(x=avgDiff, ymin=minJSDratio, ymax=maxJSDratio), width=0, alpha=0.6, size=0.5) +
  geom_errorbarh(aes(xmin=minDiff, xmax=maxDiff, y=avgJSDratio), height=0, alpha=0.6, size=0.5) +
  #geom_smooth(aes(x=UCSubtracted, y=JSDratio), method="lm", se=FALSE) +
  xlim(-0.2,0.2) +
  ylim(0.5, 5.5) +
  xlab("UC relative abundance diff") +
  ylab("JSD to parents ratio") +
  DEFAULTS.THEME_PRINT
p_JSDUC
save_plot(paste0(plotdir,"/dominanceUC.pdf"), p_JSDUC, base_width=2, base_height=2)

JSDUC <- lm(dominanceOCUC %>% filter(passage==5 & category=="undercolonizing") %>% pull(avgJSDratio) ~ 
              dominanceOCUC %>% filter(passage==5 & category=="undercolonizing") %>% pull(avgDiff))

# Analyze correlation of theoretical-actual JSD with total proportion of OC/UC ASVs. --------

# Get 1:1 theoretical-actual JSD.
theoryActualJSD <- summary_JSD_df %>% 
  filter(type=="theoretical-actual" & ratio1=="1:1") %>% 
  group_by(comparison, passage1) %>% 
  summarize(meanJSD = mean(JSD), maxJSD = max(JSD), minJSD = min(JSD)) %>% 
  rename(combo=comparison, passage=passage1)

OCUCRelAbCombined <- theoretical_mixtures_df_categorized %>% 
  filter(category %in% c("overcolonizing","undercolonizing") & ratio=="1:1" & mixtureType=="actual") %>% 
  group_by(combo, passage, line) %>% 
  summarize(totalAbundance = sum(rel_abundance_old, na.rm=TRUE)) %>% 
  ungroup() %>% 
  group_by(combo, passage) %>% 
  summarize(maxAbundance = max(totalAbundance), minAbundance = min(totalAbundance), totalAbundance = mean(totalAbundance))

theoryActualOCUC <- left_join(theoryActualJSD, OCUCRelAbCombined)

# Plot JSD of theoretical-actual comparison and total abundance of OC/UC ASVs at 1:1.
p_theoryActualJSDOCUC <- theoryActualOCUC %>% 
  filter(passage==5) %>% 
  ggplot() +
  geom_point(aes(x=totalAbundance, y=meanJSD), alpha=0.6, size=0.75) +
  geom_errorbar(aes(x=totalAbundance, ymin=minJSD, ymax=maxJSD), width=0, alpha=0.6, size=0.5) +
  geom_errorbarh(aes(xmin=minAbundance, xmax=maxAbundance, y=meanJSD), height=0, alpha=0.6, size=0.5) +
  #geom_smooth(aes(x=totalAbundance, y=meanJSD), method="lm", se=FALSE) +
  #geom_abline(intercept=0.11876, slope=-0.05236) +
  xlim(0,0.3) +
  ylim(0,0.3) +
  xlab("OC/UC relative abundance") +
  ylab("JSD to theory") +
  DEFAULTS.THEME_PRINT
p_theoryActualJSDOCUC
save_plot(paste0(plotdir,"/theoryActualJSDOCUC.pdf"), p_theoryActualJSDOCUC, base_width=2, base_height=2)

theoryActualLm <- lm(theoryActualOCUC$meanJSD ~ theoryActualOCUC$totalAbundance)
