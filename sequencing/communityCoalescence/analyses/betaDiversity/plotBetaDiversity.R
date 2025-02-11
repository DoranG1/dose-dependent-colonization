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

p_dominance_JSD_p5_2ex <- dominance_JSD_df %>%
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
  filter(combo %in% c("B1/C1","D1/D2")) %>% 
  mutate(combo = fct_relevel(combo, "B1/C1","D1/D2")) %>% 
  # Remove errorbars for theoreticals and parents (for presentation examples).
  mutate(JSD_1_max = ifelse(lineType=="theoretical" | ratio %in% c("1:0","0:1"), 0, JSD_1_max),
         JSD_1_min = ifelse(lineType=="theoretical" | ratio %in% c("1:0","0:1"), 0, JSD_1_min),
         JSD_2_max = ifelse(lineType=="theoretical" | ratio %in% c("1:0","0:1"), 0, JSD_2_max),
         JSD_2_min = ifelse(lineType=="theoretical" | ratio %in% c("1:0","0:1"), 0, JSD_2_min)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(aes(x=JSD_1_avg, y=JSD_2_avg, color=color), size=0.75) +
  geom_line(aes(x=JSD_1_avg, y=JSD_2_avg, group=color), linetype="dotted", alpha=0.5, size=0.5) +
  geom_errorbar(aes(x=JSD_1_avg, color=color, ymin=JSD_2_min, ymax=JSD_2_max), alpha=0.5, size=0.5) +
  geom_errorbarh(aes(y=JSD_2_avg, color=color, xmin=JSD_1_min, xmax=JSD_1_max), alpha=0.5, size=0.5) +
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
  #facet_wrap(~combo, nrow=2, scales="free") +
  facet_wrap(~combo, scales="free") +
  theme(strip.text = element_text(face="bold"))
p_dominance_JSD_p5_2ex
save_plot(paste0(plotdir, "/dominance_JSD_p5_all.pdf"), 
          p_dominance_JSD_p5_all, base_width=6, base_height=3)
#save_plot(paste0(plotdir, "/dominance_JSD_p5_2ex.pdf"),
#                 p_dominance_JSD_p5_2ex, base_width=3, base_height=1.5)
save_plot(paste0(plotdir, "/dominance_JSD_p5_2ex.png"), p_dominance_JSD_p5_2ex,
          base_width=3, base_height=1.5)


dominanceJSDLegend <- get_legend(p_dominance_JSD_p5_all)
save_plot(paste0(plotdir, "/dominanceJSDLegend.pdf"), dominanceJSDLegend, base_width=0.75, base_height=1.2)

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
