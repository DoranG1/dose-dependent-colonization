library(tidyverse)

# Import taxa color palette.
palette <- read.table("config/palettes/KCHcolors-gg.txt",
                      header=TRUE, stringsAsFactors = FALSE)

# Convert taxa to Silva name format.
palette <- palette %>%
  mutate(taxaSilva=gsub(";",".",taxa),
         taxaSilva=gsub("k__","",taxaSilva),
         taxaSilva=gsub("p__","",taxaSilva),
         taxaSilva=gsub("c__","",taxaSilva),
         taxaSilva=gsub("o__","",taxaSilva),
         taxaSilva=gsub("f__","",taxaSilva),
         taxaSilva=gsub("g__","",taxaSilva))
palette <- palette %>%
  mutate(taxaSilvaCorrected=taxaSilva)

# Manually correct the names of common taxa with minor spelling differences (level 2).
palette <- palette %>%
  mutate(taxaSilvaCorrected=
           ifelse(taxaSilva=="Bacteria.Bacteroidetes","Bacteria.Bacteroidota",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=
           ifelse(taxaSilva=="Bacteria.Verrucomicrobia","Bacteria.Verrucomicrobiota",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=
           ifelse(taxaSilva=="Bacteria.Fusobacteria","Bacteria.Fusobacteriota",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=
           ifelse(taxaSilva=="Bacteria.Actinobacteria","Bacteria.Actinobacteriota",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=
           ifelse(taxaSilva=="Bacteria.Acidobacteria","Bacteria.Acidobacteriota",taxaSilvaCorrected))
# Add the closest color for taxa that do not have clear analogues in Green Genes (level 2).
palette <- palette %>%
  rbind(c(2, "", "", "#8DCAB9", 0.5517241, 0.7931034, 0.7241379,
          "Bacteria.Desulfobacterota","Bacteria.Desulfobacterota")) %>%
  rbind(c(2, "", "", "#46003E", 0.2758621, 0.0000000, 0.2413793,
          "Bacteria.Campylobacterota","Bacteria.Campylobacterota")) %>%
  rbind(c(2, "", "", "#0289b6", 2/255, 137/255, 182/255,
          "Bacteria.Patescibacteria","Bacteria.Patescibacteria"))

# Manually correct the names of common taxa with minor spelling differences (level 3).
palette <- palette %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Bacteroidota.Bacteroidia",
                                   "Bacteria.Bacteroidetes.Bacteroidia",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Verrucomicrobiota.Verrucomicrobiae",
                                   "Bacteria.Verrucomicrobia.Verrucomicrobiae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Fusobacteriota.Fusobacteriia",
                                   "Bacteria.Fusobacteria.Fusobacteriia",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Actinobacteriota.Actinobacteria",
                                   "Bacteria.Actinobacteria.Actinobacteria",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Actinobacteriota.Coriobacteriia",
                                   "Bacteria.Actinobacteria.Coriobacteriia",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Verrucomicrobiota.Lentisphaeria",
                                   "Bacteria.Lentisphaerae.[Lentisphaeria]",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Verrucomicrobiota.Lentisphaeria",
                                   "Bacteria.Lentisphaerae.[Lentisphaeria]",taxaSilvaCorrected))
# Add the closest color for taxa that do not have clear analogues in Green Genes (level 3).
# Match the phylum-level colors when only one class is observed in the phylum,
# or when the phylum is extremely rare.
palette <- palette %>%
  rbind(c(3, "", "", "#8DCAB9", 0.5517241, 0.7931034, 0.7241379,
          "Bacteria.Desulfobacterota.Desulfovibrionia","Bacteria.Desulfobacterota.Desulfovibrionia")) %>%
  rbind(c(3, "", "", "#46003E", 0.2758621, 0.0000000, 0.2413793,
          "Bacteria.Campylobacterota.Campylobacteria","Bacteria.Campylobacterota.Campylobacteria")) %>%
  rbind(c(3, "", "", "#0289b6", 2/255, 137/255, 182/255,
          "Bacteria.Patescibacteria.Gracilibacteria","Bacteria.Patescibacteria.Gracilibacteria")) %>%
  rbind(c(3, "", "", "#0289b6", 2/255, 137/255, 182/255,
          "Bacteria.Patescibacteria.Saccharimonadia","Bacteria.Patescibacteria.Saccharimonadia")) %>%
  rbind(c(3, "", "", "#70cd70", 112/255, 205/255, 112/255,
          "Bacteria.Cyanobacteria.Cyanobacteriia","Bacteria.Cyanobacteria.Cyanobacteriia")) %>%
  rbind(c(3, "", "", "#376637", 55/255, 102/255, 55/255,
          "Bacteria.Cyanobacteria.Vampirivibrionia","Bacteria.Cyanobacteria.Vampirivibrionia")) %>%
  rbind(c(3, "", "", "#cb5d6b", 203/255, 93/255, 107/255,
          "Bacteria.Acidobacteriota.Subgroup_19","Bacteria.Acidobacteriota.Subgroup_19")) %>%
  rbind(c(3, "", "", "#eb93ee", 235/255, 147/255, 238/255,
          "Bacteria.Firmicutes.Negativicutes","Bacteria.Firmicutes.Negativicutes")) %>%
  rbind(c(3, "", "", "#5e3a5f", 94/255, 58/255, 95/255,
          "Bacteria.Firmicutes.Incertae_Sedis","Bacteria.Firmicutes.Incertae_Sedis")) %>%
  rbind(c(3, "", "", "#beb0bf", 190/255, 176/255, 191/255,
          "Bacteria.Firmicutes.uncultured","Bacteria.Firmicutes.uncultured"))


# Manually correct the names of common taxa with minor spelling differences (level 4).
palette <- palette %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Actinobacteriota.Actinobacteria.Actinomycetales",
                                   "Bacteria.Actinobacteria.Actinobacteria.Actinomycetales",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Actinobacteriota.Actinobacteria.Bifidobacteriales",
                                   "Bacteria.Actinobacteria.Actinobacteria.Bifidobacteriales",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Actinobacteria.Coriobacteriia.Coriobacteriales",
                                   "Bacteria.Actinobacteriota.Coriobacteriia.Coriobacteriales",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales",
                                   "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Bacteroidetes.Cytophagia.Cytophagales",
                                   "Bacteria.Bacteroidota.Bacteroidia.Cytophagales",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Firmicutes.Erysipelotrichi.Erysipelotrichales",
                                   "Bacteria.Firmicutes.Bacilli.Erysipelotrichales",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Fusobacteria.Fusobacteriia.Fusobacteriales",
                                   "Bacteria.Fusobacteriota.Fusobacteriia.Fusobacteriales",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales",
                                   "Bacteria.Verrucomicrobia.Verrucomicrobiae.Verrucomicrobiales",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Verrucomicrobiota.Verrucomicrobiae.Verrucomicrobiales",
                                   "Bacteria.Firmicutes.Bacilli.Erysipelotrichales",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Firmicutes.Erysipelotrichi.Erysipelotrichales",
                                   "Bacteria.Firmicutes.Bacilli.Erysipelotrichales",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Lentisphaerae.[Lentisphaeria].Victivallales",
                                   "Bacteria.Verrucomicrobiota.Lentisphaeria.Victivallales",taxaSilvaCorrected)) %>%
  # This next one is weird - Silva places Burkholderiales in gamma rather than betaproteobacteria.
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Proteobacteria.Betaproteobacteria.Burkholderiales",
                                   "Bacteria.Proteobacteria.Gammaproteobacteria.Burkholderiales",taxaSilvaCorrected))
# Add the closest color for taxa that do not have clear analogues in Green Genes (level 3).
# Match the phylum/class-level colors when only one class is observed in the phylum,
# or when the phylum/class is extremely rare.
palette <- palette %>%
  rbind(c(4, "", "", "#8DCAB9", 0.5517241, 0.7931034, 0.7241379,
          "Bacteria.Desulfobacterota.Desulfovibrionia.Desulfovibrionales",
          "Bacteria.Desulfobacterota.Desulfovibrionia.Desulfovibrionales")) %>%
  rbind(c(4, "", "", "#46003E", 0.2758621, 0.0000000, 0.2413793,
          "Bacteria.Campylobacterota.Campylobacteria.Campylobacterales",
          "Bacteria.Campylobacterota.Campylobacteria.Campylobacterales")) %>%
  rbind(c(4, "", "", "#0289b6", 2/255, 137/255, 182/255,
          "Bacteria.Patescibacteria.Gracilibacteria.Absconditabacteriales_.SR1.",
          "Bacteria.Patescibacteria.Gracilibacteria.Absconditabacteriales_.SR1.")) %>%
  rbind(c(4, "", "", "#0289b6", 2/255, 137/255, 182/255,
          "Bacteria.Patescibacteria.Saccharimonadia.Saccharimonadales",
          "Bacteria.Patescibacteria.Saccharimonadia.Saccharimonadales")) %>%
  rbind(c(4, "", "", "#70cd70", 112/255, 205/255, 112/255,
          "Bacteria.Cyanobacteria.Cyanobacteriia.Chloroplast",
          "Bacteria.Cyanobacteria.Cyanobacteriia.Chloroplast")) %>%
  rbind(c(4, "", "", "#376637", 55/255, 102/255, 55/255,
          "Bacteria.Cyanobacteria.Vampirivibrionia.Gastranaerophilales",
          "Bacteria.Cyanobacteria.Vampirivibrionia.Gastranaerophilales")) %>%
  rbind(c(4, "", "", "#eb93ee", 235/255, 147/255, 238/255,
          "Bacteria.Firmicutes.Negativicutes.Acidaminococcales",
          "Bacteria.Firmicutes.Negativicutes.Acidaminococcales")) %>%
  rbind(c(4, "", "", "#eb93ee", 235/255, 147/255, 238/255,
          "Bacteria.Firmicutes.Negativicutes.Veillonellales.Selenomonadales",
          "Bacteria.Firmicutes.Negativicutes.Veillonellales.Selenomonadales")) %>%
  rbind(c(4, "", "", "#0A0A0A", 0.0384615384615385, 0.0384615384615385, 0.0384615384615385,
          "Bacteria.Actinobacteriota.Actinobacteria.Corynebacteriales",
          "Bacteria.Actinobacteriota.Actinobacteria.Corynebacteriales")) %>%
  rbind(c(4, "", "", "#093E4F", 0.0344827586206897, 0.241379310344828, 0.310344827586207,
          "Bacteria.Actinobacteriota.Actinobacteria.Micrococcales",
          "Bacteria.Actinobacteriota.Actinobacteria.Micrococcales")) %>%
  rbind(c(4, "", "", "#454545", 0.269230769230769, 0.269230769230769, 0.269230769230769,
          "Bacteria.Firmicutes.Bacilli.Staphylococcales",
          "Bacteria.Firmicutes.Bacilli.Staphylococcales")) %>%
  rbind(c(4, "", "", "#4E4E4E", 0.307692307692308, 0.307692307692308, 0.307692307692308,
          "Bacteria.Firmicutes.Clostridia.Christensenellales",
          "Bacteria.Firmicutes.Clostridia.Christensenellales")) %>%
  rbind(c(4, "", "", "#A86666", 0.66, 0.4, 0.4,
          "Bacteria.Firmicutes.Clostridia.Eubacteriales",
          "Bacteria.Firmicutes.Clostridia.Eubacteriales"))  %>%
  rbind(c(4, "", "", "#66A61E", 0.4, 0.650980392156863, 0.117647058823529,
          "Bacteria.Firmicutes.Clostridia.Lachnospirales",
          "Bacteria.Firmicutes.Clostridia.Lachnospirales")) %>%
  rbind(c(4, "", "", "#8D6100", 0.551724137931034, 0.379310344827586, 0,
          "Bacteria.Firmicutes.Clostridia.Oscillospirales",
          "Bacteria.Firmicutes.Clostridia.Oscillospirales")) %>%
  rbind(c(4, "", "", "#230023", 0.137931034482759, 0, 0.137931034482759,
          "Bacteria.Firmicutes.Clostridia.Peptococcales",
          "Bacteria.Firmicutes.Clostridia.Peptococcales")) %>%
  rbind(c(4, "", "", "#FADBDB", 0.98, 0.86, 0.86,
          "Bacteria.Firmicutes.Clostridia.Peptostreptococcales.Tissierellales",
          "Bacteria.Firmicutes.Clostridia.Peptostreptococcales.Tissierellales")) %>%
  rbind(c(4, "", "", "#FADBDB", 0.98, 0.86, 0.86,
          "Bacteria.Firmicutes.Clostridia.Peptostreptococcales.Tissierellales",
          "Bacteria.Firmicutes.Clostridia.Peptostreptococcales.Tissierellales")) %>%
  rbind(c(4, "", "", "#7010e3", 112/255, 16/255, 227/255,
          "Bacteria.Firmicutes.Bacilli.Izemoplasmatales",
          "Bacteria.Firmicutes.Bacilli.Izemoplasmatales")) %>%
  rbind(c(4, "", "", "#f998cf", 249/255, 152/255, 207/255,
          "Bacteria.Firmicutes.Bacilli.RF39",
          "Bacteria.Firmicutes.Bacilli.RF39")) %>%
  rbind(c(4, "", "", "#943830", 148/255, 56/255, 48/255,
          "Bacteria.Firmicutes.Clostridia.Clostridia_UCG.014",
          "Bacteria.Firmicutes.Clostridia.Clostridia_UCG.014")) %>%
  rbind(c(4, "", "", "#122185", 18/255, 33/255, 133/255,
          "Bacteria.Firmicutes.Clostridia.Clostridia_vadinBB60_group",
          "Bacteria.Firmicutes.Clostridia.Clostridia_vadinBB60_group")) %>%
  rbind(c(4, "", "", "#87d7f6", 135/255, 215/255, 246/255,
          "Bacteria.Firmicutes.Clostridia.Monoglobales",
          "Bacteria.Firmicutes.Clostridia.Monoglobales"))

# Manually correct the names of common taxa with minor spelling differences (level 5).
palette <- palette %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Actinobacteria.Actinobacteria.Bifidobacteriales.Bifidobacteriaceae",
                                   "Bacteria.Actinobacteriota.Actinobacteria.Bifidobacteriales.Bifidobacteriaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Bacteroidaceae",
                                   "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Bacteroidaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.[Barnesiellaceae]",
                                   "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Barnesiellaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Marinilabiaceae",
                                   "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Marinifilaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Prevotellaceae",
                                   "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Prevotellaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae",
                                   "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Rikenellaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Firmicutes.Erysipelotrichi.Erysipelotrichales.Erysipelotrichaceae",
                                   "Bacteria.Firmicutes.Bacilli.Erysipelotrichales.Erysipelotrichaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Firmicutes.Clostridia.Clostridiales.Eubacteriaceae",
                                   "Bacteria.Firmicutes.Clostridia.Eubacteriales.Eubacteriaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae",
                                   "Bacteria.Firmicutes.Clostridia.Lachnospirales.Lachnospiraceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae",
                                   "Bacteria.Firmicutes.Clostridia.Oscillospirales.Ruminococcaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Firmicutes.Clostridia.Clostridiales.Peptococcaceae",
                                   "Bacteria.Firmicutes.Clostridia.Peptococcales.Peptococcaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae",
                                   "Bacteria.Firmicutes.Negativicutes.Veillonellales.Selenomonadales.Veillonellaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Fusobacteria.Fusobacteriia.Fusobacteriales.Fusobacteriaceae",
                                   "Bacteria.Fusobacteriota.Fusobacteriia.Fusobacteriales.Fusobacteriaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae",
                                   "Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacterales.Enterobacteriaceae",taxaSilvaCorrected)) %>%
  mutate(taxaSilvaCorrected=ifelse(taxaSilva=="Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae",
                                   "Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacterales.Enterobacteriaceae",taxaSilvaCorrected))
# Add the closest color for taxa that do not have clear analogues in Green Genes (level 3).
# Match the phylum/class-level colors when only one class is observed in the phylum,
# or when the phylum/class is extremely rare.
palette <- palette %>%
  rbind(c(5, "", "", "#8DCAB9", 0.5517241, 0.7931034, 0.7241379,
          "Bacteria.Desulfobacterota.Desulfovibrionia.Desulfovibrionales.Desulfovibrionaceae",
          "Bacteria.Desulfobacterota.Desulfovibrionia.Desulfovibrionales.Desulfovibrionaceae")) %>%
  rbind(c(5, "", "", "#eb93ee", 235/255, 147/255, 238/255,
          "Bacteria.Firmicutes.Negativicutes.Acidaminococcales.Acidaminococcaceae",
          "Bacteria.Firmicutes.Negativicutes.Acidaminococcales.Acidaminococcaceae")) %>%
  rbind(c(5, "", "", "#8D6100", 0.551724137931034, 0.379310344827586, 0,
          "Bacteria.Firmicutes.Clostridia.Oscillospirales..Eubacterium._coprostanoligenes_group",
          "Bacteria.Firmicutes.Clostridia.Oscillospirales..Eubacterium._coprostanoligenes_group")) %>%
  rbind(c(5, "", "", "#8D6100", 0.551724137931034, 0.379310344827586, 0,
          "Bacteria.Firmicutes.Clostridia.Oscillospirales.Butyricicoccaceae",
          "Bacteria.Firmicutes.Clostridia.Oscillospirales.Butyricicoccaceae")) %>%
  rbind(c(5, "", "", "#8D6100", 0.551724137931034, 0.379310344827586, 0,
          "Bacteria.Firmicutes.Clostridia.Oscillospirales.Ethanoligenenaceae",
          "Bacteria.Firmicutes.Clostridia.Oscillospirales.Ethanoligenenaceae")) %>%
  rbind(c(5, "", "", "#8D6100", 0.551724137931034, 0.379310344827586, 0,
          "Bacteria.Firmicutes.Clostridia.Oscillospirales.Hydrogenoanaerobacterium",
          "Bacteria.Firmicutes.Clostridia.Oscillospirales.Hydrogenoanaerobacterium")) %>%
  rbind(c(5, "", "", "#8D6100", 0.551724137931034, 0.379310344827586, 0,
          "Bacteria.Firmicutes.Clostridia.Oscillospirales.Oscillospiraceae",
          "Bacteria.Firmicutes.Clostridia.Oscillospirales.Oscillospiraceae")) %>%
  rbind(c(5, "", "", "#daa520", 218/255, 165/255, 32/255,
          "Bacteria.Firmicutes.Clostridia.Peptostreptococcales.Tissierellales.Anaerovoracaceae",
          "Bacteria.Firmicutes.Clostridia.Peptostreptococcales.Tissierellales.Anaerovoracaceae")) %>%
  rbind(c(5, "", "", "#daa520", 218/255, 165/255, 32/255,
          "Bacteria.Firmicutes.Clostridia.Peptostreptococcales.Tissierellales.Family_XI",
          "Bacteria.Firmicutes.Clostridia.Peptostreptococcales.Tissierellales.Family_XI")) %>%
  rbind(c(5, "", "", "#daa520", 218/255, 165/255, 32/255,
          "Bacteria.Firmicutes.Clostridia.Peptostreptococcales.Tissierellales.Peptostreptococcaceae",
          "Bacteria.Firmicutes.Clostridia.Peptostreptococcales.Tissierellales.Peptostreptococcaceae")) %>%
  rbind(c(5, "", "", "#123561", 0.0689655172413793, 0.206896551724138, 0.379310344827586,
          "Bacteria.Actinobacteriota.Coriobacteriia.Coriobacteriales.Eggerthellaceae",
          "Bacteria.Actinobacteriota.Coriobacteriia.Coriobacteriales.Eggerthellaceae")) %>%
  rbind(c(5, "", "", "#FFFF72", 1, 1, 0.4482759,
          "Bacteria.Firmicutes.Bacilli.Erysipelotrichales.Erysipelatoclostridiaceae",
          "Bacteria.Firmicutes.Bacilli.Erysipelotrichales.Erysipelatoclostridiaceae")) %>%
  rbind(c(5, "", "", "#B0848D", 0.689655172413793, 0.517241379310345, 0.551724137931034,
          "Bacteria.Proteobacteria.Gammaproteobacteria.Burkholderiales.Sutterellaceae",
          "Bacteria.Proteobacteria.Gammaproteobacteria.Burkholderiales.Sutterellaceae")) %>%
  rbind(c(5, "", "", "#88da1b", 136/255, 218/255, 27/255,
          "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Muribaculaceae",
          "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Muribaculaceae")) %>%
  rbind(c(5, "", "", "#9ec3f1", 158/255, 195/255, 241/255,
          "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Tannerellaceae",
          "Bacteria.Bacteroidota.Bacteroidia.Bacteroidales.Tannerellaceae")) %>%
  rbind(c(5, "", "", "#960ce3", 150/255, 12/255, 227/255,
          "Bacteria.Firmicutes.Clostridia.Monoglobales.Monoglobaceae",
          "Bacteria.Firmicutes.Clostridia.Monoglobales.Monoglobaceae")) %>%
  rbind(c(5, "", "", "#1B9E77", 0.10588235, 0.61960784, 0.46666667,
          "Bacteria.Verrucomicrobiota.Verrucomicrobiae.Verrucomicrobiales.Akkermansiaceae",
          "Bacteria.Verrucomicrobiota.Verrucomicrobiae.Verrucomicrobiales.Akkermansiaceae")) %>%
  rbind(c(5, "", "", "#6A3D9A", 106/255, 61/255, 154/255,
          "Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacterales.Hafniaceae",
          "Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacterales.Hafniaceae"))


# Export Silva palette in same format as GG palette.
# Eliminate duplicate entries to prevent formatting issues.
write.table(palette %>% mutate(taxa=taxaSilvaCorrected) %>%
              dplyr::select(level, taxa, taxashort, hex, R, G, B),
            "config/palettes/KCHcolors-Silva-partial.txt",
            row.names=FALSE, quote=TRUE)

# Export full GG to Silva conversion.
write.table(palette,"config/palettes/GGtoSilva-partial.txt",
            row.names=FALSE, quote=TRUE)
