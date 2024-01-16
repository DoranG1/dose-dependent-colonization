# dose-dependent-colonization

The code in this repository performs the analyses described in the manuscript "Competition for shared resources increases dependence on initial population size during coalescence of gut microbial communities," available here: https://www.biorxiv.org/content/10.1101/2023.11.29.569120v2

The code here processes raw 16S sequencing reads, categorizes ASV colonization behavior, and performs additional custom analyses.

Each sequencing run is located in a separate directory:

•	 community-coalescence: the original eight pairwise community coalescence experiments. The raw data from these mixtures is located at: 210827-e0017-16S/analysis/mixtureDataframe.zip.

•	 followup-mixtures-1: mixtures of L. garvieae/L. lactis and A1/L.lactis. The raw data from these mixtures is located at: 16S-e0043/data/ps_all.txt.gz.

•	 followup-mixtures-2: all other pairwise,  three-strain, and strain-by-community mixtures. The raw data from  these mixtures is located at: 16S-e049-e0050-doseDependence/workflow/out/e0049-e0050/DADA2_output/ps_all.txt.gz.
