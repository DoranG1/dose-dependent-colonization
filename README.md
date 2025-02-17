# dose-dependent-colonization

The code in this repository performs the analyses described in the manuscript "Competition for shared resources increases dependence on initial population size during coalescence of gut microbial communities," available here: https://www.biorxiv.org/content/10.1101/2023.11.29.569120v2

The code here performs analysis on the sequencing, growth curve, and metabolomic data associated with the manuscript and runs the described consumer-resource model. Individual published figure panels are included for reference alongside each analysis script, as well as relevant data files.

For sequencing data, cleaned and consolidated data text files corresponding to published figures are located in the following directories:\n
For parent communities: sequencing/parentCommunities/data/intermediates/e0012-mixtureDataframe.txt\n
For community coalescence: sequencing/communityCoalescence/data/intermediates/mixtureDataframe.zip\n
For follow-up mixtures (including pairwise, 3-strain, and strain-community mixtures): sequencing/followUpMixtures/combinedCleanedData/combinedCleanedMixtureDataframe.txt

Note the following correspondence of published strain isolate and in vitro community names to internal names used in analysis scripts:\n
XBA:  A1\n
XBB:  A2\n
XCA:  B1\n
XCB:  B2\n
XDA:  C1\n
XDB:  C2\n
XFA:  D1\n
XFB:  D2\n
EnteC2 (E2):  E. faecalis\n
EnteC3 (E3):  E. casseliflavus\n
Strep7 (S7):  L. garvieae\n
Strep17 (S17):  L. lactis\n
Bacte0126 (B126): B. fragilis\n
Tanne0007 (T7): P. goldsteinii\n
