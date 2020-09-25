This project contains scripts for the analysis of screening data. The screening was performed on the E. coli KEIO collection of deletions mutants. Each mutant was grown in the presence or absence of the antimicrobial peptide (AMP) TAT-RasGAP.

The dataset of screening results (dataset1) is used as input for the analysis. It contains the OD590 measurements of the strains after 0h (T0), 1h30 (T1), 3h (T2), 6h (T3) and 24h (T4) incubation with or without TAT-RasGAP. 

The screening_selection.R script contains the information about the calculation of the normalized OD590 and the Fold change calculation. It also contains the script to plot the distribution of fold change and the selection of resistant, hypersensitive and more sensitive strains.

The screening_GO.R and screening_KEGG.R scripts contain the information about GO and KEGG enrichment analysis for the set of hypersensitive genes. 

The results of this analysis are in a manuscript on bioRxiv andsubmitted for publication.



