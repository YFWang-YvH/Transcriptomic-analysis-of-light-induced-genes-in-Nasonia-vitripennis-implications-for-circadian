# Transcriptomic-analysis-of-light-induced-genes-in-Nasonia-vitripennis-implications-for-circadian

This folder contains all the scripts needed to reproduce the analysis of RNAseq data from Nasonia vitripennis, published in (journal) - Transcriptomic analysis of light-induced genes in Nasonia vitripennis: implications for circadian light entrainment pathways.

All used packages and software are reported in the manuscript.

The raw RNAseq reads can be found on the European Nucleotide Archive (ENA) under accession no. PRJEB57723. All processed data used to run the analysis can be found on (data archive), including the final gene counts matrix, sample information file, GO annotation table, and all supplementary data files (RNAseq data quality, Genome mapping stats, DEGs list, GO analysis outputs, KEGG analysis ouput, and motif analysis output).

Author
Yifan Wang (YF Wang), ORCID ID: 0000-0002-6541-7435

Step1: Preprocessing RNAseq read data

File for linux to run the bioinformatics analysis process for RNAseq processing including QC, trimming, mapping, and transcript quantification.

Step2: Functional annotation

File for linux to run functional annotation for further process RNAseq data including blast search, interproscan, Pannzer, and uniprot.

Step3: Preprocessing RNAseq count matrix

R script to preprocess RNAseq count matrix, including filtering the data, visualization of the data, PCA analysis, and model comparison.

Step4: Differential expression analysis

R script for differential expression analysis using R package DEseq2.

Step5: Clustering analysis for DEGs 

R script for clustering analysis for DEGs, including hierarchical clustering and time series clustering analysis.

Step6: Functional annotation in R

R script to further tidy up functional annotation from Step 2 in R.

Step7: GO and KEGG functional gene set analysis 

R script for GO and KEGG overrepresentation analysis on DEGs, using R package TopGO.

Step8: Motif analysis

File for linux to run motif analysis for DEGs

Step9: Figures of manuscript 

R script to reproduce all figures reported in the manuscript.
