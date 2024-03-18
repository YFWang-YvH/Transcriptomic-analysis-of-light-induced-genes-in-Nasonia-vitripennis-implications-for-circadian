# Transcriptomic-analysis-of-light-induced-genes-in-Nasonia-vitripennis-implications-for-circadian

## Information

This folder contains all the scripts needed to reproduce the analysis of RNAseq data from Nasonia vitripennis, published in the special issue Biological Rhythms and Molecular Clockworks in Physiology and Pathology in Biology - Transcriptomic analysis of light-induced genes in Nasonia vitripennis: implications for circadian light entrainment pathways.

All used packages and software are reported in the manuscript.

The raw RNAseq reads can be found on the European Nucleotide Archive (ENA) under accession no. PRJEB57723. All processed data used to run the analysis can be found on (data archive), including the final gene counts matrix, sample information file, GO annotation table, and all supplementary data files (RNAseq data quality, Genome mapping stats, DEGs list, GO analysis outputs, KEGG analysis output, and motif analysis output).

Author:

Yifan Wang (YF Wang), ORCID ID: 0000-0002-6541-7435


## Step1: Preprocessing RNAseq read data 
[Preprocessing RNAseq read data](https://github.com/YFWang-YvH/Transcriptomic-analysis-of-light-induced-genes-in-Nasonia-vitripennis-implications-for-circadian/blob/main/1_preprocessing%20RNAseq%20read%20data)


File for linux to run the bioinformatics analysis process for RNAseq processing including QC, trimming, mapping, and transcript quantification.

## Step2: Functional annotation
[Functional annotation](https://github.com/YFWang-YvH/Transcriptomic-analysis-of-light-induced-genes-in-Nasonia-vitripennis-implications-for-circadian/blob/main/2_functional%20annotation)


File for linux to run functional annotation for further process RNAseq data including blast search, **interproscan**, **Pannzer**, and **uniprot**.

## Step3: Preprocessing RNAseq count matrix
[Preprocessing RNAseq count matrix](https://github.com/YFWang-YvH/Transcriptomic-analysis-of-light-induced-genes-in-Nasonia-vitripennis-implications-for-circadian/blob/main/3_preprocessing%20RNAseq%20count%20matrix%20R%20script.R)


R script to preprocess RNAseq count matrix, including filtering the data, visualization of the data, PCA analysis, and model comparison.

## Step4: Differential expression analysis
[Differential expression analysis](https://github.com/YFWang-YvH/Transcriptomic-analysis-of-light-induced-genes-in-Nasonia-vitripennis-implications-for-circadian/blob/main/4_differential%20expression%20analysis%20R%20script.R)


R script for differential expression analysis using R package **DEseq2**.

## Step5: Clustering analysis for DEGs 
[Clustering analysis for DEGs](https://github.com/YFWang-YvH/Transcriptomic-analysis-of-light-induced-genes-in-Nasonia-vitripennis-implications-for-circadian/blob/main/5_clustering%20analysis%20for%20DEGs%20R%20script.R)


R script for clustering analysis for DEGs, including hierarchical clustering and time series clustering analysis.

## Step6: Functional annotation in R
[Functional annotation in R](https://github.com/YFWang-YvH/Transcriptomic-analysis-of-light-induced-genes-in-Nasonia-vitripennis-implications-for-circadian/blob/main/6_fuctional%20annotation%20R%20script.R)


R script to further tidy up functional annotation from **Step 2** in R.

## Step7: GO and KEGG functional gene set analysis 
[GO and KEGG functional gene set analysis](https://github.com/YFWang-YvH/Transcriptomic-analysis-of-light-induced-genes-in-Nasonia-vitripennis-implications-for-circadian/blob/main/7_GO%20and%20KEGG%20functional%20gene%20set%20analysis%20R%20script.R)


R script for GO and KEGG overrepresentation analysis on DEGs, using R package **TopGO**.

## Step8: Motif analysis
[Motif analysis](https://github.com/YFWang-YvH/Transcriptomic-analysis-of-light-induced-genes-in-Nasonia-vitripennis-implications-for-circadian/blob/main/8_motif%20analysis)


File for linux to run motif analysis for DEGs

## Step9: Figures of manuscript 
[Figures of manuscript](https://github.com/YFWang-YvH/Transcriptomic-analysis-of-light-induced-genes-in-Nasonia-vitripennis-implications-for-circadian/blob/main/9_figures%20of%20manuscript%20R%20script.R)


R script to reproduce all figures reported in the manuscript.
