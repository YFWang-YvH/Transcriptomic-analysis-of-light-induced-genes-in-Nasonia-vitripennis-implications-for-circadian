###Using StringTie with DESeq2 to analyze Nasonia RNAseq data
rm(list=ls()) #clear all objects
gc()

#dataset: 
#Nasonia RNAseq data
#experimental light treatment and dark control
#4 different time point: 0.5h 1h 2h 4h
#8 sample groups
#3 replicates per experimental sample group, 2 replicates per control sample group
#20 samples in total

#load packages 
library("vctrs")
library("ellipsis")
library(DESeq2)
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library("ashr")
library(ggplotify)
library(ggplot2)

#load gene/transcript count matrix and labels
countData <- as.matrix(read.csv("gene_count_matrix_annotatedstringtie.csv", row.names = "gene_id"))
colData <- read.csv("sampleinfo.csv", row.names = 1)

#adjust colData to factors
colData$LaneID <- as.factor(colData$LaneID)
colData$Batch <- as.factor(colData$Batch)
colData$RNAextractionround <- as.factor(colData$RNAextractionround)
colData$Treatment <- as.factor(colData$Treatment)
colData$TreatmentDuration <- as.factor(colData$TreatmentDuration)
colData$Generationnumber <- as.factor(colData$Generationnumber)

########differential expression analysis########
ddsbatch <- DESeqDataSetFromMatrix(countData = countData,
                                   colData = colData,
                                   design= ~ Batch + Treatment*TreatmentDuration)

nrow(ddsbatch) ##15904

ddsbatch <- DESeq(ddsbatch) #15904
ddsbatchClean <- ddsbatch[which(mcols(ddsbatch)$betaConv),]

ddsbatch <- ddsbatchClean
res <- results(ddsbatch, alpha=0.05,cooksCutoff = FALSE)

res
sum(res$padj < 0.05, na.rm=TRUE) #71

levels(ddsbatch$Treatment)
levels(ddsbatch$TreatmentDuration)

#manually do the DESeq contrasts and for an imbalanced design (with less replicates for some samples)
#DESeq defines the model matrix
#1)get the model matrix
mod_mat <- model.matrix(design(ddsbatch), colData(ddsbatch))
mod_mat

#2)define coefficient vectors for each condition
LT_0.5H <- colMeans(mod_mat[ddsbatch$Treatment == "LT" & ddsbatch$TreatmentDuration == "0.5",])
LT_1H <- colMeans(mod_mat[ddsbatch$Treatment == "LT" & ddsbatch$TreatmentDuration == "1",])
LT_2H <- colMeans(mod_mat[ddsbatch$Treatment == "LT" & ddsbatch$TreatmentDuration == "2",])
LT_4H <- colMeans(mod_mat[ddsbatch$Treatment == "LT" & ddsbatch$TreatmentDuration == "4",])
DC_0.5H <- colMeans(mod_mat[ddsbatch$Treatment == "DC" & ddsbatch$TreatmentDuration == "0.5",])
DC_1H <- colMeans(mod_mat[ddsbatch$Treatment == "DC" & ddsbatch$TreatmentDuration == "1",])
DC_2H <- colMeans(mod_mat[ddsbatch$Treatment == "DC" & ddsbatch$TreatmentDuration == "2",])
DC_4H <- colMeans(mod_mat[ddsbatch$Treatment == "DC" & ddsbatch$TreatmentDuration == "4",])
LT <- colMeans(mod_mat[ddsbatch$Treatment == "LT",])
DC <- colMeans(mod_mat[ddsbatch$Treatment == "DC",])

#3)define any contrast of interest from these vectors
#LT VS DC (AT EACH TIME POINTS) 
#OUTLIERS, a diagnostic test for outliers called Cook's distance, measure of how much a single sample is influencing the fitted coefficients for a gene,
#and a large value of Cook's distance is intended to indicate an outlier count. at least 3 replicates are required for flagging, as it is difficult to judge which sample might be an outlier with only 2 replicates
#this filtering can be turned off with results(dds, cooksCutoff=FALSE)
#when there are very many outliers, it might make more sense to turn off the outlier filtering and perform manual inspection
res0.5H <- results(ddsbatch, contrast = (LT_0.5H - DC_0.5H), alpha = 0.05, cooksCutoff=FALSE) 
res1H <- results(ddsbatch,contrast = (LT_1H - DC_1H), alpha = 0.05, cooksCutoff=FALSE)
res2H <- results(ddsbatch, contrast = (LT_2H - DC_2H), alpha = 0.05, cooksCutoff=FALSE)
res4H <- results(ddsbatch, contrast = (LT_4H - DC_4H), alpha = 0.05, cooksCutoff=FALSE)

#how many differentially expressed genes are there at FDR < 0.05

sum(res0.5H$padj < 0.05, na.rm = TRUE) #277 genes
sum(res1H$padj < 0.05, na.rm = TRUE) #381 genes
sum(res2H$padj < 0.05, na.rm = TRUE) #1432 genes
sum(res4H$padj < 0.05, na.rm = TRUE) #147 genes

summary(res0.5H) 
summary(res1H) 
summary(res2H)
summary(res4H)

#It is more useful visualize the MAplot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low countgenes
#without requiring arbitrary filtering thresholds on low count genes

res0.5Hashr <- lfcShrink(ddsbatch,contrast = (LT_0.5H - DC_0.5H), type = "ashr",alpha = 0.05, res=res0.5H)
res1Hashr <- lfcShrink(ddsbatch,contrast = (LT_1H - DC_1H), type = "ashr",alpha = 0.05, res=res1H)
res2Hashr <- lfcShrink(ddsbatch,contrast = (LT_2H - DC_2H), type = "ashr",alpha = 0.05, res=res2H)
res4Hashr <- lfcShrink(ddsbatch,contrast = (LT_4H - DC_4H), type = "ashr",alpha = 0.05, res=res4H)
par(mfrow=c(1,4), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res0.5Hashr, xlim=xlim, ylim=ylim, main="LT0.5HvsDC0.5H")
plotMA(res1Hashr, xlim=xlim, ylim=ylim, main="LT1HvsDC1H")
plotMA(res2Hashr, xlim=xlim, ylim=ylim, main="LT2HvsDC2H")
plotMA(res4Hashr, xlim=xlim, ylim=ylim, main="LT4HvsDC4H")

summary(res0.5Hashr) 
summary(res1Hashr) 
summary(res2Hashr)
summary(res4Hashr)

#how many differentially expressed genes are there at FDR < 0.05

sum(res0.5Hashr$padj < 0.05, na.rm = TRUE) #277 genes
sum(res1Hashr$padj < 0.05, na.rm = TRUE) #381 genes
sum(res2Hashr$padj < 0.05, na.rm = TRUE) #1432 genes
sum(res4Hashr$padj < 0.05, na.rm = TRUE) #147 genes

##################################
#data filtering
# first remove the filtered genes (FDR=NA) and create a -log10(FDR) column
dres0.5Hashrfilt <- dres0.5Hashr %>% 
  filter(!is.na(FDR)) %>% 
  mutate(`-log10(FDR)` = -log10(FDR))
dres0.5Hashrfilt <- dres0.5Hashrfilt[order(dres0.5Hashrfilt$FDR),] #12712

dres1Hashrfilt <- dres1Hashr %>% 
  filter(!is.na(FDR)) %>% 
  mutate(`-log10(FDR)` = -log10(FDR))
dres1Hashrfilt <- dres1Hashrfilt[order(dres1Hashrfilt$FDR),]#13308

dres2Hashrfilt <- dres2Hashr %>% 
  filter(!is.na(FDR)) %>% 
  mutate(`-log10(FDR)` = -log10(FDR))
dres2Hashrfilt <- dres2Hashrfilt[order(dres2Hashrfilt$FDR),]#12116

dres4Hashrfilt <- dres4Hashr %>% 
  filter(!is.na(FDR)) %>% 
  mutate(`-log10(FDR)` = -log10(FDR))
dres4Hashrfilt <- dres4Hashrfilt[order(dres4Hashrfilt$FDR),]#13010

#filter only significant genes out
x <- -log10(0.05)
dres0.5Hashrsig <- dres0.5Hashrfilt %>% 
  filter(`-log10(FDR)` > x)
dres1Hashrsig <- dres1Hashrfilt %>% 
  filter(`-log10(FDR)` > x)
dres2Hashrsig <- dres2Hashrfilt %>% 
  filter(`-log10(FDR)` > x)
dres4Hashrsig <- dres4Hashrfilt %>% 
  filter(`-log10(FDR)` > x)

#pull genesID together and make list for venn diagram

#figure 2 in the manuscript
LTvsDC0.5H <- dres0.5Hashrsig$GeneID
LTvsDC1H <- dres1Hashrsig$GeneID
LTvsDC2H <- dres2Hashrsig$GeneID
LTvsDC4H <- dres4Hashrsig$GeneID

venn_list <- list("30 mins" = LTvsDC0.5H, "4 hrs"=LTvsDC4H, "1 hr" = LTvsDC1H, "2 hrs" = LTvsDC2H)

library("VennDiagram")
# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

# Four dimension Venn plot
# Change category names
# Change fill color
# Further customization
# Export graph into ppt
p1 <- as.ggplot(~display_venn(
  venn_list,
  #category.names = c("30 mins" , "4 hrs" , "1 hr", "2 hrs"),
  # Circles
  lwd = 2,
  #lty = 'blank',
  fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
  # Numbers
  cex = 1.5,
  fontface = "italic",
  # Set names
  cat.cex = 2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.dist = c(0.055, 0.055, 0.1, 0.1)
))

#making the DEG list
DEgenelist <- rbind(dres0.5Hashrsig,
                    dres1Hashrsig,
                    dres2Hashrsig,
                    dres4Hashrsig)
sum(duplicated(DEgenelist$GeneID) == TRUE) #351 duplicates
DEgenelist3 <- DEgenelist[!duplicated(DEgenelist$GeneID), ] #1886 DEGs intotal
#write.csv(DEgenelist3, "DEgenelist3_23012022.csv")

sessionInfo()







