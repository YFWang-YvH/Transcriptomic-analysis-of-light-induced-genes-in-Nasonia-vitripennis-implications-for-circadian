##Using StringTie with DESeq2 to analysis Nasonia RNAseq data
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
library(pheatmap)
library(RColorBrewer)
library(ggrepel) #to add text in the graph
library(ggfortify) #another type of PCA
library(glmpca) #generalized PCA

#### Preprocessing RNAseq count matrix ####

### Step 1 load gene count data and sample information ###

#load gene/transcript count matrix and labels
countData <- as.matrix(read.csv("gene_count_matrix_annotatedstringtie.csv", row.names = "gene_id"))
colData <- read.csv("sampleinfo.csv", row.names = 1)

#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

# check dimension of count matrix
dim(countData) #count matrix 3 15904 20
head(countData)
head(rownames(countData))
countData

#adjust colData to factors
colData$LaneID <- as.factor(colData$LaneID)
colData$Batch <- as.factor(colData$Batch)
colData$RNAextractionround <- as.factor(colData$RNAextractionround)
colData$Treatment <- as.factor(colData$Treatment)
colData$TreatmentDuration <- as.factor(colData$TreatmentDuration)
colData$Generationnumber <- as.factor(colData$Generationnumber)
str(colData)
levels(colData$Treatment)
levels(colData$TreatmentDuration)

#Create a DESeqDataSet from count matrix and labels
#can add different factors into the model
#such as generation information and rna extraction information, but model matrix not full rank so ignored
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData, design = ~ Batch +Treatment + TreatmentDuration + Treatment:TreatmentDuration)

### Step 2 Exploratory analysis and visualization ###

## A. pre-filtering the dataset ##
#we decided to not set threshold to filter out low expression genes as we suspect that light effects may be small on gene expression
nrow(dds) #new data 21289 15904
ddsClean <- dds[which(mcols(dds)$betaConv),]
dds <- ddsClean

## B. Quality assessment - important to assess the quality of our data before moving on to do the actual differential expression analysis
#data transformations and visualization
#Differential expression calculations with DESeq2 uses raw read counts as input, 
#but for visualization purposes we use transformed counts.

hist(countData[,1])
summary(dds)
summary(countData)

# few outliers affect distribution visualization
p1 <- boxplot(countData, main='Raw counts', las=2)
p1

# Raw counts mean expression Vs standard Deviation (SD)
p2 <- plot(rowMeans(countData), rowSds(countData), 
           main='Raw counts: sd vs mean', 
           xlim=c(0,10000),
           ylim=c(0,5000))  #zoom in to get rid of outliers
abline(lm(rowSds(countData) ~ rowMeans(countData)), col="green") #add regression line

## C. data transformation and different types of visualization ##

#to avoid problems posed by raw counts, they can be transformed, several transformation methods exist to limit the dependence of variance 
#on mean gene expression

#use rlog transformation
rld <- rlog(dds, blind = FALSE)

head(rld, 3)
head(vsd,3)
colData(vsd)

#rlog distribution
str(colData)
levels(colData$Sample)
rlogcounts <- as.matrix(assay(rld))
str(rlogcounts)
samplecol <- match(colData$Sample, c("DC_05H_1","DC_05H_2","DC_1H_1","DC_1H_2","DC_2H_1","DC_2H_2","DC_4H_1","DC_4H_2" , "LT_05H_1" ,"LT_05H_2" ,"LT_05H_3", "LT_1H_1" , "LT_1H_2" , "LT_1H_3", "LT_2H_1",  "LT_2H_2" , "LT_2H_3" , "LT_4H_1" , "LT_4H_2",  "LT_4H_3" )) + 1 #+1 to avoid color 1 which is black
boxplot(rlogcounts,
        las=2,
        xlab="",
        ylab="rlog counts",
        col=samplecol,
        main="rlog counts distribution")
abline(h=median(rlogcounts), col="blue") #median logcounts h=horizontal line

#rlog counts sd vs mean
rlogmeans <- rowMeans(rlogcounts)
rlogsds <- rowSds(rlogcounts)
plot(x=rlogmeans, y=rlogsds,
     main="rlog Counts: sd vs mean")

abline(lm(rlogsds ~ rlogmeans), col="green")

meanSdPlot(assay(rld))

# sample distances
sampleDists <- dist(t(assay(rld)))
sampleDists

#visualize the distances in a heatmap figure

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Sample, sep = "_")
colnames(sampleDistMatrix) <- paste(rld$Sample, sep = "_")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         col = colors,
         display_numbers = T,
         number_format = '%.0f',
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists)

#data quality assessment by sample clustering and visualization
#heatmap of the count matrix, top 20
dds <- DESeq(dds) 

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Treatment","TreatmentDuration")])

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# PCA plot-another way to visualize sample-to-sample distances
#PCA with rld data
pcaData <- plotPCA(rld, intgroup = c( "Treatment", "Batch"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = Treatment, shape = Batch)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with RLD data")+
  geom_text_repel(label = pcaData$name, box.padding = 0.8)

#run PCA with prcomp
pcDAT <- prcomp(t(rlogcounts)) #transpose the matrix
#plot PCA
autoplot(pcDAT,
         data=colData,
         colour="Treatment",
         shape="TreatmentDuration",
         size=5) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Sample), box.padding = 0.8)

#PCA plot using generalized PCA
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$Treatment <- dds$Treatment
gpca.dat$TreatmentDuration <- dds$TreatmentDuration
gpca.dat$Sample <- dds$Sample

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Treatment, shape = TreatmentDuration)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")+
  geom_text_repel(label = gpca.dat$Sample, box.padding = 0.8)

# MDS plot
#similar to the PCA plot, using the multidimensional scaling funtion
#plot a matrix of distances

#mds with RLD data
mds <- as.data.frame(colData(rld)) %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Treatment, shape = TreatmentDuration)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with RLD data")+
  geom_text_repel(label = gpca.dat$Sample, box.padding = 0.8)

### Step 3 Compare models for next step analysis ###

dds1 <- DESeq(dds, test="LRT", reduced = ~ Treatment+TreatmentDuration+ TreatmentDuration:Treatment)
resTC <- results(dds1)
head(resTC[order(resTC$padj),], 4)
sum(resTC$padj < 0.05, na.rm=TRUE)
#332genes affected in the new data

dds2 <- DESeq(dds, test="LRT", reduced = ~ Treatment:TreatmentDuration)
resTC <- results(dds2)
head(resTC[order(resTC$padj),], 4)
sum(resTC$padj < 0.05, na.rm=TRUE)

dds3 <- DESeq(dds, test="LRT", reduced = ~ Batch)
resTC <- results(dds3)
head(resTC[order(resTC$padj),], 4)
sum(resTC$padj < 0.05, na.rm=TRUE)

# we compared different models and decided to use the original model of dds that fits our question better
sessioninfo()
