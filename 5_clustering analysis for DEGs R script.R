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

#### clustering analysis to explore the DEGs ####

# load packages
library(SummarizedExperiment)
library("grid")
library("ComplexHeatmap")
library("circlize")
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(stringr)
library(reshape2)
library(dendextend)
library(ggpubr)

#zscore and rlog transformation for all genes all samples
plotDat <- rlog(ddsbatch) %>% 
  assay()
z.mat_rlog <- t(scale(t(plotDat), center=TRUE, scale=TRUE))
z.mat_rlog <- as.data.frame(z.mat_rlog) %>%  rownames_to_column("GeneID") 
plotDat <- as.data.frame(plotDat) %>%  rownames_to_column("GeneID") 
#fwrite(plotDat, "rlogdata_allgenesallsamples_30092022.csv")

#zscore and rlog transformation for DEGs
DEG <- read.csv("DEgenelist3_23012022.csv")
sigGenes <- DEG %>% pull("GeneID") #geneid for DEGs
plotDat <- rlog(ddsbatch)[sigGenes,] %>% 
  assay()
z.mat_rlog <- t(scale(t(plotDat), center=TRUE, scale=TRUE))
z.mat_rlog <- as.data.frame(z.mat_rlog) %>%  rownames_to_column("GeneID") 
#fwrite(z.mat_rlog, "z.matrlog_allsamples_20052022.csv")

### Step 1 heatmap hierarchical clustering of DEGs ###
hcDat <- hclust(dist(z.mat_rlog))
cutGroups <- cutree(hcDat, h=4)

ha1 = HeatmapAnnotation(df = colData(ddsbatch)[,c("Treatment", "TreatmentDuration")])

# add group annotation
anno<-data.frame(row.names=c("DC_05H_1","DC_05H_2","LT_05H_1","LT_05H_2","LT_05H_3",
                             "DC_1H_1","DC_1H_2","LT_1H_1","LT_1H_2","LT_1H_3",
                             "DC_2H_1","DC_2H_2","LT_2H_1","LT_2H_2","LT_2H_3",
                             "DC_4H_1","DC_4H_2","LT_4H_1","LT_4H_2","LT_4H_3"),
                 Treatment=c("DC","DC","LT","LT","LT",
                             "DC","DC","LT","LT","LT",
                             "DC","DC","LT","LT","LT",
                             "DC","DC","LT","LT","LT"),
                 TreatmentDuration=c("0.5","0.5","0.5","0.5","0.5",
                                     "1","1","1","1","1","2","2","2","2","2",
                                     "4","4","4","4","4"))

#for rows colouring by cluster
my_hclust_gene <- hclust(dist(z.mat_rlog), method = "complete")

as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)

#form cluster of 3
my_gene_row <- cutree(tree = as.dendrogram(my_hclust_gene), k = 3)
my_gene_row <- data.frame(cluster = ifelse(test = my_gene_row == 1, yes = "cluster 1",
                                           ifelse(test = my_gene_row == 2, yes = "cluster 2", no = "cluster 3")))

# colour palette
annotation_colors = list(
  Treatment = c(DC="dimgrey", LT="yellow"),
  TreatmentDuration = c("0.5"="bisque","1"="burlywood","2"="lightsalmon","4"="plum"),
  cluster = c("cluster 1"= "grey70", "cluster 2"="grey40","cluster 3"="grey30"))

#FOR export graph into ppt
p1 <- as.ggplot(~pheatmap::pheatmap(z.mat_rlog, cluster_rows=TRUE, cluster_cols=FALSE, #cellwidth = 65, cellheight = 17,
                                    # main = "Hierarchical clustering analysis of all DEGs", 
                                    color=colorRampPalette(c("blue3", "white", "red3"))(20),
                                    annotation_col = anno, annotation_colors = annotation_colors, 
                                    gaps_col = cumsum(c(5,5,5,5)),
                                    show_rownames = FALSE, annotation_row = my_gene_row,
                                    cutree_rows=3,
                                    legend_breaks = c(min(z.mat_rlog)+0.5,0, max(z.mat_rlog)-0.5), legend_labels = c("low","","high"),
                                    clustering_distance_rows = "correlation", #pearson correlation
                                    show_colnames = FALSE,
                                    fontsize=12))
p1

### Step 2 time series clustering analysis for the DEGs ###
#calculate differences between LT and DC at each time point, count value
d <- setDT(plotDat, keep.rownames = TRUE)[]
d <- melt(data = plotDat, id.vars = 'rn')
d$variable <- as.factor(d$variable)
levels(d$variable) <- c('DC_05H','DC_05H','DC_1H','DC_1H','DC_2H','DC_2H','DC_4H','DC_4H',
                        'LT_05H','LT_05H','LT_05H','LT_1H','LT_1H','LT_1H','LT_2H','LT_2H','LT_2H','LT_4H','LT_4H','LT_4H')
d$rn <- as.factor(d$rn)
d$variable <- as.character(d$variable)
x <- as.data.frame(str_split_fixed(d$variable, "_", 2))
d$group <- x$V1
d$duration <- x$V2
d$sd <- NULL
d$variable <-NULL

d$rn <- as.factor(d$rn)
d$duration <- as.factor(d$duration)
d$group <- as.factor(d$group)
d<- dcast(d, rn + duration ~ group, fun.aggregate=mean)
d$results <- d$LT - d$DC
dx <- d
dx$DC <- NULL
dx$LT <- NULL
dx<- dcast(dx, rn ~ duration, value='results')
colnames(dx)[1] <- "GeneID"
DEG$GeneID <- as.factor(DEG$GeneID)

dx <- as.data.table(dx)
DEG <- as.data.table(DEG)
dx <- dx[DEG, on = 'GeneID']

dx$GeneID <- as.character(dx$GeneID)
dx <- as.data.frame(dx)
rownames(dx) <- dx$GeneID
dx$GeneID <- NULL
colnames(dx)[1:4] <- c('30 mins', "1 hr" , "2 hrs", "4 hrs")
dx$X <- NULL
dx <- as.matrix(dx)
z.mat_heatmap <- t(scale(t(dx), center=TRUE, scale=TRUE))

#clustering analysis with TCseq
set.seed(120)

tca1 <- timeclust(z.mat_heatmap, algo = "km", k = 10, dist= "euclidean", standardize = TRUE)
p <- timeclustplot(tca1, categories = "treatment duration (h)", value = "rlog transformed count data", cols = 5)
p

datatca1 <- clustData(tca1)
dtca1 <- clustCluster(tca1)
dtca1 <- as.data.frame(dtca1)
#write.csv(dtca1,"dtca1_20220118.csv")
dtca1 %>% count(dtca1)

sessionInfo()
