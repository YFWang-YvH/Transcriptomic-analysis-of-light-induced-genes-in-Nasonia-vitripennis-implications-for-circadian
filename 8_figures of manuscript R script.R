##Using StringTie with DESeq2 for Nasonia RNAseq data
rm(list=ls()) #clear all objects
gc()

#dataset: 
#Nasonia RNAseq data
#experimental light treatment and dark control
#4 different time point: 0.5h 1h 2h 4h
#8 sample groups
#3 replicates per experimental sample group, 2 replicates per control sample group
#20 samples in total

#### Figures plotting for the RNAseq manuscripts ####

#load packages 
library("vctrs")
library("ellipsis")
library(DESeq2)
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library("VennDiagram")
library(ggplotify)
library(SummarizedExperiment)
library("grid")
library("ComplexHeatmap")
library("circlize")
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(dendextend)
library(dplyr)
library(cowplot)
library(ggpubr)
library(clusterProfiler)
library(R.utils)
library(forcats)
library(enrichplot)
library(ggplotify)
library("TCseq")
library(universalmotif)
library(ggseqlogo)

### Figure 2 in the manuscript ### 
#load differential expression analysis results from 4_differential expression analysis R script

LTvsDC0.5H <- dres0.5Hashrsig$GeneID
LTvsDC1H <- dres1Hashrsig$GeneID
LTvsDC2H <- dres2Hashrsig$GeneID
LTvsDC4H <- dres4Hashrsig$GeneID

venn_list <- list("30 mins" = LTvsDC0.5H, "4 hrs"=LTvsDC4H, "1 hr" = LTvsDC1H, "2 hrs" = LTvsDC2H)

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
# export graph into ppt

p2 <- as.ggplot(~display_venn(
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

### Figure 3 in the manuscript ###

#Fig. 3a heatmap for all samples hieracial clustering
z.mat_rlog <- as.matrix(read.csv("z.matrlog_allsamples_20052022.csv", row.names = "GeneID"))
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
p3a <- as.ggplot(~pheatmap::pheatmap(z.mat_rlog, cluster_rows=TRUE, cluster_cols=FALSE, #cellwidth = 65, cellheight = 17,
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
p3a

#Fig. 3b GO plot FOR EACH TIME POINT up or down regulated
GOsig <- read.csv("allRes_timepoint_updownGO.csv")
GOsig$Percentage <- GOsig$Significant/GOsig$Annotated
GOsig <- GOsig %>% filter(elimFisher <=0.01)
GOsig$category <- as.factor(GOsig$category)
#fwrite(GOsig, "GOsig_timepointupdown.csv")

G0.5up <- subset(GOsig, number ==1) 
G0.5down <- subset(GOsig, number==2)
G01up <- subset(GOsig, number==3)
G01down <- subset(GOsig, number==4)
G02up <- subset(GOsig, number==5)
G02down <- subset(GOsig, number==6)
G04up <- subset(GOsig, number==7)
G04down <- subset(GOsig, number==8)

G0.5up$Term <- factor(G0.5up$Term, levels = G0.5up$Term)
G0.5down$Term <- factor(G0.5down$Term, levels = G0.5down$Term)
G01up$Term <- factor(G01up$Term, levels = G01up$Term)
G01down$Term <- factor(G01down$Term, levels = G01down$Term)
G02up$Term <- factor(G02up$Term, levels = G02up$Term)
G02down$Term <- factor(G02down$Term, levels = G02down$Term)
G04up$Term <- factor(G04up$Term, levels = G04up$Term)
G04down$Term <- factor(G04down$Term, levels = G04down$Term)

pGO0.5up <- ggplot(G0.5up, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="bisque") +
  #scale_y_continuous(sec.axis = sec_axis(~ ., name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes")
pGO0.5up

pGO0.5down <- ggplot(G0.5down, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="bisque") +
  #scale_y_continuous(sec.axis = sec_axis(~ ., name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes")
pGO0.5down
pGO1up <- ggplot(G01up, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity",fill="burlywood") +
  # scale_y_continuous(sec.axis = sec_axis(~ ., name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes")
pGO1down <- ggplot(G01down, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity",fill="burlywood") +
  # scale_y_continuous(sec.axis = sec_axis(~ ., name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes")
pGO2up <- ggplot(G02up, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="lightsalmon") +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes")
pGO2down <- ggplot(G02down, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="lightsalmon") +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes")
pGO4up <- ggplot(G04up, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="plum") +
  #scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes")
pGO4down <- ggplot(G04down, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="plum") +
  #scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes")

p <- ggarrange(pGO0.5up,
               pGO1up+rremove('y.title'),
               pGO2up+rremove('y.title'),
               pGO4up+rremove('y.title'),
               pGO0.5down,
               pGO1down+rremove('y.title'),
               pGO2down+rremove('y.title'),
               pGO4down+rremove('y.title'),
               ncol=4,
               nrow=2,
               align = "hv",
               # labels=c("(b)","(c)","(d)","(e)"),
               common.legend=T,
               legend = 'right')

p

pp <- ggarrange(p1,p,
                nrow=2,ncol=1)
pp

### Figure 4 in the manuscript ###
#kegg for each time point
DEG <- read.csv("DEgenelist3_23012022.csv")

d0.5genelist <- DEG %>% dplyr::filter(time=="0.5") %>% drop_na(FDR) %>% filter(FDR <0.05) #277
d1genelist <- DEG %>% dplyr::filter(time=="1") %>% drop_na(FDR) %>% filter(FDR <0.05) #381
d2genelist <- DEG %>% dplyr::filter(time=="2") %>% drop_na(FDR) %>% filter(FDR <0.05) #1432
d4genelist <- DEG %>% dplyr::filter(time=="4") %>% drop_na(FDR) %>% filter(FDR <0.05) #147

d0.5G <- d0.5genelist[!d0.5genelist$Genename=="gene-",]
d0.5G <- d0.5G[, c(11,20)]
d0.5G<- d0.5G[!duplicated(d0.5G),] 
d0.5G <- filter(d0.5G)$NCBI.GeneID %>% unique %>% 
  as.character()#270

d1G <- d1genelist[!d1genelist$Genename=="gene-",]
d1G <- d1G[, c(11,20)]
d1G<- d1G[!duplicated(d1G),] 
d1G <- filter(d1G)$NCBI.GeneID %>% unique %>% 
  as.character()#363

d2G <- d2genelist[!d2genelist$Genename=="gene-",]
d2G <- d2G[, c(11,20)]
d2G<- d2G[!duplicated(d2G),]  
d2G <- filter(d2G)$NCBI.GeneID %>% unique %>% 
  as.character()#1386

d4G <- d4genelist[!d4genelist$Genename=="gene-",]
d4G <- d4G[, c(11,20)]
d4G<- d4G[!duplicated(d4G),] 
d4G <- filter(d4G)$NCBI.GeneID %>% unique %>% 
  as.character()#145

#KEGG pathway over-representaion analysis
#do this to successfully download from kegg database
R.utils::setOption("clusterProfiler.download.method","auto")

search_kegg_organism('nvi', by='kegg_code')

kk0.5 <- enrichKEGG(gene=d0.5G, organism="nvi", 
                    pvalueCutoff =0.05, 
                    pAdjustMethod="BH",
                    qvalueCutoff = 0.05)
kk0.5
head(kk0.5)
kk0.5 <- as.data.frame(kk0.5)
kk0.5$timepoint <- 0.5

#kegg module over-representaion analysis
mkk0.5 <- enrichMKEGG(gene=d0.5G, organism="nvi", 
                      pvalueCutoff =0.05,
                      qvalueCutoff = 0.05)
head(mkk0.5)  
mkk0.5 <- as.data.frame(mkk0.5)
mkk0.5$timepoint <- 0.5

kk_final <- rbind(kk0.5,kk1,kk2,kk4)
mkk_final <- rbind(mkk0.5, mkk1,mkk2,mkk4)

#biological theme comparison
DEGSaple <- list("0.5"=d0.5G, "1"=d1G,"2"=d2G,"4"=d4G)
xx <- compareCluster(DEGSaple, fun="enrichKEGG",
                     organism="nvi", pvalueCutoff=0.05)
xx <- pairwise_termsim(xx)                     

p4 <- as.ggplot(~dotplot(xx))
p4.1 <-  cnetplot(xx)

### Figure 5 in the manuscript ###
#Fig 5a clustering analysis with TCseq
set.seed(120)

tca1 <- timeclust(z.mat_heatmap, algo = "km", k = 10, dist= "euclidean", standardize = TRUE)
p <- timeclustplot(tca1, categories = "treatment duration (h)", value = "rlog transformed count data", cols = 5)

p_1 <- print(p[[8]])
p_1<-p_1+ggtitle('Cluster 1(N=76)')+ylim(-1.5,1.5)
p_2 <- print(p[[7]])
p_2<-p_2+ggtitle('Cluster 2(N=116)')+ylim(-1.5,1.5)
p_3 <- print(p[[1]])
p_3<-p_3+ggtitle('Cluster 3(N=293)')+ylim(-1.5,1.5)
p_4 <- print(p[[3]])
p_4<-p_4+ggtitle('Cluster 4(N=177)')+ylim(-1.5,1.5)
p_5 <- print(p[[6]])
p_5<-p_5+ggtitle('Cluster 5(N=209)')+ylim(-1.5,1.5)
p_6 <- print(p[[9]])
p_6<-p_6+ggtitle('Cluster 6(N=184)')+ylim(-1.5,1.5)
p_7 <- print(p[[5]])
p_7<-p_7+ggtitle('Cluster 7(N=134)')+ylim(-1.5,1.5)
p_8 <- print(p[[4]])
p_8<-p_8+ggtitle('Cluster 8(N=292)')+ylim(-1.5,1.5)
p_9 <- print(p[[10]])
p_9<-p_9+ggtitle('Cluster 9(N=223)')+ylim(-1.5,1.5)
p_10 <- print(p[[2]])
p_10<-p_10+ggtitle('Cluster 10(N=182)')+ylim(-1.5,1.5)
p <- ggarrange(p_1+ rremove('x.title')+rremove("x.text")+rremove('y.title'),p_2+ rremove('y.title')+rremove('x.title')+rremove('y.text')+rremove("x.text"),
               p_3+ rremove('y.title')+rremove('x.title')+rremove('y.text')+rremove("x.text"),p_4+ rremove('y.title')+rremove('x.title')+rremove('y.text')+rremove("x.text"),
               p_5 + rremove('y.title')+rremove('x.title')+rremove('y.text')+rremove("x.text"),
               p_6+ rremove('x.title')+rremove('y.title'),p_7+ rremove('y.title')+rremove('x.title')+rremove('y.text'),
               p_8+ rremove('y.title')+rremove('x.title')+rremove('y.text'),p_9+ rremove('y.title')+rremove('x.title')+rremove('y.text'),
               p10 + rremove('y.title')+rremove('x.title')+rremove('y.text'),
               ncol=5,
               nrow=2,
               align = "hv",
               labels=c("A","B","C","D","E","F","G","H","I","J"),
               common.legend=T,
               legend = 'right')
p5a <- annotate_figure(p, left = text_grob("rlog transformed count data",rot=90, size = 18, face = "bold"),
                        bottom=text_grob("treatment duration (h)", size = 18, face = "bold"))

p5a

#Fig 5b GO plot for each cluster
GOsig <- read.csv("allRessig_clusterGO.csv")
GOsig$category <- as.factor(GOsig$category)

for (i in 1:10) {
  Gi <- subset(GOsig,  Module == i)
  Gi$Term <- factor(Gi$Term, levels = Gi$Term)
  Gi$Percentage <- Gi$Significant/Gi$Annotated
  assign(paste("Gi", i, sep=""), Gi)
}

pGO2 <- ggplot(Gi2, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="grey") +
  #scale_y_continuous(sec.axis = sec_axis(~ ., name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes", title = "Cluster 2 (N=122)")+ylim(0,45)
pGO3 <- ggplot(Gi3, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity",fill="grey") +
  # scale_y_continuous(sec.axis = sec_axis(~ ., name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes",title = "Cluster 3 (N=91)")+ylim(0,45)
pGO4 <- ggplot(Gi4, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="grey") +
  # scale_y_continuous(sec.axis = sec_axis(~ . , name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes",title = "Cluster 4 (N=264)")+ylim(0,45)
pGO5 <- ggplot(Gi5, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="grey") +
  #scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes",title = "Cluster 5 (N=264)")+ylim(0,45)
pGO6 <- ggplot(Gi6, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="grey") +
  #scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes", ,title = "Cluster 6 (N=184)")+ylim(0,45)
pGO7 <- ggplot(Gi7, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="grey") +
  #scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes", ,title = "Cluster 7 (N=134)")+ylim(0,45)
pGO8 <- ggplot(Gi8, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="grey") +
  #scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes", ,title = "Cluster 8 (N=292)")+ylim(0,45)
pGO9 <- ggplot(Gi9, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="grey") +
  #scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes", ,title = "Cluster 9 (N=223)")+ylim(0,45)
pGO10 <- ggplot(Gi10, aes(x=Term, y=Percentage*100)) +
  geom_bar(stat="identity", fill="grey") +
  #scale_y_continuous(sec.axis = sec_axis(~ . * 10, name = "Significant"))+
  theme_light()+ background_grid(major = 'none', minor = "none") +
  theme(axis.text.x = element_text(angle = 80,hjust=1,vjust=1),
        axis.title.x=element_blank())+#hjust to leave some space between the labels and the axis #vjust to center them in this case
  labs(y="percentage of genes", ,title = "Cluster 10 (N=182)")+ylim(0,45)

p <- ggarrange(pGO2+rremove('y.title'), pGO2+rremove('y.title')+rremove("y.text"),
               pGO3+rremove('y.title')+rremove("y.text"),
               pGO4+rremove('y.title')+rremove("y.text"),
               pGO5+rremove('y.title')+rremove("y.text"),
               pGO6+rremove('y.title'),
               pGO7+rremove('y.title')+rremove("y.text"),
               pGO8+rremove('y.title')+rremove("y.text"),
               pGO9+rremove('y.title')+rremove("y.text"),
               pGO10+rremove('y.title')+rremove("y.text"),
               ncol=5,
               nrow=2,
               align = "hv",
               labels=c("(a)","(b)","(c)","(d)","(e)",
                        "(f)","(g)","(h)","(i)","(j)"),
               common.legend=T,
               legend = 'right')

p
p5b <- annotate_figure(p, left = text_grob("percentage of genes",rot=90, size = 18, face = "bold"))

### Figure 6 in the manuscript ###
#circadian related interesting DEG for plotting
circadian <- fread("circadianGenes_manuscript.csv")
circadian <- melt(circadian, id.vars = c('GeneID',"GeneName","Category","DEGs","Value"),variable.name = 'samplename') #transpose the dataset for further plotting
circadian1 <- as.data.frame(str_split_fixed(circadian$samplename, "_", 3))
circadian2 <- cbind(circadian,circadian1)
#fwrite(circadian2, "circadianGenes_manuscript_R.csv")

circadian <- fread("circadianGenes_manuscript_R.csv")
circadian <- subset(circadian, circadian$Value =="rlog")
circadian1 <- subset(circadian, circadian$DEGs=="1")
circadian0 <- subset(circadian, circadian$DEGs=="0")
circadianmean1 <- setDT(circadian1)[,.(mean=mean(value), sd = sd(value)), by = c('Treatment', 'Treatmentduration',"GeneName")]
circadianmean0<- setDT(circadian0)[,.(mean=mean(value), sd = sd(value)), by = c('Treatment', 'Treatmentduration',"GeneName")]


circadian1$GeneName = factor(circadian1$GeneName, 
                                 levels = c("OpBlue", "cryptochrome 2", "neuronal PAS domain-containing protein 2", 
                                           "beta-TrCP","tra","juvenile hormone acid O-methyltransferase",
                                           "circadian clock-controlled protein-like"))
circadian1$Treatment = factor(circadian1$Treatment, levels=c("LT","DC"))

circadianmean1$GeneName = factor(circadianmean1$GeneName, 
                             levels = c("OpBlue", "cryptochrome 2", "neuronal PAS domain-containing protein 2", 
                                        "beta-TrCP","tra","juvenile hormone acid O-methyltransferase",
                                        "circadian clock-controlled protein-like"))
circadianmean1$Treatment = factor(circadianmean1$Treatment, levels=c("LT","DC"))

p6 <- ggplot(circadian1, mapping=aes(x=Treatmentduration, y=value, colour=Treatment,linetype=Treatment))+
  #geom_errorbar(aes(ymin=meanzscore-sd, ymax=meanzscore+sd), width=0) +
  geom_line(circadianmean1, mapping=aes(x=Treatmentduration, y=mean, colour=Treatment,linetype=Treatment), size=1.1) + 
  geom_point(size=2)+
  facet_wrap(~GeneName, ncol = 3,scales = "free")+
  theme_classic()+
  scale_color_brewer(palette="Paired")+
  scale_x_continuous(breaks = c(0.5,1,2,4))+
  ylab("rlog transformed gene expression") +xlab("treatment duration (hr)")+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold") )+
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) +theme(legend.position = "none")+#change legend text fo+
  theme(strip.background = element_rect(colour="light grey", fill="light grey", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold"))#change the apperance of the rectangle around facet label

p6

#supplementary S1
circadian0$GeneName = factor(circadian0$GeneName, levels = c("OpUV","OpLW","rhodopsin","rhodopsin-like",
                                                                     "period","cycle","clock1","clock2",
                                                                     "timeless homolog","circadian clock-controlled protein"))
circadian0$Treatment = factor(circadian0$Treatment, levels=c("LT","DC"))

circadianmean0$GeneName = factor(circadianmean0$GeneName, levels = c("OpUV","OpLW","rhodopsin","rhodopsin-like",
                                                                     "period","cycle","clock1","clock2",
                                                                     "timeless homolog","circadian clock-controlled protein"))
circadianmean0$Treatment = factor(circadianmean0$Treatment, levels=c("LT","DC"))

pS1 <- ggplot(circadian0, mapping=aes(x=Treatmentduration, y=value, colour=Treatment,linetype=Treatment))+
  #geom_errorbar(aes(ymin=meanzscore-sd, ymax=meanzscore+sd), width=0) +
  geom_line(circadianmean0, mapping=aes(x=Treatmentduration, y=mean, colour=Treatment,linetype=Treatment), size=1.1) + 
  geom_point(size=2)+
  facet_wrap(~GeneName, ncol = 3,scales = "free")+
  theme_classic()+
  scale_color_brewer(palette="Paired")+
  scale_x_continuous(breaks = c(0.5,1,2,4))+
  ylab("rlog transformed gene expression") +xlab("treatment duration (hr)")+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold") )+
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) +theme(legend.position = "none")+#change legend text fo
  theme(strip.background = element_rect(colour="light grey", fill="light grey", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold"))#change the apperance of the rectangle around facet label

pS1


### Figure 7 in the manuscript ###
#interesting glutmate related DEGs
glutmate <- fread("glutamate for manuscript.csv")
glutmate <- melt(glutmate, id.vars = c('GeneID',"GeneName","CommonName","Value"),variable.name = 'samplename') #transpose the dataset for further plotting
glutmate1 <- as.data.frame(str_split_fixed(glutmate$samplename, "_", 3))
glutmate2 <- cbind(glutmate,glutmate1)
#fwrite(glutmate2, "glutmateGenes_manuscript_R.csv")

glutmate <- fread("glutmateGenes_manuscript_R.csv")
glutmate <- subset(glutmate, glutmate$Value == "rlog")
glutmatemean <- setDT(glutmate)[,.(meanzscore=mean(value), sd = sd(value)), by = c('Treatment', 'Treatmentduration',"GeneName","CommonName")]

glutmate$Treatment = factor(glutmate$Treatment, levels=c("LT","DC"))
glutmatemean$Treatment = factor(glutmatemean$Treatment, levels=c("LT","DC"))

p7 <- ggplot(glutmate, mapping=aes(x=Treatmentduration, y=value, colour=Treatment,linetype=Treatment))+
  #geom_errorbar(aes(ymin=meanzscore-sd, ymax=meanzscore+sd), width=0) +
  geom_line(glutmatemean, mapping=aes(x=Treatmentduration, y=meanzscore, colour=Treatment,linetype=Treatment), size=1.1) + 
  geom_point(size=2)+
  facet_wrap(~CommonName, ncol = 3,scales = "free")+
  theme_classic()+
  scale_color_brewer(palette="Paired")+
  scale_x_continuous(breaks = c(0.5,1,2,4))+
  ylab("rlog transformed gene expression") +xlab("treatment duration (hr)")+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold") )+
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) +theme(legend.position = "none")+#change legend text fo+
  theme(strip.background = element_rect(colour="light grey", fill="light grey", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold"))#change the apperance of the rectangle around facet label

p7

### Figure 8 in the manuscript ###
#interesting DEG for plotting of creb and ap1 pathways
TF <- fread("transcription factor for manuscript.csv")
TF <- melt(TF, id.vars = c('GeneID',"GeneName","CommonName","Value"),variable.name = 'samplename') #transpose the dataset for further plotting
TF1 <- as.data.frame(str_split_fixed(TF$samplename, "_", 3))
TF2 <- cbind(TF,TF1)
#fwrite(TF2, "TF_manuscript_R.csv")

TF <- fread("TF_manuscript_R.csv")
TF <- subset(TF, TF$Value =="rlog")
TF$Treatment = factor(TF$Treatment, levels=c("LT","DC"))


CREB <- subset(TF, TF$CommonName=="CREB" |
                  TF$CommonName=="CREBBP1" |
                  TF$CommonName=="CREBBP2")
CREB_mean <- setDT(CREB)[,.(meanzscore=mean(value), sd = sd(value)), by = c('Treatment', 'Treatmentduration',"GeneName","CommonName")]

p8_1 <- ggplot(CREB, mapping=aes(x=Treatmentduration, y=value, colour=Treatment,linetype=Treatment))+
  #geom_errorbar(aes(ymin=meanzscore-sd, ymax=meanzscore+sd), width=0) +
  geom_line(CREB_mean, mapping=aes(x=Treatmentduration, y=meanzscore, colour=Treatment,linetype=Treatment), size=1.1) +
  geom_point(size=2)+
  theme_classic()+ facet_wrap(~CommonName, ncol = 3,scales = "free")+
  scale_color_brewer(palette="Paired")+
  scale_x_continuous(breaks = c(0.5,1,2,4))+
  ylab("rlog transformed gene expression") +xlab("treatment duration (hr)")+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold") )+
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))+theme(legend.position = "none")+
  theme(strip.background = element_rect(colour="light grey", fill="light grey", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold.italic"))#change the apperance of the rectangle around facet label
#import matrix and plot motif
CREB_motifmatrix <- read_homer(file="cre_c1.motif")
CREB1_motifmatrix<- read_homer(file="cre_c1c5.motif")
CREB2_motifmatrix<- read_homer(file="CREB_c2.motif")

CREB_motif<- convert_type(CREB_motifmatrix, "PPM")
CREB1_motif<- convert_type(CREB1_motifmatrix, "PPM")
CREB2_motif<- convert_type(CREB2_motifmatrix, "PPM")
## Only need the matrix itself
CREB_motif <- CREB_motif["motif"]
CREB1_motif <- CREB1_motif["motif"]
CREB2_motif <- CREB2_motif["motif"]
## ggseqlogo:
p8_2<- ggplot()+
  geom_logo(data=CREB_motif, method="probability",seq_type = "DNA",stack_width = 0.95,
            rev_stack_order = F)+theme_classic()+
  scale_x_continuous(breaks=seq(1,12,1))+
  scale_y_continuous(breaks=seq(0,1,0.25))+
  geom_segment(aes(x=1,xend=12,y=-Inf,yend=-Inf))+
  geom_segment(aes(y=0,yend=1,x=-Inf,xend=-Inf))+
  ylab("probability") +xlab("position")+
  ggtitle("c-Jun-CRE(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer")+
  theme(axis.line=element_blank(),
        axis.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold"))+
  theme(plot.title = element_text(size=12, color="black",
                                  face="bold"))#change the apperance of the rectangle around facet label
p8_3<- ggplot()+
  geom_logo(data=CREB1_motif, method="probability",seq_type = "DNA",stack_width = 0.95,
            rev_stack_order = F)+theme_classic()+
  scale_x_continuous(breaks=seq(1,12,1))+
  scale_y_continuous(breaks=seq(0,1,0.25))+
  geom_segment(aes(x=1,xend=12,y=-Inf,yend=-Inf))+
  geom_segment(aes(y=0,yend=1,x=-Inf,xend=-Inf))+
  ylab("probability") +xlab("position")+
  ggtitle("CRE(bZIP)/Promoter/Homer")+
  theme(axis.line=element_blank(),
        axis.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold"))+
  theme(plot.title = element_text(size=12, color="black",
                                  face="bold"))#change the apperance of the rectangle around facet label
p8_4<- ggplot()+
  geom_logo(data=CREB2_motif, method="probability",seq_type = "DNA",stack_width = 0.95,
            rev_stack_order = F)+theme_classic()+
  scale_x_continuous(breaks=seq(1,12,1))+
  scale_y_continuous(breaks=seq(0,1,0.25))+
  geom_segment(aes(x=1,xend=12,y=-Inf,yend=-Inf))+
  geom_segment(aes(y=0,yend=1,x=-Inf,xend=-Inf))+
  ylab("probability") +xlab("position")+
  ggtitle("CREB5(bZIP)/LNCaP-CREB5.V5-ChIP-Seq(GSE137775)/Homer")+
  theme(axis.line=element_blank(),
        axis.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold"))+
  theme(plot.title = element_text(size=12, color="black",
                                  face="bold"))#change the apperance of the rectangle around facet label

####################
AP1 <- subset(TF, TF$CommonName=="AP1B1" |
                 TF$CommonName=="AP1M1" )
AP1_mean <- setDT(AP1)[,.(meanzscore=mean(value), sd = sd(value)), by = c('Treatment', 'Treatmentduration',"GeneName","CommonName")]

p8_5 <- ggplot(AP1, mapping=aes(x=Treatmentduration, y=value, colour=Treatment,linetype=Treatment))+
  #geom_errorbar(aes(ymin=meanzscore-sd, ymax=meanzscore+sd), width=0) +
  geom_line(AP1_mean, mapping=aes(x=Treatmentduration, y=meanzscore, colour=Treatment,linetype=Treatment), size=1.1) + 
  geom_point(size=2)+
  theme_classic()+ facet_wrap(~CommonName, ncol = 3,scales = "free")+
  scale_color_brewer(palette="Paired")+
  scale_x_continuous(breaks = c(0.5,1,2,4))+
  ylab("rlog transformed gene expression") +xlab("treatment duration (hr)")+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold") )+
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))+theme(legend.position = "none")+
  theme(strip.background = element_rect(colour="light grey", fill="light grey", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold.italic"))#change the apperance of the rectangle around facet label

#import matrix and plot motif
AP1_motifmatrix <- read_homer(file="ap-1_c2c3.motif")
AP2_motifmatrix<- read_homer(file="ap1_c4c9.motif")

AP1_motif<- convert_type(AP1_motifmatrix, "PPM")
AP2_motif<- convert_type(AP2_motifmatrix, "PPM")
## Only need the matrix itself
AP1_motif <- AP1_motif["motif"]
AP2_motif <- AP2_motif["motif"]
## ggseqlogo:
p8_6<- ggplot()+
  geom_logo(data=AP1_motif, method="probability",seq_type = "DNA",stack_width = 0.95,
            rev_stack_order = F)+theme_classic()+
  scale_x_continuous(breaks=seq(1,10,1))+
  scale_y_continuous(breaks=seq(0,1,0.25))+
  geom_segment(aes(x=1,xend=10,y=-Inf,yend=-Inf))+
  geom_segment(aes(y=0,yend=1,x=-Inf,xend=-Inf))+
  ylab("probability") +xlab("position")+
  ggtitle("AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer")+
  theme(axis.line=element_blank(),
        axis.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold"))+
  theme(plot.title = element_text(size=12, color="black",
                                  face="bold"))#change the apperance of the rectangle around facet label
p8_7<- ggplot()+
  geom_logo(data=AP2_motif, method="probability",seq_type = "DNA",stack_width = 0.95,
            rev_stack_order = F)+theme_classic()+
  scale_x_continuous(breaks=seq(1,20,1))+
  scale_y_continuous(breaks=seq(0,1,0.25))+
  geom_segment(aes(x=1,xend=20,y=-Inf,yend=-Inf))+
  geom_segment(aes(y=0,yend=1,x=-Inf,xend=-Inf))+
  ylab("probability") +xlab("position")+
  ggtitle("NFAT:AP1(RHD,bZIP)/Jurkat-NFATC1-ChIP-Seq(Jolma_et_al.)/Homer")+
  theme(axis.line=element_blank(),
        axis.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold"))+
  theme(plot.title = element_text(size=12, color="black",
                                  face="bold"))#change the apperance of the rectangle around facet label

pp1 <- ggarrange(p8_1, p8_5,
                 ncol=2,
                 nrow=1,
                 align = "h",
                 labels=c("(a)","(e)"))
pp1


pp2 <- ggarrange(p8_2+rremove("x.title"), p8_6+rremove("x.title"),
                 p8_3+rremove("x.title"), p8_7,
                 p8_4,
                 ncol=2,
                 nrow=3,
                 align = "hv",
                 labels=c("(b)","(f)","(c)","(g)","(d)"))
pp2


p8 <- ggarrange(pp1,pp2, ncol=1,nrow = 2,
                  heights=c(1/4,3/4))

p8

### Figure 9 in the manuscript ###
#plotting tbx20 and znf467 trabscruotuib factir

#tbx20 transcription factor
TBX20 <- subset(TF, TF$CommonName=="TBX20")
TBX20_mean <- setDT(TBX20)[,.(meanzscore=mean(value), sd = sd(value)), by = c('Treatment', 'Treatmentduration',"GeneName","CommonName")]

p9_1 <- ggplot(TBX20, mapping=aes(x=Treatmentduration, y=value, colour=Treatment,linetype=Treatment))+
  #geom_errorbar(aes(ymin=meanzscore-sd, ymax=meanzscore+sd), width=0) +
  geom_line(TBX20_mean, mapping=aes(x=Treatmentduration, y=meanzscore, colour=Treatment,linetype=Treatment), size=1.1) + 
  geom_point(size=2)+
  theme_classic()+
  facet_wrap(~CommonName, ncol = 3,scales = "free")+
  scale_color_brewer(palette="Paired")+
  scale_x_continuous(breaks = c(0.5,1,2,4))+
  ylab("rlog transformed gene expression") +xlab("treatment duration (hr)")+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold") )+
  guides(color = guide_legend(override.aes = list(size = 1))) + theme(legend.position = "none")+#change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))+#change legend text fo
  theme(strip.background = element_rect(colour="light grey", fill="light grey", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold.italic"))#change the apperance of the rectangle around facet label


#import matrix and plot motif
tbx20_motifmatrix <- read_homer(file="TBX20.motif")

tbx20_motif <- convert_type(tbx20_motifmatrix, "PPM")
## Only need the matrix itself
tbx20_motif <- tbx20_motif["motif"]
## ggseqlogo:
p9_2 <- ggplot()+
  geom_logo(data=tbx20_motif, method="probability",seq_type = "DNA",stack_width = 0.95,
            rev_stack_order = F)+theme_classic()+
  scale_x_continuous(breaks=seq(1,12,1))+
  scale_y_continuous(breaks=seq(0,1,0.25))+
  geom_segment(aes(x=1,xend=12,y=-Inf,yend=-Inf))+
  geom_segment(aes(y=0,yend=1,x=-Inf,xend=-Inf))+
  ylab("probability") +xlab("position")+
  ggtitle("Tbx20(T-box)/Heart-Tbx20-ChIP-Seq(GSE29636)/Homer")+
  theme(axis.line=element_blank(),
        axis.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold"))+
  theme(plot.title = element_text(size=12, color="black",
                                  face="bold"))#change the apperance of the rectangle around facet label


pp3 <- ggarrange(p9_1, p9_2,
                 ncol=1,
                 nrow=2,
                 align = "v",
                 labels=c("(a)","(b)"))
pp3

####################
ZNF467 <- subset(TF, TF$CommonName=="ZNF467")
ZNF467_mean <- setDT(ZNF467)[,.(meanzscore=mean(value), sd = sd(value)), by = c('Treatment', 'Treatmentduration',"GeneName","CommonName")]

p9_3 <- ggplot(ZNF467, mapping=aes(x=Treatmentduration, y=value, colour=Treatment,linetype=Treatment))+
 # geom_errorbar(aes(ymin=meanzscore-sd, ymax=meanzscore+sd), width=0) +
  geom_line(ZNF467_mean, mapping=aes(x=Treatmentduration, y=meanzscore, colour=Treatment,linetype=Treatment), size=1.1) + 
  geom_point(size=2)+
  theme_classic()+
  scale_color_brewer(palette="Paired")+
  scale_x_continuous(breaks = c(0.5,1,2,4))+
  facet_wrap(~CommonName, ncol = 3,scales = "free")+
  ylab("rlog transformed gene expression") +xlab("treatment duration (hr)")+
  theme(axis.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
  theme(axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold") )+
  guides(color = guide_legend(override.aes = list(size = 1))) + #change the size of legend
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))+theme(legend.position = "none")+#change legend text fo
  theme(strip.background = element_rect(colour="light grey", fill="light grey", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size=12, color="black",
                                    face="bold.italic"))#change the apperance of the rectangle around facet label

#import matrix and plot motif
ZNF467_motifmatrix <- read_homer(file="ZNF467.motif")

ZNF467_motif <- convert_type(ZNF467_motifmatrix, "PPM")
## Only need the matrix itself
ZNF467_motif <- ZNF467_motif["motif"]
## ggseqlogo:
p9_4<- ggplot()+
  geom_logo(data=ZNF467_motif, method="probability",seq_type = "DNA",stack_width = 0.95,
            rev_stack_order = F)+theme_classic()+
  scale_x_continuous(breaks=seq(1,12,1))+
  scale_y_continuous(breaks=seq(0,1,0.25))+
  geom_segment(aes(x=1,xend=12,y=-Inf,yend=-Inf))+
  geom_segment(aes(y=0,yend=1,x=-Inf,xend=-Inf))+
  ylab("probability") +xlab("position")+
  ggtitle("ZNF467(Zf)/HEK293-ZNF467.GFP-ChIP-Seq(GSE58341)/Homer")+
  theme(axis.line=element_blank(),
        axis.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text=element_text(size=11),plot.title = element_text(hjust = 0.5,size = 14, face = "bold"))+
  theme(plot.title = element_text(size=12, color="black",
                                  face="bold"))#change the apperance of the rectangle around facet label

pp4 <- ggarrange(p9_3, p9_4,
                 ncol=1,
                 nrow=2,
                 align = "v",
                 labels=c("(c)","(d)"))

pp4

#final plot
p9 <- ggarrange(pp3,pp4, ncol=2,nrow = 1)
p9

############################################################################################
#create a new powerpoint document
setwd("RNAseq")
library(officer)
library(rvg)
doc <- read_pptx("final results.pptx")
doc <- add_slide(doc, 'Title and Content', 'Office Theme')

fig <- dml(ggobj =ppp4)
#add the plot
#doc <- ph_with(doc, pp, location = ph_location_fullsize())
doc <- ph_with(doc, fig, location = ph_location(width=13,height=6))
#write the document to a file
print(doc, target = "final results.pptx")
#############################################################################################

sessionInfo()