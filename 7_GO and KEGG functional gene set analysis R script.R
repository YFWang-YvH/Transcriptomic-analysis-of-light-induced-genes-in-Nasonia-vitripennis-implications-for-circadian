##Using StringTie with DESeq2 analyzing Nasonia RNAseq data
rm(list=ls()) #clear all objects
gc()

#dataset: 
#Nasonia RNAseq data
#experimental light treatment and dark control
#4 different time point: 0.5h 1h 2h 4h
#8 sample groups
#3 replicates per experimental sample group, 2 replicates per control sample group
#20 samples in total

#### GO analysis with TopGo ####
#over representation analysis 
#we determine the significance of the overlap between a list of genes of interest and a gene set by a statistical hypothesis test (Fisher's Exact Test)

# load packages
library("topGO")
library(ALL)
library(reshape2)
library(plyr)
library(ggplot2)
library(GO.db)
library(AnnotationDbi)
library(BiocGenerics)
library(ViSEAGO)
library("plotly")
library(clusterProfiler)
library(R.utils)
library(forcats)
library(enrichplot)
library(ggplotify)

# load files of GO annotation for topGO
GOTERMS14 <- read.csv("GOTERMS14_10122021.csv", header = T)
geneID2GO <- readMappings("GOTERMS14_10122021.txt")
geneNames <- names(geneID2GO)

DEG <- read.csv("DEgenelist3_23012022.csv")
DEG <- as.data.table(DEG)

#filter genes into up and down regulate genes for each time point
d0.5genelistup <- DEG %>% dplyr::filter(time=="0.5") %>% drop_na(FDR) %>% filter(FDR <0.05) %>% 
  filter(logFC > 0)
d0.5genelistup <- d0.5genelistup[!is.na(d0.5genelistup$NCBI.GeneID)]#270
d0.5genelistup <- as.data.frame(d0.5genelistup)
d0.5genelistup$number <- 1

d0.5genelistdown <- DEG %>% dplyr::filter(time=="0.5") %>% drop_na(FDR) %>% filter(FDR <0.05) %>% 
  filter(logFC < 0)
d0.5genelistdown <- d0.5genelistdown[!is.na(d0.5genelistdown$NCBI.GeneID)]#270
d0.5genelistdown <- as.data.frame(d0.5genelistdown)
d0.5genelistdown$number <- 2

d1genelistup <- DEG %>% dplyr::filter(time=="1") %>% drop_na(FDR) %>% filter(FDR <0.05)%>% 
  filter(logFC > 0)
d1genelistup <- d1genelistup[!is.na(d1genelistup$NCBI.GeneID)]#363
d1genelistup <- as.data.frame(d1genelistup)
d1genelistup$number <- 3

d1genelistdown <- DEG %>% dplyr::filter(time=="1") %>% drop_na(FDR) %>% filter(FDR <0.05)%>% 
  filter(logFC < 0)
d1genelistdown <- d1genelistdown[!is.na(d1genelistdown$NCBI.GeneID)]#363
d1genelistdown <- as.data.frame(d1genelistdown)
d1genelistdown$number <- 4

d2genelistup <- DEG %>% dplyr::filter(time=="2") %>% drop_na(FDR) %>% filter(FDR <0.05)%>% 
  filter(logFC > 0)
d2genelistup <- d2genelistup[!is.na(d2genelistup$NCBI.GeneID)]#1386
d2genelistup <- as.data.frame(d2genelistup)
d2genelistup$number <- 5

d2genelistdown <- DEG %>% dplyr::filter(time=="2") %>% drop_na(FDR) %>% filter(FDR <0.05)%>% 
  filter(logFC < 0)
d2genelistdown <- d2genelistdown[!is.na(d2genelistdown$NCBI.GeneID)]#1386
d2genelistdown <- as.data.frame(d2genelistdown)
d2genelistdown$number <- 6


d4genelistup <- DEG %>% dplyr::filter(time=="4") %>% drop_na(FDR) %>% filter(FDR <0.05) %>% 
  filter(logFC > 0)
d4genelistup <- d4genelistup[!is.na(d4genelistup$NCBI.GeneID)]#145
d4genelistup <- as.data.frame(d4genelistup)
d4genelistup$number <- 7

d4genelistdown <- DEG %>% dplyr::filter(time=="4") %>% drop_na(FDR) %>% filter(FDR <0.05) %>% 
  filter(logFC < 0)
d4genelistdown <- d4genelistdown[!is.na(d4genelistdown$NCBI.GeneID)]#145
d4genelistdown <- as.data.frame(d4genelistdown)
d4genelistdown$number <- 8


DEG_Fig3 <- rbind(d0.5genelistup,d0.5genelistdown,
                  d1genelistup,d1genelistdown,
                  d2genelistup,d2genelistdown,
                  d4genelistup,d4genelistdown
                  )

## Step 1 GO analysis for each time point and for up and down regulation separately
# do it in a loop
# gene universe = all genes in dataset
universalgene <- read.csv("gene_count_matrix_annotatedstringtie1_23012022.csv")
universalgenelist <- universalgene[!universalgene$Genename=="gene-",]#15696, some has no annotation
universalgenelist <- universalgenelist[,c(4)]
universalgenelist <- as.list(unique(universalgenelist))#15236

# Explore topGO environments
BPterms <- ls(GOBPTerm) # list of GO terms from BP ontology
head(BPterms)

# Choose ontology
ont <- "MF" #ontologies Biological Processes (MF), Molecular Function (MF), Cellular Component (CC)

# create lists to store results in #
GOdfs_list <- as.list(unique(DEG_Fig3$number))
names(GOdfs_list) <- unique(DEG_Fig3$number)

FisherRes_list <- as.list(unique(DEG_Fig3$number))
names(FisherRes_list) <- unique(DEG_Fig3$number)

# LOOP THROUGH MODULES = 8x ####
for (mod in 1:length(GOdfs_list)){
  module <- names(GOdfs_list)[mod]
  res <- filter(DEG_Fig3, number==module)$Symbol %>% unique
  
  # Note for each gene in the universe whether it's a DEG or not
  table(universalgenelist %in% res)
  geneList <- factor(as.integer(universalgenelist %in% res))
  names(geneList) <- c(universalgenelist)
  str(geneList) # named gene list DEG yes or no
  table(geneList)
  
  sampleGOdata <- new("topGOdata",
                      description = module, # optional 
                      ontology = ont, # GO graph
                      allGenes = geneList, # gene universe as logical DEG yes or no
                      #geneSel = topDiffGenes, # function to select genes of interest based on score for GSEA
                      nodeSize = 5, #prune GO hierarchy from terms with less than X annotated genes (5-10 recommended)
                      annot = annFUN.gene2GO, gene2GO=geneID2GO) #function used to extract gene-to-GO mappings from affyLib object
  GOdfs_list[[mod]] <- sampleGOdata
  GOdfs_list[[mod]] # see summary of object
  
  # Overrepresentation analysis ####
  resultFisher <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher") 
  FisherRes_list[[mod]] <- resultFisher # see summary of object
  FisherRes_list[[mod]]
  
  # gather results in a table #### gather all GO results even if it's not significant
  allRes <- GenTable(sampleGOdata, elimFisher = resultFisher,
                     orderBy = "elimFisher", ranksOf = "elimFisher",topNodes = length( usedGO(sampleGOdata)) )
  allRes$number <- module
  colnames(allRes)[6] <- "elimFisher"
  allRes
  
  # Expand table with gene level info ####
  # pull GO2Gene annotation for significant GOterms
  
  if(length(subset(allRes, elimFisher<=0.01)$GO.ID)>0){
    # retrieve genes2GO list from the "expanded" annotation in GOdata
    allGO <- genesInTerm(sampleGOdata) # pull the GO2Gene annotation
    head(allGO)
    
    # pull GO2Gene annotation for DEGs of set
    DEGs_GOannot <- lapply(allGO,function(x) x[x %in% names(geneList[geneList==1])])
    DEGs_GOannot <- DEGs_GOannot[lapply(DEGs_GOannot,length)>0] # remove GO terms without hits
    head(DEGs_GOannot)
    
    DEGs_GOannot[[allRes$GO.ID[1]]]
    DEGs_signGO <- DEGs_GOannot[subset(allRes, elimFisher<=0.01)$GO.ID]
    DEGs_signGO <- DEGs_signGO[lapply(DEGs_signGO,length)>0] # remove GO terms without hits
    head(DEGs_signGO)
    
    gen <- lapply(DEGs_signGO, as.data.frame) # format as data frame with one row per gene
    gen <- mapply(cbind, gen, "GO.ID"=names(gen), SIMPLIFY=F)
    gen <- lapply(names(gen), function(x) setNames(gen[[x]], c("Symbol", "GO.ID")) )
    gen <- data.table::rbindlist(gen)
    head(gen)
    head(allRes)
    DEG1 <- DEG[,c(11,12,13,18)]
    DEG1 <- as.data.table(DEG1)
    DEG1 <- unique(DEG1)
    allRes_ext <- merge(allRes, gen, by="GO.ID") # complete table with one row per gene
    allRes_ext <- merge(allRes_ext, DEG1, by="Symbol") # add annotation to the genes
    head(allRes_ext)
    
    
  } else{
    
    allRes_ext <- character(0)
  }
  
  
  # Store results
  head(allRes)
  head(allRes_ext)
  
  # end loop
  assign(paste("allRes", mod, sep=""), allRes)
  assign(paste("allRes_ext", mod, sep=""), allRes_ext)
}

head(GOdfs_list)

allRes_MF <- bind_rows(allRes1,allRes2,allRes3,allRes4,allRes5,
                       allRes6,allRes7,allRes8)
allRes_MF$category <- "MF"
allRes_ext_MF<- bind_rows(allRes_ext1, allRes_ext3,allRes_ext5,
                          allRes_ext6,allRes_ext7,allRes_ext8)
allRes_ext_MF$category <- "MF"

# do the GO analysis for BP, MP, CC then combine everything together
allRes_final <- rbind(allRes_BP, allRes_MF,allRes_CC)
allRes_ext_final <- rbind(allRes_ext_BP, allRes_ext_MF,allRes_ext_CC)
#fwrite(allRes_ext_final,"allRes_ext_timepoint_updownGO.csv")

## Step 2 Additionally, we also performed GO analysis based on time course clusters from TCseq, we adjusted the above model and redid GO analysis

## Step 3 GO analysis with ViSEAGO 
x <- as.list(GO.db[GOTERMS14.2$GOID])
x <- as.data.frame(x)
x<-AnnotationDbi::select(GO.db,columns=columns(GO.db),keys=keys(GO.db))

GOTERMS14.1<- GOTERMS14 %>%
  mutate(GO=strsplit(GO, ",")) %>%
  unnest(GO)
GOTERMS14.1$taxid <- "N.vitripennis"
GOTERMS14.1  <-  GOTERMS14.1 %>% dplyr::rename(gene_id=ref_gene_id, GOID=GO) %>%
  dplyr::mutate(gene_symbol=gene_id, evidence="IEA") %>%
  relocate(taxid, gene_id, gene_symbol, GOID, evidence)


GOTERMS14.2 <- GOTERMS14.1[!duplicated(GOTERMS14.1[,3:4]),]
sum(duplicated(GOTERMS14.2[,3:4]))
str(GOTERMS14.2)
GOTERMS14.2$GOID <- as.factor(GOTERMS14.2$GOID)

x$GOID <- as.factor(x$GOID)
x <- as.data.table(x)
GOTERMS14.2 <- as.data.table(GOTERMS14.2)
GOTERMS14.3 <- x[GOTERMS14.2, on = 'GOID']
sum(is.na(GOTERMS14.3$DEFINITION))
GOTERMS14.3 <- GOTERMS14.3[complete.cases(GOTERMS14.3),]
#fwrite(GOTERMS14.3,"GOTERMS14.3complete.txt", sep = "\t")

# connect to Custom file
Custom<-ViSEAGO::Custom2GO(file="GOTERMS14.3complete.txt")
str(Custom)

# Display table of available organisms with Custom
ViSEAGO::available_organisms(Custom)

# load GO annotations from Custom
myGENE2GO<-ViSEAGO::annotate(
  "N.vitripennis",
  Custom
)

# the following analysis is done for each subset of DEGs (either by timepoints or by clusters)
# example of the analysis we performed

#d0.5h dataset
d0.5gene <- read.csv("dres0.5Hashrfilt.csv")
d0.5gene <- d0.5gene %>% dplyr::rename(MSTRG=GeneID)
dID <- strsplit(d0.5gene$MSTRG, ' ')
dID <- as.data.frame((dID))
dID <- as.data.frame(t(dID))
d0.5gene$GeneID <- dID$V1
d0.5gene$Locus <- dID$V2

#d1h dataset
d1gene <- read.csv("dres1Hashrfilt.csv")#13308
d1gene <- d1gene %>% dplyr::rename(MSTRG=GeneID)
dID <- strsplit(d1gene$MSTRG, ' ')
dID <- as.data.frame((dID))
dID <- as.data.frame(t(dID))
d1gene$GeneID <- dID$V1
d1gene$Locus <- dID$V2

#d2h dataset
d2gene <- read.csv("dres2Hashrfilt.csv")#12116
d2gene <- d2gene %>% dplyr::rename(MSTRG=GeneID)
dID <- strsplit(d2gene$MSTRG, ' ')
dID <- as.data.frame((dID))
dID <- as.data.frame(t(dID))
d2gene$GeneID <- dID$V1
d2gene$Locus <- dID$V2

#d4h dataset
d4gene<- read.csv("dres4Hashrfilt.csv")#13010
d4gene <- d4gene %>% dplyr::rename(MSTRG=GeneID)
dID <- strsplit(d4gene$MSTRG, ' ')
dID <- as.data.frame((dID))
dID <- as.data.frame(t(dID))
d4gene$GeneID <- dID$V1
d4gene$Locus <- dID$V2

#named vector of type factor 
#0.5
allGenes0.5 <- d0.5gene$Locus
#named vector of type numeric
#feature2:numeric vector
geneList0.5 <- as.numeric(d0.5gene$FDR)
#feature2:name vector
names(geneList0.5) <- as.character(d0.5gene$Locus)

topDiffGenes <- function(allScore){return(allScore < 0.05)}
x <- topDiffGenes(geneList0.5)
sum(x)
head(geneList0.5)
str(geneList0.5)

selection0.5 <- d0.5genelist$Locus

# 1.Create topGOData object

BP_0.5<-ViSEAGO::create_topGOdata(
  geneSel=selection0.5,
  allGenes=allGenes0.5,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

BP_1<-ViSEAGO::create_topGOdata(
  geneSel=selection1,
  allGenes=allGenes1,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

BP_2<-ViSEAGO::create_topGOdata(
  geneSel=selection2,
  allGenes=allGenes2,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

BP_4<-ViSEAGO::create_topGOdata(
  geneSel=selection4,
  allGenes=allGenes4,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

# Create topGOData object
GOdata_0.5BP <- new("topGOdata",
                    description = "LT vs DC 0.5hr BP", ontology = "BP",
                    allGenes = geneList0.5,geneSel=topDiffGenes,nodeSize=5,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_1BP <- new("topGOdata",
                  description = "LT vs DC 1hr BP", ontology = "BP",
                  allGenes = geneList1,geneSel=topDiffGenes,nodeSize=5,
                  annot =annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_2BP <- new("topGOdata",
                  description = "LT vs DC 2hrs BP", ontology = "BP",
                  allGenes = geneList2,geneSel=topDiffGenes,nodeSize=5,
                  annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_4BP <- new("topGOdata",
                  description = "LT vs DC 4hrs BP", ontology = "BP",
                  allGenes = geneList4,geneSel=topDiffGenes,nodeSize=5,
                  annot = annFUN.gene2GO, gene2GO = geneID2GO)

#description(GOdata_0.5)
#graph(GOdata_0.5)
result_D0.5BP <- runTest(GOdata_0.5BP, algorithm = "elim", statistic = "fisher")
result_D1BP  <- runTest(GOdata_1BP, algorithm = "elim", statistic = "fisher")
result_D2BP  <- runTest(GOdata_2BP, algorithm = "elim", statistic = "fisher")
result_D4BP  <- runTest(GOdata_4BP, algorithm = "elim", statistic = "fisher")

# merge results from topGO
BP_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    "LT VS DC 0.5 hr"=c("BP_0.5","result_D0.5BP"),
    "LT VS DC 1 hr"=c("BP_1","result_D1BP"), 
    "LT VS DC 2 hrs"=c("BP_2","result_D2BP"),
    "LT VS DC 4 hrs"=c("BP_4","result_D4BP")))

#2.combine enriched GO terms
# display the merged table
ViSEAGO::show_table(BP_sResults)

# print the merged table in a file
ViSEAGO::show_table(
  BP_sResults,
  "BP_sResults.xls"
)

#3.graphs of GO enrichment tests
# count significant (or not) pvalues by condition
ViSEAGO::GOcount(BP_sResults)

# display interactions
ViSEAGO::Upset(
  BP_sResults,
  file="BPOLexport.xls"
)

#4.GO terms semantic similarity
# initialyse 
myGOs<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults
)

# compute all available Semantic Similarity (SS) measures
myGOs<-ViSEAGO::compute_SS_distances(
  myGOs,
  distance="Wang"
)

#5.visualization and interpretation of enriched GO TERMS
#5.1 Multi dimensional scaling of GO terms
# display MDSplot
ViSEAGO::MDSplot(myGOs)


#5.2 clustering heatmap of GO terms
# GOterms heatmap with the default parameters
Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=TRUE,
  showGOlabels=TRUE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =2
      )
    )
  ),
  samples.tree=NULL
)

# Display the clusters-heatmap
p1 <- ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms"
)

p1
orca(p1,"VISEAGO_BP.svg", width = 10 * 96, height = 8 * 96)

ViSEAGO::show_table(Wang_clusters_wardD2)

# Print the clusters-heatmap table
ViSEAGO::show_table(
  Wang_clusters_wardD2,
  "cluster_heatmap_Wang_wardD2_BP.xls"
)

#5.3Multi dimensional scaling of GO terms
# display colored MDSplot
ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOterms"
)

#6. visualization and interpretation of GO clusters
#6.1 compute semantic similarity between GO clusters
# calculate semantic similarites between clusters of GO terms
Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
  Wang_clusters_wardD2,
  distance=c("max", "avg","rcmax", "BMA")
)
# build and highlight in an interactive MDSplot grouped clusters for one distance object
ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOclusters"
)

#6.2 GO clusters semantic similarities heatmap
# GOclusters heatmap
Wang_clusters_wardD2<-ViSEAGO::GOclusters_heatmap(
  Wang_clusters_wardD2,
  tree=list(
    distance="BMA",
    aggreg.method="ward.D2"
  )
)
# sisplay the GOClusters heatmap
p2 <- ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOclusters"
)
p2
orca(p2,"VISEAGO_BP_heatmap2.svg", width = 10 * 96, height = 8 * 96)

#### KEGG pathway analysis ####
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
# kk_final1 <- subset(kk_final, kk_final$qvalue<0.01)
# mkk_final1 <- subset(mkk_final, mkk_final$qvalue<0.01)
#fwrite(kk_final, "KEGG_mkk_cluster_analysis_timepoint.csv")

#biological theme comparison
DEGSaple <- list("0.5"=d0.5G, "1"=d1G,"2"=d2G,"4"=d4G)
xx <- compareCluster(DEGSaple, fun="enrichKEGG",
                     organism="nvi", pvalueCutoff=0.05)
xx <- pairwise_termsim(xx)                     

p1 <- emapplot(xx, legend_n=2) 

kk_final1 <- mutate(kk_final, geneRatio = parse_ratio(GeneRatio)) %>%
  arrange(desc(geneRatio))
kk_final1 <- mutate(kk_final1, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ggplot(kk_final1,
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("rich factor") +
  ylab(NULL) + 
  ggtitle("Enriched Disease Ontology")
p2 <- as.ggplot(~dotplot(xx))
p3 <-  cnetplot(xx)

#Additionally we also performed KEGG pathway analysis with the 10 time series clusters, however, we did not find any significant KEGG pathway

sessionInfo()

