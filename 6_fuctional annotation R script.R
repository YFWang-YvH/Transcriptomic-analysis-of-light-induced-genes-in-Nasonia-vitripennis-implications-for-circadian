# Tidy up blasting and GO annotation files for RNAseq
rm(list=ls()) #clear all objects

# load packages
library(tidyverse)
library(dplyr)
library(data.table)
library(stringr)

### Step 1 We combined blast results and pulled additional annotation for DEGs only ###

# read sprot and trembl and nr blast file
sprotblast <-read.csv("uniprot_sprot.blast.tsv", sep = "\t", header = FALSE) #884045
sprotblast <- sprotblast %>% dplyr::rename(qseqid=V1, sseqid=V2, qframe=V3, qlen=V4, stitle=V5, salltitles=V6, evalue=V7)
sum(is.na(sprotblast$qframe)) #39737
sprotblast <- subset(sprotblast, !is.na(sprotblast$qframe)) #844308

tremblblast <-read.csv("uniprot_trembl.blast.tsv", sep = "\t", header = FALSE) #1603799
tremblblast <- tremblblast %>% dplyr::rename(qseqid=V1, sseqid=V2, qframe=V3, qlen=V4, stitle=V5, salltitles=V6, evalue=V7)
sum(is.na(tremblblast$qframe)) #0

NRblast <- read.csv("allhits.blast.tsv", sep = "\t", header = FALSE) #1615323
NRblast <- NRblast %>% dplyr::rename(qseqid=V1, sseqid=V2, qframe=V3, qlen=V4, stitle=V5, salltitles=V6, evalue=V7, qtitle=V8)
sum(is.na(NRblast$qframe)) #0

# merge blast files and make a new dataset with annotation from blast results for DEG genes only
sprotblast$blast <- 'sprotblast'
tremblblast$blast <- 'tremblblast'
NRblast$blast <- 'NRblast'
NRblast$qtitle <- NULL
sprotblast$qlen <- as.integer(sprotblast$qlen)
dblast <- bind_rows(sprotblast, tremblblast, NRblast)

dblast <- as.data.table(dblast)
DEgenes_annotation$qseqid <- DEgenes_annotation$ref_id
DEgenes_fullannotation <- dblast[DEgenes_annotation, on = 'qseqid']

sum(is.na(DEgenes_fullannotation$sseqid)) #122
subset <- subset(DEgenes_fullannotation, is.na(DEgenes_fullannotation$sseqid)==TRUE)
#fwrite(DEgenes_fullannotation, file=paste('DEgenes_fullannotation_20210614.csv', sep="")

### Step 2 Tidy up additional GO terms from several database ###

#GO annotation from interproscan
interproscan0 <- read.csv("myseq0.tsv", sep = "\t", header = FALSE) 
interproscan10000 <- read.csv("myseq100000.tsv", sep = "\t", header = FALSE)
interproscan20000 <- read.csv("myseq200000.tsv", sep = "\t", header = FALSE)
interproscan30000 <- read.csv("myseq300000.tsv", sep = "\t", header = FALSE)
interproscan <- rbind(interproscan0, interproscan10000,interproscan20000,interproscan30000)
interproscan <- as.data.frame(interproscan) %>% dplyr::rename("qpid"="V1", "Sequence_MD5_digest"="V2", "Sequence_length"="V3",
                                                              "Analysis"="V4", "Signature_accession"="V5", "Signature_description"="V6",
                                                              "start_location"="V7", "stop_location"="V8", "score"="V9", "status"="V10", "Data"="V11",
                                                              "Interproscan_accession"="V12", "Interproscan_description"="V13",
                                                              "GO"="V14", "Pathway"="V15")
interproscan <- interproscan[!(is.na(interproscan$GO) | interproscan$GO=="" | interproscan$GO == "-"), ]
#write.csv(interproscan, "interproscan.csv")

#GO annotation from pannzer2
GO_myseq0 <- read.csv("GO_myseq0.txt", header = TRUE, sep = "\t", colClasses = c("character"))
DE_myseq0 <- read.csv("DE_myseq0.txt", header = TRUE, sep = "\t", colClasses = c("character"))
GO_myse50000 <- read.csv("GO_myseq50000.txt",header = TRUE, sep = "\t", colClasses = c("character"))
DE_myseq50000 <- read.csv("DE_myseq50000.txt", header = TRUE, sep = "\t", colClasses = c("character"))
GO_myse100000 <- read.csv("GO_myseq100000.txt",header = TRUE, sep = "\t", colClasses = c("character"))
DE_myseq100000 <- read.csv("DE_myseq100000.txt", header = TRUE, sep = "\t", colClasses = c("character"))
GO_myse150000 <- read.csv("GO_myseq150000.txt",header = TRUE, sep = "\t", colClasses = c("character"))
DE_myseq150000 <- read.csv("DE_myseq150000.txt", header = TRUE, sep = "\t", colClasses = c("character"))
GO_myse200000 <- read.csv("GO_myseq200000.txt",header = TRUE, sep = "\t", colClasses = c("character"))
DE_myseq200000 <- read.csv("DE_myseq200000.txt", header = TRUE, sep = "\t", colClasses = c("character"))
GO_myse250000 <- read.csv("GO_myseq250000.txt",header = TRUE, sep = "\t", colClasses = c("character"))
DE_myseq250000 <- read.csv("DE_myseq250000.txt", header = TRUE, sep = "\t", colClasses = c("character"))
GO_myse300000 <- read.csv("GO_myseq300000.txt",header = TRUE, sep = "\t", colClasses = c("character"))
DE_myseq300000 <- read.csv("DE_myseq300000.txt", header = TRUE, sep = "\t", colClasses = c("character"))
GO_myse350000 <- read.csv("GO_myseq350000.txt",header = TRUE, sep = "\t", colClasses = c("character"))
DE_myseq350000 <- read.csv("DE_myseq350000.txt", header = TRUE, sep = "\t", colClasses = c("character"))
GO_myseq <- rbind(GO_myseq0, GO_myse100000,GO_myse150000,GO_myse200000,GO_myse250000,GO_myse300000,GO_myse350000)
DE_myseq <- rbind(DE_myseq0,DE_myseq100000,DE_myseq150000,DE_myseq200000,DE_myseq250000,DE_myseq300000,DE_myseq350000)
GO_myseq <- as.data.table(GO_myseq) #119555
DE_myseq <- as.data.table(DE_myseq)#19581
GO_DE_myseq <- GO_myseq[DE_myseq, on="qpid"]#138305
GO_DE_myseq$go <- "GO:"
GO_DE_myseq <- GO_DE_myseq %>% mutate(paste=paste(GO_DE_myseq$go, GO_DE_myseq$goid))
x <- GO_DE_myseq$paste
GO_DE_myseq$GO <- str_replace_all(x, regex("\\s*"), "")
GO_DE_myseq$goid <- NULL
GO_DE_myseq$go <- NULL
GO_DE_myseq$paste <- NULL
interproscan <- as.data.table(interproscan)#210975
#GO_annotation <- merge(interproscan, GO_DE_myseq, by="qpid", all=T,allow.cartesian=TRUE )

#extra GO annotation from uniprot retriving ID mapping
uniprot_entrezid <- read.csv("uniprot-Entrezid.csv", header=T) #11817
uniprot_entrezid1 <-  uniprot_entrezid[!(is.na(uniprot_entrezid$Gene.ontology.IDs) | uniprot_entrezid$Gene.ontology.IDs==""), ] #8424

uniprot_Refseq_proteinid <- read.csv("uniprot-Refseq protein id.csv", header = T) #7696
uniprot_Refseq_proteinid1 <-  uniprot_Refseq_proteinid[!(is.na(uniprot_Refseq_proteinid$Gene.ontology.IDs) | uniprot_Refseq_proteinid$Gene.ontology.IDs==""), ]#5454

#check which elements of first vector exist in second vector 
uniprot_Refseq_proteinid1$Entry %in% uniprot_entrezid1$Entry 
#find elements that exist in both first and second 
list <- intersect(uniprot_Refseq_proteinid1$Entry, uniprot_entrezid1$Entry )  
#find elements that exist only in first, but not in second vector
list_refseq <- setdiff(uniprot_Refseq_proteinid1$Entry, uniprot_entrezid1$Entry)   
list_entrezid <- setdiff(uniprot_entrezid1$Entry, uniprot_Refseq_proteinid1$Entry)   

#extract elements only in refseq 590 GO to combine with entrezid
uniprot_Refseq_proteinid2 <- subset(uniprot_Refseq_proteinid1, Entry %in% list_refseq)
list_protein <- uniprot_Refseq_proteinid2$Protein.product

allgenesfromRNAseq <- read.csv("allgenesfromRNAseq_02122021.csv", header = T)
allgenes <- allgenesfromRNAseq[, c(3,6)]
allgenes <- subset(allgenes, Protein.product %in% list_protein)
allgenes <- allgenes[!duplicated(allgenes$Entrezid),]

allgenes <- as.data.table(allgenes)
uniprot_Refseq_proteinid2 <- as.data.table(uniprot_Refseq_proteinid2)
uniprot_Refseq_proteinid2 <- uniprot_Refseq_proteinid2[allgenes, on="Protein.product"]
uniprot_Refseq_proteinid2 <- uniprot_Refseq_proteinid2 %>% relocate(Entrezid)
uniprot_Refseq_proteinid2 <- uniprot_Refseq_proteinid2[, -c(2)]

uniprot_annotation <- rbind(uniprot_entrezid1, uniprot_Refseq_proteinid2)
#fwrite(uniprot_annotation, "uniprot_annotation_10122021.csv")

uniprot_GO <- uniprot_annotation[, c(1, 10)]

uniprot_GO1 <- uniprot_GO[c(1:8424),]
uniprot_GO2 <- uniprot_GO[c(8425:9014),]
uniprot_GO2 <- uniprot_GO2 %>% dplyr::rename("ref_gene_id"="Entrezid","GO"="Gene.ontology.IDs")
uniprot_GO1$LOC <- "LOC"
uniprot_GO1 <- uniprot_GO1 %>% mutate(paste=paste(uniprot_GO1$LOC, uniprot_GO1$Entrezid))

uniprot_GO <- rbind(uniprot_GO1,uniprot_GO2)

uniprot_GO3 <- uniprot_GO %>%
  separate(GO, c("col1", "col2", "col3","col4","col5","col6","col7","col8","col9","col10","col11","col12","col13","col14","col15"), ";")

uniprot_GO3 <- as.data.table(uniprot_GO3)
uniprot_GO4 <- melt(uniprot_GO3, id.vars = "ref_gene_id" )
uniprot_GO5 <- uniprot_GO4 %>% drop_na(value) %>% rename("GO"="value")
uniprot_GO5 <- uniprot_GO5[,-2]

# tidy up GO terms from interproscan, Pannzer2, and uniprot
interproscan2 <- interproscan[!duplicated(interproscan[, c("qpid", "GO")],)]
interproscan3 <- subset(interproscan2, select=c("qpid", "GO"))

GO_DE_myseq2 <- GO_DE_myseq[!duplicated(GO_DE_myseq[, c("qpid", "GO")],)]#79606
GO_DE_myseq3 <- subset(GO_DE_myseq2, select=c("qpid", "GO"))

GOTERMS <- rbind(GO_DE_myseq3,interproscan3)
GOTERMS2 <- GOTERMS[!duplicated(GOTERMS[, c("qpid", "GO")],)]#126424

GOTERMS3 <- GOTERMS2 %>% separate(qpid, c("A", "B", "C","D"))
GO_DE_sub1 <- subset(GOTERMS3, GOTERMS3$A == "MSTRG")
GO_DE_sub1 <- GO_DE_sub1 %>% mutate(paste=paste(GO_DE_sub1$A, GO_DE_sub1$B,GO_DE_sub1$C, sep="."))
GO_DE_sub1 <- GO_DE_sub1 %>% dplyr::rename(qpid=paste) %>% subset(select=c("qpid","GO"))#38979

GO_DE_sub2 <- subset(GOTERMS3, GOTERMS3$A == "gene")
GO_DE_sub2 <- GO_DE_sub2 %>% mutate(paste=paste(GO_DE_sub2$A, GO_DE_sub2$B,sep="-"))
GO_DE_sub2 <- GO_DE_sub2 %>% dplyr::rename(qpid=paste) %>% subset(select=c("qpid","GO"))#330

GO_DE_sub3 <- subset(GOTERMS3, GOTERMS3$A == "XM"|GOTERMS3$A =="XR"|GOTERMS3$A == "NM"|GOTERMS3$A == "NR")#87115
GO_DE_sub3 <- GO_DE_sub3 %>% mutate(paste=paste(GO_DE_sub3$A, GO_DE_sub3$B,sep="_"))
GO_DE_sub3 <- GO_DE_sub3 %>% dplyr::rename(Q=paste) 
GO_DE_sub3 <- GO_DE_sub3 %>% mutate(paste=paste(GO_DE_sub3$Q, GO_DE_sub3$C,sep=".")) 
GO_DE_sub3 <- GO_DE_sub3 %>% dplyr::rename(qpid=paste) %>% subset(select=c("qpid","GO"))

GOTERMS5 <- rbind(GO_DE_sub1, GO_DE_sub2,GO_DE_sub3)
GOTERMS5$qry_id <- GOTERMS5$qpid
GOTERMS5 <- as.data.table(GOTERMS5)
merged.stringtie <- as.data.table(merged.stringtie)

GOTERMS6 <- merged.stringtie[GOTERMS5, on="qry_id"]
GOTERMS6 <- GOTERMS6 %>% subset(select=c("ref_gene_id","ref_id","qpid","GO"))
GOTERMS6 <- GOTERMS6[!duplicated(GOTERMS6[, c("ref_gene_id", "GO")],)]#42199
sum(GOTERMS6$GO =="GO:NA") #234
GOTERMS7 <- GOTERMS6[!GOTERMS6$GO=="GO:NA"] #41965

GOTERMS8 <- GOTERMS7 %>%
  group_by(ref_gene_id, grp = rleid(ref_gene_id)) %>% 
  summarise(GO = str_c(GO, collapse='|')) %>%
  ungroup %>%
  dplyr::select(-grp)#8602

GOTERMS9 <- GOTERMS8 %>%
  separate(GO, c("col1", "col2", "col3","col4","col5","col6","col7","col8","col9","col10","col11","col12","col13","col14","col15","col16","col17",
                 "col18","col19","col20","col21","col22","col23","col24","col25","col26","col27","col28","col29","col30",
                 "col31","col32","col33","col34","col35","col36","col37","col38","col39","col40","col41","col42","col43","col44","col45"), ",")

GOTERMS9 <- as.data.table(GOTERMS9)
GOTERMS10 <- melt(GOTERMS9, id.vars = "ref_gene_id" )
GOTERMS11 <- GOTERMS10 %>% drop_na(value) %>% rename("GO"="value")
GOTERMS11 <- GOTERMS11[,-2]
GOTERMS12 <- rbind(GOTERMS11, uniprot_GO5)#73566
GOTERMS13 <- GOTERMS12[!duplicated(GOTERMS12[, c("ref_gene_id", "GO")],)] #51737
list <- unique(GOTERMS13$ref_gene_id)#9817

GOTERMS14 <- GOTERMS13 %>%
  group_by(ref_gene_id) %>% 
  summarise(GO = paste0(GO, collapse=',')) %>% 
  ungroup()#9817
#write.csv(GOTERMS14, "GOTERMS14_10122021.csv")
