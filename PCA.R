rm(list=ls())

source("https://bioconductor.org/biocLite.R")

library(FactoMineR)
library(DESeq2)


setwd("/Users/xinling_li/Documents")

a1=read.table("count.csv",sep=',',header=T,row.names=1,check=F,comment.char="")

dim(a1)

colnames(a1)<-sub("*_1_1","",colnames(a1))
colnames(a1)


countdataraw=round(a1[rowSums(a1)>0,]) 

countdata<-as.matrix(countdataraw)
dim(countdata)


colnames(countdata)
condition<-factor(c("DSG_Tx", "DSG_Tx", "DSG_Tx", "DSG_Tx", "DSG_Tx", "DSG_Veh", "DSG_Veh", "WT_Veh", "DSG_Veh", "WT_Veh", "WT_Veh", "WT_Veh", "DSG_Veh", "DSG_Veh"))
#condition<-factor(c("GOOD", "GOOD", "GOOD", "GOOD", "GOOD", "GOOD", "GOOD", "GOOD", "GOOD", "GOOD", "POOR", "POOR", "POOR", "POOR", "POOR", "POOR", "POOR", "POOR", "POOR", "POOR", "POOR", "POOR"))

coldata <- data.frame(row.names=colnames(countdata), condition)

coldata

coldata$condition<-as.factor(as.character(coldata$condition))

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

dds=DESeq(dds,minReplicatesForReplace = 50)

rld=rlog(dds, blind=TRUE)

name<-"Plot"
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
dat.norm<-t(assay(rld)[select, ])
dat.norm<-data.frame(dat.norm,condition)
res.pca=PCA(dat.norm,ncp=5,scale.unit=T,graph=F,quali.sup=ncol(dat.norm))
pdf(paste(name,"_PCA.pdf",sep=""),16,12)
plot.PCA(res.pca,axes=c(1,2),habillage=ncol(dat.norm),cex=1, label='ind')
plot.PCA(res.pca,axes=c(1,2),habillage=ncol(dat.norm),cex=2,label='none')
plot.PCA(res.pca,axes=c(1,3),habillage=ncol(dat.norm),cex=1, label='ind')
plot.PCA(res.pca,axes=c(1,3),habillage=ncol(dat.norm),cex=2,label='none')
plot.PCA(res.pca,axes=c(2,3),habillage=ncol(dat.norm),cex=1, label='ind')
plot.PCA(res.pca,axes=c(2,3),habillage=ncol(dat.norm),cex=2,label='none')
dev.off()