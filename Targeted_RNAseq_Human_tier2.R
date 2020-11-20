library(DESeq2)
library(plotly)
library("FactoMineR")
library("rmarkdown")
library("vegan")
library("manhattanly")
library("DT")
library("heatmaply")
library(gplots)
rm(list=ls())
args=commandArgs(TRUE)
count_file <- args[1]
sample_info <-args[2]
comparison <- args[3]
project <- args[4]
a1=read.table(count_file, sep=',',header=T,row.names=1,check=F,comment.char="")
a2=read.table(sample_info,sep=',',header=T)
a3=read.table(comparison,sep=',',header=F)
comps <- as.matrix(a3) 
a2 <- a2[order( a2[,3], a2[,2] ),]

a1<-a1[,as.vector(a2[,1])]

colnames(a1) <- as.matrix(a2)[,2]
countdataraw=round(a1[rowSums(a1)>0,])
countdata<-as.matrix(countdataraw)
dim(countdata)
temp <- data.frame(countdata)
temp$Gene <- row.names(temp)
library("AnnotationDbi")
library("org.Hs.eg.db")
temp$symbol = mapIds(org.Hs.eg.db,
                          keys=temp$Gene, 
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
row.names(temp) <- paste(temp$Gene,temp$symbol, sep = "_")
n <- dim(a2)[1]
countdata <- as.matrix(temp[,1:n])
condition<-factor(as.matrix(a2)[,3])
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds=DESeq(dds,minReplicatesForReplace = 50)
norm=counts(dds,normalized = T)

rld=rlog(dds, blind=TRUE)
rv <- rowVars(assay(rld)) 
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
dat.norm<-t(assay(rld)[select, ])
dat.norm<-data.frame(dat.norm,condition)
res.pca=PCA(dat.norm,ncp=5,scale.unit=T,graph=F,quali.sup=ncol(dat.norm))

pdf(paste(project,"_PCA_all_samples.pdf",sep=""),16,12)
plot.PCA(res.pca,axes=c(1,2),habillage=ncol(dat.norm),cex=1)
plot.PCA(res.pca,axes=c(1,2),habillage=ncol(dat.norm),cex=2,label='none')
plot.PCA(res.pca,axes=c(1,3),habillage=ncol(dat.norm),cex=1)
plot.PCA(res.pca,axes=c(1,3),habillage=ncol(dat.norm),cex=2,label='none')
plot.PCA(res.pca,axes=c(2,3),habillage=ncol(dat.norm),cex=1)
plot.PCA(res.pca,axes=c(2,3),habillage=ncol(dat.norm),cex=2,label='none')
dev.off()

for (i in 1:dim(comps)[1]) {
  name <- paste0(comps[i,][2], "_vs_", comps[i,][3])  #"M_vs_MC"
  dir.create(name)
  out_folder <- paste0(getwd(),"/",name,"/")
  pattern <- c(comps[i,][2], comps[i,][3])
  countdata_1 <- countdata[, grep(paste(pattern, collapse="|"), colnames(countdata), value = TRUE)]
  condition<-factor(as.matrix(a2)[,3])
  condition <- factor(grep(paste(pattern, collapse="|"), condition, value = TRUE))
  coldata_1 <- data.frame(row.names=colnames(countdata_1), condition)
  dds_1 <- DESeqDataSetFromMatrix(countData=countdata_1, colData=coldata_1, design=~condition)
  dds_1=DESeq(dds_1,minReplicatesForReplace = 50)
  norm_1=counts(dds_1,normalized = T)
  rld_1=rlog(dds_1, blind=TRUE)
  rv_1 <- rowVars(assay(rld_1))
  select_1 <- order(rv_1, decreasing = TRUE)[seq_len(min(500,length(rv_1)))]
  dat.norm_1<-t(assay(rld_1)[select_1, ])
  dat.norm_1<-data.frame(dat.norm_1,condition)
  res.pca_1=PCA(dat.norm_1,ncp=5,scale.unit=T,graph=F,quali.sup=ncol(dat.norm_1))
  pdf(paste(out_folder,name,"_PCA_top500.pdf",sep=""),16,12)
  plot.PCA(res.pca_1,axes=c(1,2),habillage=ncol(dat.norm_1),cex=1)
  plot.PCA(res.pca_1,axes=c(1,2),habillage=ncol(dat.norm_1),cex=2,label='none')
  plot.PCA(res.pca_1,axes=c(1,3),habillage=ncol(dat.norm_1),cex=1)
  plot.PCA(res.pca_1,axes=c(1,3),habillage=ncol(dat.norm_1),cex=2,label='none')
  plot.PCA(res.pca_1,axes=c(2,3),habillage=ncol(dat.norm_1),cex=1)
  plot.PCA(res.pca_1,axes=c(2,3),habillage=ncol(dat.norm_1),cex=2,label='none')
  dev.off()

  assign(paste("res", i, sep = ""), results(dds_1, contrast = comps[i,]))
  assign(paste0("resdata",i), merge(as.data.frame(eval(parse(text=paste0("res", i)))), as.data.frame(norm_1), by="row.names", sort=FALSE))
  assign(paste("resdata", i, sep = ""), eval(parse(text=paste0("resdata", i)))[,-c(4,5)]) 
  write.csv(eval(parse(text=paste0("resdata", i))), file=paste(out_folder,name,"_DEGs_all.csv",sep=''),row.names=F)
  n1 <- dim(subset(eval(parse(text=paste0("resdata", i))), padj<0.05))[1]
  n2 <- dim(subset(eval(parse(text=paste0("resdata", i))), padj<0.01))[1]
  if (n1 < 2000 & n1 >= 200) {
    resSig=subset(eval(parse(text=paste0("resdata", i))), padj<0.05)
    write.csv(resSig,file=paste(out_folder,name,"_DEGs_padj0.05.csv",sep=''),row.names=F)
    sig_DEGs <- "p_adj < 0.05"
  }
  if (n1 >= 2000) {
    if (n2 < 2000) {
      resSig=subset(eval(parse(text=paste0("resdata", i))), padj<0.01)
      write.csv(resSig,file=paste(out_folder,name,"_DEGs_padj0.01.csv",sep=''),row.names=F)
      sig_DEGs <- "p_adj < 0.01"
    } else {
      resSig=subset(eval(parse(text=paste0("resdata", i))), padj<0.01 & abs(log2FoldChange)>1)
      write.csv(resSig,file=paste(out_folder,name,"_DEGs_padj0.01_FC2.csv",sep=''),row.names=F)
      sig_DEGs <- "p_adj < 0.01 and |FC| > 2"
    }
  }
  if (n1 < 200 & n1 >= 20) {
    resSig=subset(eval(parse(text=paste0("resdata", i))), padj<0.1)
    write.csv(resSig,file=paste(out_folder,name,"_DEGs_padj0.1.csv",sep=''),row.names=F)
    sig_DEGs <- "p_adj < 0.1"
  }
  if (n1 < 20) {
    resSig=subset(eval(parse(text=paste0("resdata", i))), pvalue<0.05)
    write.csv(resSig,file=paste(out_folder,name,"_DEGs_p0.05.csv",sep=''), row.names=F)
    sig_DEGs <- "p < 0.05"
  }
  rmarkdown::render("Interactive_report_RNAseq.Rmd", params = list(data = a1, info = a2, comparison = a3, project = project), output_file = paste0(out_folder,name, ".html"), output_dir=out_folder)
  i=i+1 
};rm(i)

