library(scran)
library(dplyr)
library(CONICSmat)


glioma<- Read10X(data.dir = "/scRNA_Seq/10x/cnv/glioma/GRCh38/")
glioma<-CreateSeuratObject(counts = glioma, project = "10X")

glioma<- NormalizeData(glioma, normalization.method = "LogNormalize", scale.factor = 100000)

regions=read.table("/database/CONICSmat/chromosome_arm_positions_grch38.txt",sep="\t",row.names = 1,header = T)
dim(glioma_CNV)
gene_pos=getGenePositions(rownames(glioma_CNV))
glioma_CNV=filterMatrix(glioma_CNV,gene_pos[,"hgnc_symbol"],minCells=10)
normFactor=calcNormFactors(glioma_CNV)
l=plotAll(glioma_CNV,normFactor,regions,gene_pos,"glioma_CNVs")

lrbic=read.table("glioma_CNVs_BIC_LR.txt",sep="\t",header=T,row.names=1,check.names=F)
colnames(lrbic)
candRegions=rownames(lrbic)[which(lrbic[,"BIC difference"]>300 & lrbic[,"LRT adj. p-val"]<0.001)]

plot.new()
png("glioma_CNVs_plotHistogram_update.png",width = 1000,height = 600)
hi=plotHistogram(l[,candRegions],glioma_CNV,clusters=4,zscoreThreshold=4)
dev.off()

normal= which(hi==3)
tumor=which(hi!=3)
redu=plotAll(glioma_CNV,normFactor,regions[candRegions,],gene_pos,"glioma_CNVs_with_info.pdf",normal=normal,tumor=tumor)

plot.new()
png("glioma_CNVs_plotChromosomeHeatmap.png",width = 1000,height = 600)
plotChromosomeHeatmap(glioma_CNV,normal = normal, plotcells = 1:(length(normal)+length(tumor)), gene_pos = gene_pos, windowsize = 121,thresh = 0.4,expThresh = 0.1,chr=T)##colo=c(rep("blue",length(normal)),rep("red",length(tumor))
dev.off()
  

