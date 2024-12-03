rm(list=ls())

library(wesanderson)
library(BiocManager)
library(topGO)
library(scales)
library("superheat")
library("RColorBrewer")
library("factoextra")
library("FactoMineR")
library(RColorBrewer)
library("UpSetR")
library("gplots")
library("tidyverse")
library("reshape2")
library(GGally)
library("microbiome")
library(vegan)
library("superheat")
library(DESeq2)
library(dplyr)
library(readr)
library(magrittr)
library(tximport)
library(jsonlite)
library(pheatmap)
library(viridis)

getwd()
setwd("~/Dobslab Dropbox/Adam Dobson/students/LlanwarneFrances/frances_manuscript")
getwd()
list.files()
OTUs<-read.table("FeatureTableLevel7_no_endos_taxons.txt", sep="\t", header=T)
OTUs
head(OTUs) 
str(OTUs)
dim(OTUs)
colnames(OTUs)
otuIndex <- OTUs$OTUID
rownames(OTUs)
#OTUs$OTU.ID

	#remove unassigned reads
OTUs <- OTUs[1:nrow(OTUs)-1,]

	#output the OTU index
write.table(otuIndex, file="otuIndex.txt", quote=F, row.names=F, col.names=F)

	#input OTU decode
otuDecode <- read.csv("otuDecode.csv")

rownames(OTUs)<-otuDecode$OTU_lowest

head(OTUs)
colnames(OTUs)
colnames(OTUs[,2:7])
OTUs<-OTUs[,2:7]
head(OTUs)
str(OTUs)



	#upset
OTU_bin <- data.frame(ifelse(OTUs>0, 1, 0))
upset(OTU_bin, nsets=6, matrix.color=1)

	#couccurence heatmap
heatmap.2(as.matrix(OTU_bin), trace="none", col=c("white", "grey"), key=F, dendro="none")
coreSpp <- which(rowSums(OTU_bin)==6)
otuIndex[coreSpp]
coreOTUs <- OTUs[coreSpp,]
#rownames(coreOTUs) <- c("g__Leuconostoc", "g__Lactococcus", "f__Acetobacteraceae", "g__Gluconobacter", "f__Enterobacteriaceae", "g__Brenneria")

	#are there correlated OTUs?
ggcorr(t(OTUs), method = c("pairwise", "spearman"), nbreaks = 11, digits=1, low="blue", mid="white", high="red")
ggcorr(t(coreOTUs), method = c("pairwise", "spearman"), nbreaks = 11, digits=1, low="blue", mid="white", high="red")
cor(t(coreOTUs))

	#abundance heatmaps
superheat(log(data.frame(coreOTUs)+1), 
	heat.lim=log(c(1e-10, 35000)+1),
	heat.col.scheme="red",
	grid.hline.col="black", 
	grid.vline.col="black", 
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	left.label.size=1.5,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white"
	)

superheat(log(OTUs+1), 
	pretty.order.rows=T,
	pretty.order.cols=T,
	row.dendrogram=T,
	col.dendrogram=T,
	heat.lim=log(c(1e-10, 35000)+1),
	heat.na.col="lightgrey",
	heat.col.scheme="red",
	grid.vline.col="black", 
	smooth.heat=T,
	grid.vline.size=0.1,
	grid.hline=F,
	left.label="none",
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white"
	)

	#alpha diversity

aDiv <- alpha(OTUs, index='shannon', zeroes=T)
aDiv

	#beta diversity
betaDiv <- betadiver(t(OTUs))
plot(betaDiv)

	#what is the distribution of total reads?
plot(density(log(rowSums(OTUs)+1)))
table(rowSums(OTUs)>10)
table(rowSums(OTUs)>7)
table(rowSums(OTUs)>3)
table(rowSums(OTUs)>4)

xOTU <- log(t(OTUs[rowSums(OTUs)>4,]+1))
#xxOTU <- xOTU / rowSums(xOTU)

## use factominer to make a pca, then use dimdesc to describe dimensions
dev.new()
thePca <- PCA(xOTU, scale.unit=T)
summary(thePca)
	#use the following in DESeq
fviz_pca_ind(thePca)
thePca$ind$coord
write.csv(thePca$ind$coord, "otuPCsToUseInDESeq.csv", row.names=T)
pcaDesc <- dimdesc(thePca, axes=1:5, proba=0.1)

str(pcaDesc, max.level=2)
str(pcaDesc$Dim.1$quanti)
pcaDesc$Dim.1$quanti
pcaDesc$Dim.2$quanti

	#so e.g. OTU 5 (log+1) should correlate PC1
plot(thePca$ind$coord[,1] ~ as.numeric(log(OTUs[rownames(OTUs)==rownames(pcaDesc$Dim.1$quanti)[1],]+1)))
	
	#have a look at PC2, where there is less erratic OTU variation
plot(thePca$ind$coord[,2] ~ as.numeric(log(OTUs[rownames(OTUs)==rownames(pcaDesc$Dim.2$quanti)[1],]+1)), ylab="PC 2", xlab="OTU count")	#yes, much more sensible

	#AD PCA plots

pdf("PCA.percentage.explained.pdf")
par(las=2)
bp <- barplot(c(0,thePca$eig[,2]), col=1, ylim=c(0,100), ylab="%age variance explained")
lines(x=bp, y=c(0,thePca$eig[,3]), col="red", lwd=2)
dev.off()

pdf("PCA.pairs.pdf")
pairs(thePca$ind$coord, pch=16, col=brewer.pal(6, "Dark2"), cex=2, las=1, lower.panel=NULL)
dev.off()
pdf("PCA.legend.pdf")
plot(1,1,col=0, axes=F)
legend("bottomleft", legend=rownames(thePca$ind$coord), col=brewer.pal(6, "Dark2"), pch=16)
dev.off()


	#PC1 (by AD)
pdf("otu_pc1.pdf")
superheat(log(OTUs[rownames(OTUs) %in% rownames(pcaDesc$Dim.1$quanti),]+1),
          yt=thePca$ind$coord[,1], 
          order.cols=order(thePca$ind$coord[,1]),
          yt.axis.name="PC1", 
          yt.plot.type="line",
          yt.num.ticks=5,
          row.dendrogram = T,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,8),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")
dev.off()          

	#PC2 (by AD)
pdf("otu_pc2.pdf")
superheat(log(OTUs[rownames(OTUs) %in% rownames(pcaDesc$Dim.2$quanti),]+1),
          yt=thePca$ind$coord[,2], 
          order.cols=order(thePca$ind$coord[,2]),
          yt.axis.name="PC2", 
          yt.plot.type="line",
          yt.num.ticks=5,
          row.dendrogram = T,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,8),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")
dev.off()

	#PC3 (by AD)
pdf("otu_pc3.pdf")
superheat(log(OTUs[rownames(OTUs) %in% rownames(pcaDesc$Dim.3$quanti),]+1),
          yt=thePca$ind$coord[,3], 
          order.cols=order(thePca$ind$coord[,3]),
          yt.axis.name="PC3", 
          yt.plot.type="line",
          yt.num.ticks=5,
          row.dendrogram = T,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,8),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")
dev.off()

	#PC4 (by AD)
pdf("otu_pc4.pdf")
superheat(log(OTUs[rownames(OTUs) %in% rownames(pcaDesc$Dim.4$quanti),]+1),
          yt=thePca$ind$coord[,4], 
          order.cols=order(thePca$ind$coord[,4]),
          yt.axis.name="PC4", 
          yt.plot.type="line",
          yt.num.ticks=5,
          row.dendrogram = T,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,8),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")
dev.off()          

	#PC5 (by AD)
pdf("otu_pc5.pdf")
superheat(log(OTUs[rownames(OTUs) %in% rownames(pcaDesc$Dim.5$quanti),]+1),
          yt=thePca$ind$coord[,5], 
          order.cols=order(thePca$ind$coord[,5]),
          yt.axis.name="PC5", 
          yt.plot.type="line",
          yt.num.ticks=5,
          row.dendrogram = T,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,8),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")
dev.off()


	#with p-values
superheat(log(OTUs[rownames(OTUs) %in% rownames(pcaDesc$Dim.1$quanti),]+1),
          yt=thePca$ind$coord[,1], 
          yr=(-log10(pcaDesc$Dim.1$quanti[,2])),
          order.cols=order(thePca$ind$coord[,1]),
          order.rows=order(-log10(pcaDesc$Dim.1$quanti[,2])),
          yt.axis.name="PC1", 
          yt.plot.type="line",
          yr.plot.type="bar",
          yr.axis.name="-log10(p)",
          yt.num.ticks=5,
          row.dendrogram = F,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,8),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")

superheat(log(OTUs[rownames(OTUs) %in% rownames(pcaDesc$Dim.2$quanti),]+1),
          yt=thePca$ind$coord[,2], 
          yr=(-log10(pcaDesc$Dim.2$quanti[,2])),
          order.cols=order(thePca$ind$coord[,2]),
          order.rows=order(-log10(pcaDesc$Dim.2$quanti[,2])),
          yt.axis.name="PC2", 
          yt.plot.type="line",
          yr.plot.type="bar",
          yr.axis.name="-log10(p)",
          yt.num.ticks=5,
          row.dendrogram = F,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,8),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")

superheat(log(OTUs[rownames(OTUs) %in% rownames(pcaDesc$Dim.3$quanti),]+1),
          yt=thePca$ind$coord[,3], 
          yr=(-log10(pcaDesc$Dim.3$quanti[,2])),
          order.cols=order(thePca$ind$coord[,3]),
          order.rows=order(-log10(pcaDesc$Dim.3$quanti[,2])),
          yt.axis.name="PC3", 
          yt.plot.type="line",
          yr.plot.type="bar",
          yr.axis.name="-log10(p)",
          yt.num.ticks=5,
          row.dendrogram = F,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,8),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")

superheat(log(OTUs[rownames(OTUs) %in% rownames(pcaDesc$Dim.4$quanti),]+1),
          yt=thePca$ind$coord[,4], 
          yr=(-log10(pcaDesc$Dim.4$quanti[,2])),
          order.cols=order(thePca$ind$coord[,4]),
          order.rows=order(-log10(pcaDesc$Dim.4$quanti[,2])),
          yt.axis.name="PC4", 
          yt.plot.type="line",
          yr.plot.type="bar",
          yr.axis.name="-log10(p)",
          yt.num.ticks=5,
          row.dendrogram = F,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,8),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")

superheat(log(OTUs[rownames(OTUs) %in% rownames(pcaDesc$Dim.5$quanti),]+1),
          yt=thePca$ind$coord[,5], 
          yr=(-log10(pcaDesc$Dim.5$quanti[,2])),
          order.rows=order(-log10(pcaDesc$Dim.5$quanti[,2])),
          order.cols=order(thePca$ind$coord[,5]),
          yt.axis.name="PC5", 
          yt.plot.type="line",
          yr.plot.type="bar",
          yr.axis.name="-log10(p)",
          yt.num.ticks=5,
          row.dendrogram = F,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,8),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")


##############
##	DESeq2	##
##############

# Getting metadata into R and manipulating to make my sample table metadata file
#sample_table <- read.csv("SraRunTable.csv") %>%
 # select(`Sample.Name`)

# Prep for 'files' and running tximport
# gene_map for tx2gene was created in terminal - see scripts

#sample_files <- paste0(
##  pull(sample_table, `Sample.Name`),
#  '/quant.sf')

#names(sample_files) <- pull(sample_table, `Sample.Name`)

gene_map <- read.csv("gene_map.csv",
                     col.names = c('FBtrid', 'FBgnid'))

#count_data <- tximport(files = sample_files, 
#         type = "salmon",
#         tx2gene = gene_map,)

#count_data[['counts']] # this is the data for the 'num reads' column in the quant.sf file

#count_data$abundance # this is the data for the 'abundance' column in the quant.sf file

#colnames(gene_map)

#write.csv(unique(gene_map$FBgnid), file = "Unique_gene_IDs.csv")

## Start DESeq2 ##

# Made new metadata file containing the PCs1-4 and their values for each sample
# followed the same steps as before to do the tximport

sample_table_PCs <- data.frame(thePca$ind$coord)
sample_table_PCs$SampleName <- rownames(thePca$ind$coord)
colnames(sample_table_PCs) <- gsub("Dim.", "PC", colnames(sample_table_PCs))
sample_table_PCs <- sample_table_PCs %>%
  select(`SampleName`, `PC1`, `PC2`, `PC3`, `PC4`)

	#check if scaling required for DESeq
apply(sample_table_PCs[,2:5], 2, sd)
apply(sample_table_PCs[,2:5], 2, mean)
sample_table_PCs_scaled <- apply(sample_table_PCs[,2:5], 2, scale, center=T, scale=T)
sample_table_PCs_scaled <- data.frame(SampleName=sample_table_PCs$SampleName, sample_table_PCs_scaled)

	#and how do they look?
pairs(sample_table_PCs_scaled[,2:5], pch=16, col=brewer.pal(6, "Dark2"), cex=2, las=1, lower.panel=NULL)

sample_files_PCs <- paste0(
  pull(sample_table_PCs, `SampleName`),
  '/quant.sf')

names(sample_files_PCs) <- pull(sample_table_PCs, `SampleName`)

count_data_PCs <- tximport(files = sample_files_PCs, type = "salmon", tx2gene = gene_map)

head(count_data_PCs[['counts']])

	# Making a deseq dataset with the data from tximport
deseq_dataset <- DESeqDataSetFromTximport(
	txi = count_data_PCs, 
	colData = sample_table_PCs_scaled, # colData is the metadata associated with the count data table
	design = ~ PC1 + PC2 + PC3 + PC4)

	#how many genes expressed per sample?
genesBeforeFilter <- apply(assay(deseq_dataset), 2, function(x){table(x==0)})

	#sample depth
sampleDepthBeforeFilter <- colSums(assay(deseq_dataset))

	#how many genes had one read detected?
allZero <- apply(assay(deseq_dataset), 1, function(x){all(x==0)})
table(allZero)

	#filter the dataset
notAllZero <- apply(assay(deseq_dataset), 1, function(x){all(x>0)})
table(notAllZero)

deseq_dataset <- deseq_dataset[notAllZero,]
graphics.off()

	#how many genes expressed per sample after filtering?
genesAfterFilter <- apply(assay(deseq_dataset), 2, function(x){table(x==0)})

	#sample depth
sampleDepthAfterFilter <- colSums(assay(deseq_dataset))

	#table of library properties
libProperties <- data.frame(
	data.frame(t(genesBeforeFilter)), 
	sampleDepthBeforeFilter, 
	genesAfterFilter, 
	sampleDepthAfterFilter)

libProperties
colnames(libProperties)[1:2] <- c("genesDetectedBeforeFiltering", "zeroesBeforeFiltering")
libProperties
write.csv(libProperties, "libraryProperties.csv", row.names=T, quote=F)

deseq_dataset # shows you the data object stats that you have created with colData etc

counts(deseq_dataset)

counts(deseq_dataset) [1:6, 1:6] 

	#are there variable transcripts in these data?
rlReads <- rlog(deseq_dataset)
coeffVars <- rowVars(assay(rlReads)) / rowMeans(assay(rlReads))
plot(density(coeffVars+1))
graphics.off()
coeffVarsRank <- order(coeffVars, decreasing=T)
scaleRl <- t(apply(assay(rlReads), 1, scale, center=T, scale=T))
topTranscriptsScaled <- scaleRl[coeffVarsRank,][1:1000,]

pdf("variableTranscripts_1k.pdf")
superheat(
	topTranscriptsScaled,
	col.dendrogram = T,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white")
dev.off()

CVs1k <- coeffVars[coeffVarsRank][1:1000]
ord <- hclust(dist(topTranscriptsScaled))$order
pdf("variableTranscripts_1k_withCVsidePlot.pdf")
superheat(
	topTranscriptsScaled[ord,],
	yr=CVs1k[ord],
	yr.axis.name="Coeff. of variation",
	col.dendrogram = T,
	yr.plot.type="line",
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white")
dev.off()

deseq_dataset_normalised <- estimateSizeFactors(deseq_dataset)

normalizationFactors(deseq_dataset_normalised)

counts(deseq_dataset_normalised, normalized = TRUE) [1:6, 1:6]

counts(deseq_dataset) [1:6, 1:6] # can see here that the numbers are different from before normalisation

boxplot(counts(deseq_dataset)) # box plot of counts from deseq dataset

	#cutoff for PADJ
pthresh <- 0.01

	#function to return DE genes
giveGene <- function(x, thresh){
	rownames(x)[which(x$padj <= pthresh)]
}

	## PC1 ##
pc1Deseq <- DESeq(deseq_dataset, test="LRT", reduced= ~ PC2 + PC3 + PC4)
respc1Deseq <- results(pc1Deseq)
table(is.na(respc1Deseq$padj))
table(is.na(p.adjust(respc1Deseq$pvalue, method="fdr")))
respc1Deseq$padj <- p.adjust(respc1Deseq$pvalue, method="fdr")
sum(respc1Deseq$padj <= pthresh, na.rm = TRUE)
pc1_DE <- giveGene(respc1Deseq, pThresh)
pc1_DE

	## PC2 ##
pc2Deseq <- DESeq(deseq_dataset, test="LRT", reduced= ~ PC1 + PC3 + PC4)
respc2Deseq <- results(pc2Deseq)
table(is.na(respc2Deseq$padj))
table(is.na(p.adjust(respc2Deseq$pvalue, method="fdr")))
respc2Deseq$padj <- p.adjust(respc2Deseq$pvalue, method="fdr")
sum(respc2Deseq$padj <= pthresh, na.rm = TRUE)
pc2_DE <- giveGene(respc2Deseq, pThresh)
pc2_DE

	#how do the p-vals look?
par(mfrow=c(1,2))
plot(density(-log10(respc2Deseq$padj)), "-log10 FDR")
plot(density(-log10(respc2Deseq$pvalue)), "-log10 pval")
graphics.off()

	## PC3 ##
pc3Deseq <- DESeq(deseq_dataset, test="LRT", reduced= ~ PC1 + PC2 + PC4)
respc3Deseq <- results(pc3Deseq)
table(is.na(respc3Deseq$padj))
table(is.na(p.adjust(respc3Deseq$pvalue, method="fdr")))
respc3Deseq$padj <- p.adjust(respc3Deseq$pvalue, method="fdr")
sum(respc3Deseq$padj <= pthresh, na.rm = TRUE)
pc3_DE <- giveGene(respc3Deseq, pThresh)
pc3_DE

	## PC4 ##
pc4Deseq <- DESeq(deseq_dataset, test="LRT", reduced= ~ PC1 + PC2 + PC3)
respc4Deseq <- results(pc4Deseq)
table(is.na(respc4Deseq$padj))
table(is.na(p.adjust(respc4Deseq$pvalue, method="fdr")))
respc4Deseq$padj <- p.adjust(respc4Deseq$pvalue, method="fdr")
sum(respc4Deseq$padj <= pthresh, na.rm = TRUE)
pc4_DE <- giveGene(respc4Deseq, pThresh)
pc4_DE

		#mean fdr values
mean(subset(respc1Deseq, padj <= pthresh)$padj)
mean(subset(respc2Deseq, padj <= pthresh)$padj)
mean(subset(respc3Deseq, padj <= pthresh)$padj)
mean(subset(respc4Deseq, padj <= pthresh)$padj)

	#pull correlated transcripts
pc1genes <- assay(rlReads[rownames(rlReads) %in% pc1_DE])
pc2genes <- assay(rlReads[rownames(rlReads) %in% pc2_DE])
pc3genes <- assay(rlReads[rownames(rlReads) %in% pc3_DE])
pc4genes <- assay(rlReads[rownames(rlReads) %in% pc4_DE])


	#gene decode
genes <- read.csv("FBgn-Sym-dm6.28.csv")
head(genes)
rownames(pc1genes) <- genes$gene_name[match(rownames(pc1genes), genes$gene_id)]
rownames(pc2genes) <- genes$gene_name[match(rownames(pc2genes), genes$gene_id)]
rownames(pc3genes) <- genes$gene_name[match(rownames(pc3genes), genes$gene_id)]
rownames(pc4genes) <- genes$gene_name[match(rownames(pc4genes), genes$gene_id)]
pc1genes
pc2genes
pc3genes
pc4genes

#write to csv
respc1Deseq$gene <- genes$gene_name[match(rownames(respc1Deseq), genes$gene_id)]
respc2Deseq$gene <- genes$gene_name[match(rownames(respc2Deseq), genes$gene_id)]
respc3Deseq$gene <- genes$gene_name[match(rownames(respc3Deseq), genes$gene_id)]
respc4Deseq$gene <- genes$gene_name[match(rownames(respc4Deseq), genes$gene_id)]
write.csv(respc1Deseq, "respc1Deseq.csv", quote=F, row.names=T)
write.csv(respc1Deseq, "respc2Deseq.csv", quote=F, row.names=T)
write.csv(respc1Deseq, "respc3Deseq.csv", quote=F, row.names=T)
write.csv(respc1Deseq, "respc4Deseq.csv", quote=F, row.names=T)



pdf("transcript_pc1.pdf")
superheat(
	pc1genes,
	yt=thePca$ind$coord[,1], 
	order.cols=order(thePca$ind$coord[,1]),
	yt.axis.name="PC1",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	heat.lim=c(0, 15),
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"))
dev.off()
	
pdf("transcript_pc2.pdf")
superheat(
	pc2genes,
	yt=thePca$ind$coord[,2], 
	order.cols=order(thePca$ind$coord[,2]),
	yt.axis.name="PC2",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	heat.lim=c(0, 15),
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"))
dev.off()
	
pdf("transcript_pc3.pdf")
superheat(
	pc3genes,
	yt=thePca$ind$coord[,3], 
	order.cols=order(thePca$ind$coord[,3]),
	yt.axis.name="PC3",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	heat.lim=c(0, 15),
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"))
dev.off()
	
pdf("transcript_pc4.pdf")
superheat(
	pc4genes,
	yt=thePca$ind$coord[,4], 
	order.cols=order(thePca$ind$coord[,4]),
	yt.axis.name="PC4",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,	
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	heat.lim=c(0, 15),
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"))
dev.off()
	
pc1genesScaled <- t(apply(pc1genes, 1, scale, center=T, scale=T))
pc2genesScaled <- t(apply(pc2genes, 1, scale, center=T, scale=T))
pc3genesScaled <- t(apply(pc3genes, 1, scale, center=T, scale=T))
pc4genesScaled <- t(apply(pc4genes, 1, scale, center=T, scale=T))

dimnames(pc1genesScaled)[[2]] <- dimnames(pc2genesScaled)[[2]] <- dimnames(pc3genesScaled)[[2]] <- dimnames(pc4genesScaled)[[2]] <- colnames(pc1genes)

pdf("transcript_scaled_pc1.pdf")
superheat(
	pc1genesScaled,
	yt=thePca$ind$coord[,1], 
	order.cols=order(thePca$ind$coord[,1]),
	yt.axis.name="PC1",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white",
	heat.lim=c(-2,2))
dev.off()

pdf("transcript_scaled_pc2.pdf")
superheat(
	pc2genesScaled,
	yt=thePca$ind$coord[,2], 
	order.cols=order(thePca$ind$coord[,2]),
	yt.axis.name="PC2",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white",
	heat.lim=c(-2,2))
dev.off()

pdf("transcript_scaled_pc3.pdf")
superheat(
	pc3genesScaled,
	yt=thePca$ind$coord[,3], 
	order.cols=order(thePca$ind$coord[,3]),
	yt.axis.name="PC3",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white",
	heat.lim=c(-2,2))
dev.off()

pdf("transcript_scaled_pc4.pdf")
superheat(
	pc4genesScaled,
	yt=thePca$ind$coord[,4], 
	order.cols=order(thePca$ind$coord[,4]),
	yt.axis.name="PC4",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,	
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white",
	heat.lim=c(-2,2))
dev.off()

##############################
##	PC1 without sample 6	##
##############################

pdf("transcript_scaled_pc1b.pdf")
superheat(
	pc1genesScaled[,1:5],
	yt=thePca$ind$coord[1:5,1], 
	order.cols=order(thePca$ind$coord[1:5,1]),
	yt.axis.name="PC1",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white",
	heat.lim=c(-2,2))
dev.off()

######################
##	GO enrichment	##
######################

	#GO terms
goAnnots <- list(BP=readMappings(file="./GO_annots/goBP_dm6.19.txt"), MF=readMappings(file="./GO_annots/goMF_dm6.19.txt"), CC=readMappings(file="./GO_annots/goCC_dm6.19.txt"))

	#function for GO enrichment
GOfunk3 <- function(datRef, theGenes, dagFP, resFP, annotations, cutoff, nodeSize){
	geneVec <- factor(as.integer(datRef %in% theGenes))
	names(geneVec) <- datRef

	for(GOclass in c("BP", "MF", "CC")){
	go <- new("topGOdata", ontology=GOclass, allGenes=geneVec, annot = annFUN.gene2GO, gene2GO = annotations[GOclass][[1]], nodeSize=nodeSize)
	res <- runTest(go, "classic", "fisher")
	printGraph(go, res, firstSigNodes=5, useInfo='all', fn.prefix=paste(dagFP, GOclass, sep=""), pdfSW=T)

	resS <- subset(GenTable(go, P=res, orderBy="P", topNodes=100, numChar=10000), P<=cutoff)

	write.table(resS, file=paste(resFP, GOclass, ".txt", sep=""), row.names=F, quote=F, sep="\t")
	
	resS$OE <- with(resS, Significant/Expected)
	resS$negLog <- -log10(as.numeric(resS$P))
	resS <- resS[order(resS$negLog), ]
	resS$Term <- with(resS, ifelse(nchar(Term)>ncharsToInclude, paste(substr(Term, 1, ncharsToInclude), "...", sep=""), Term))
	resS$Term <- as.factor(resS$Term)
	resS$Significant
	
goFigure <- 
	ggplot(resS, aes(x = negLog, y = reorder(Term, negLog))) + 
	geom_vline(xintercept = 1.3, col=2) +
	geom_point(aes(size=Significant, fill=OE),shape=21, alpha=0.75) +
	scale_y_discrete("GO categories") +
	scale_x_continuous(expression(paste("P-value (", -log[10],")"))) +
	scale_size_continuous("N. genes", breaks=waiver()) +
	scale_fill_gradientn(colours = pal) + 
	labs(fill="Obs. / Exp.")
  
	ggsave(paste(dagFP, GOclass, "_dotPlot.pdf", sep=""), plot = goFigure + theTheme)
	}
}

	#options for GO plot

pal <- wes_palette("Zissou1", 100, type = "continuous")
Smallfont <- 10
Margin <- 10
ncharsToInclude <- 60

theTheme <- theme(
	panel.grid.major.y = element_line(colour = grey(0.45), 
	linetype = "dotted", size = 0.2),
	panel.background = element_blank(),
	axis.title.x = element_text(size=Smallfont,colour="black"),
	axis.title.y = element_text(size=Smallfont,colour="black"),
	axis.line.x = element_line(colour="black",size=0.75),
	axis.line.y = element_line(colour="black",size=0.75),
	axis.ticks.x = element_line(size = 0.75),
	axis.ticks.y = element_line(size = 0.75),
	axis.text.x = element_text(size=Smallfont,colour="black"),
	axis.text.y = element_text(size=Smallfont-2,colour="black"),
	legend.key.height = unit(0.4, "cm"),
	legend.key.width= unit(0.6, "cm"),
	legend.margin=margin(b=0, unit='cm'),
	legend.title = element_text(face="bold",size=Smallfont), 
	legend.key = element_rect(colour = 'white', fill = "white", linetype='dashed'),
	legend.text = element_text(size=Smallfont-2),
	legend.background = element_rect(fill=NA),
	strip.text.x = element_text(size =Smallfont, colour = "black",face="italic"),
	strip.text.y = element_text(
		size = Smallfont, 
		colour = "black", 
		margin = margin(t = 3, r = 1, b = 3, l = 1)),
	strip.background = element_rect(fill=NA, colour="black"),
	strip.placement="outside")

	#options for GOfunk3
universe <- rownames(deseq_dataset)
dFP <- "./GOresults/graphs"
tFP <- "./GOresults/tables"
dir.create(dFP)
dir.create(tFP)

	#run functions
		#first analysis included nodes (nodeSize) requiring ≥3 genes
ns <- 3
CO <- 0.05
GOfunk3(datRef=universe, theGenes=pc1_DE, dagFP=paste(dFP, "pc1",sep="/"), resFP=paste(tFP, "pc1",sep="/"), annotations= goAnnots, cutoff=CO, nodeSize=ns)
GOfunk3(datRef=universe, theGenes=pc2_DE, dagFP=paste(dFP, "pc2",sep="/"), resFP=paste(tFP, "pc2",sep="/"), annotations= goAnnots, cutoff=CO, nodeSize=ns)
GOfunk3(datRef=universe, theGenes=pc3_DE, dagFP=paste(dFP, "pc3",sep="/"), resFP=paste(tFP, "pc3",sep="/"), annotations= goAnnots, cutoff=CO, nodeSize=ns)
GOfunk3(datRef=universe, theGenes=pc4_DE, dagFP=paste(dFP, "pc4",sep="/"), resFP=paste(tFP, "pc4", sep="/"), annotations= goAnnots, cutoff=CO, nodeSize=ns)

	#for the whole set of genes
GOfunk3(datRef=universe, theGenes=unique(c(pc1_DE, pc2_DE, pc3_DE, pc4_DE)), dagFP=paste(dFP, "all",sep="/"), resFP=paste(tFP, "all",sep="/"), annotations= goAnnots, cutoff=CO, nodeSize=ns)

	#replot the all-gene association
BP <- read.table("GOresults/tables/allBP.txt", header=T, sep="\t")
MF <- read.table("GOresults/tables/allMF.txt", header=T, sep="\t")
CC <- read.table("GOresults/tables/allCC.txt", header=T, sep="\t")

BP$OE <- with(BP, Significant/Expected)
BP$negLog <- -log10(BP$P)
MF$OE <- with(MF, Significant/Expected)
MF$negLog <- -log10(MF$P)
CC$OE <- with(CC, Significant/Expected)
CC$negLog <- -log10(CC$P)

truncateLabels <- function(x){
	x$Term <- with(x, ifelse(nchar(Term)>ncharsToInclude, paste(substr(Term, 1, ncharsToInclude), "...", sep=""), Term))
	x
}

replotGO <- function(x){
	ggplot(x, aes(x = negLog, y = reorder(Term, negLog))) + 
	geom_vline(xintercept = 1.3, col=2) +
	geom_point(aes(size=Significant, fill=OE),shape=21, alpha=0.75) +
	scale_y_discrete("GO categories") +
	scale_x_continuous(expression(paste("P-value (", -log[10],")")), limits=c(1,4)) +
	scale_size_continuous("N. genes", breaks=seq(0,40,5), limits=c(0,40)) +
	scale_fill_gradientn(colours = pal, limits = c(1,70)) + 
	labs(fill="Obs. / Exp.") +
	theTheme
}

BP <- truncateLabels(BP)
MF <- truncateLabels(MF)
CC <- truncateLabels(CC)

BPfig <- replotGO(BP)
MFfig <- replotGO(MF)
CCfig <- replotGO(CC)

BPfig
MFfig
CCfig

ggsave("GOresults/graphs/allBP_dotPlot.pdf", plot = BPfig)
ggsave("GOresults/graphs/allMF_dotPlot.pdf", plot = MFfig)	
ggsave("GOresults/graphs/allCC_dotPlot.pdf", plot = CCfig)

	#GO with facets
	#replot the all-gene association
BP1 <- data.frame(read.table("GOresults/tables/pc1BP.txt", header=T, sep="\t"), PC=1, ontology="BP")
BP2 <- data.frame(read.table("GOresults/tables/pc2BP.txt", header=T, sep="\t"), PC=2, ontology="BP")
BP3 <- data.frame(read.table("GOresults/tables/pc3BP.txt", header=T, sep="\t"), PC=3, ontology="BP")
BP4 <- data.frame(read.table("GOresults/tables/pc4BP.txt", header=T, sep="\t"), PC=4, ontology="BP")
MF1 <- data.frame(read.table("GOresults/tables/pc1MF.txt", header=T, sep="\t"), PC=1, ontology="MF")
MF2 <- data.frame(read.table("GOresults/tables/pc2MF.txt", header=T, sep="\t"), PC=2, ontology="MF")
MF3 <- data.frame(read.table("GOresults/tables/pc3MF.txt", header=T, sep="\t"), PC=3, ontology="MF")
MF4 <- data.frame(read.table("GOresults/tables/pc4MF.txt", header=T, sep="\t"), PC=4, ontology="MF")
CC1 <- data.frame(read.table("GOresults/tables/pc1CC.txt", header=T, sep="\t"), PC=1, ontology="CC")
CC2 <- data.frame(read.table("GOresults/tables/pc2CC.txt", header=T, sep="\t"), PC=2, ontology="CC")
CC3 <- data.frame(read.table("GOresults/tables/pc3CC.txt", header=T, sep="\t"), PC=3, ontology="CC")
#CC4 <- data.frame(read.table("GOresults/tables/pc4CC.txt", header=T, sep="\t"), PC=4, ontology="CC")

allGO <- rbind(BP1, BP2, BP3, BP4, MF1, MF2, MF3, MF4, CC1, CC2, CC3)
allGO$OE <- with(allGO, Significant/Expected)
allGO$negLog <- -log10(allGO$P)
allGO <- truncateLabels(allGO)

theTheme2 <- theme(
	panel.grid.major.y = element_line(colour = "black", 
	linetype = "solid", size = 0.025),
	panel.background = element_blank(),
	axis.title.x = element_text(size=Smallfont,colour="black"),
	axis.title.y = element_text(size=Smallfont,colour="black"),
	axis.line.x = element_line(colour="black",size=0.75),
	axis.line.y = element_line(colour="black",size=0.75),
	axis.ticks.x = element_line(size = 0.75),
	axis.ticks.y = element_line(size = 0.75),
	axis.text.x = element_text(size=Smallfont,colour="black"),
	axis.text.y = element_text(size=2,colour="black"),
	legend.position="top",
	legend.direction="horizontal",
	legend.box="horizontal",
	legend.key.height = unit(0.4, "cm"),
	legend.key.width= unit(0.6, "cm"),
	legend.margin=margin(b=0, unit='cm'),
	legend.title = element_text(face="bold",size=Smallfont), 
	legend.key = element_rect(colour = 'white', fill = "white", linetype='dashed'),
	legend.text = element_text(size=Smallfont-2),
	legend.background = element_rect(fill=NA),
	strip.text.x = element_text(size =Smallfont, colour = "black"),
	strip.text.y = element_text(
		size = Smallfont, 
		colour = "black", 
		margin = margin(t = 3, r = 1, b = 3, l = 1)),
	strip.background = element_rect(fill=NA, colour="black"),
	strip.placement="outside")

allGOplot <- ggplot(allGO, aes(x = negLog, y = reorder(Term, negLog))) + 
	geom_vline(xintercept = 1.3, col=2) +
	geom_point(aes(size=Significant, fill=OE),shape=21, alpha=0.75) +
	scale_y_discrete("GO categories") +
	scale_x_continuous(expression(paste("P-value (", -log[10],")"))) +
	scale_size_continuous("N. genes") +
	scale_fill_gradientn(colours = pal) + 
	labs(fill="Obs. / Exp.") +
	facet_grid(ontology ~ PC, scale="free_y", space="free") +
	theTheme2

ggsave("./GOresults/graphs/allGO_dotPlot_facets.pdf", 
	plot = allGOplot,
	units="mm",
	height=180,
	width=160)

	#look at co-occurence by upset
BP1.digital <- data.frame(Term=BP1$Term, PC1=1)
BP2.digital <- data.frame(Term=BP2$Term, PC2=1)
BP3.digital <- data.frame(Term=BP3$Term, PC3=1)
BP4.digital <- data.frame(Term=BP4$Term, PC4=1)

GOsets.BP <- merge(BP1.digital, BP2.digital, by="Term", all=T)
GOsets.BP <- merge(GOsets.BP, BP3.digital, by="Term", all=T)
GOsets.BP <- merge(GOsets.BP, BP4.digital, by="Term", all=T)
GOsets.BP <- as.data.frame(ifelse(is.na(GOsets.BP[,2:5]), 0, 1), row.names=GOsets.BP$Term)

MF1.digital <- data.frame(Term=MF1$Term, PC1=1)
MF2.digital <- data.frame(Term=MF2$Term, PC2=1)
MF3.digital <- data.frame(Term=MF3$Term, PC3=1)
MF4.digital <- data.frame(Term=MF4$Term, PC4=1)

GOsets.MF <- merge(MF1.digital, MF2.digital, by="Term", all=T)
GOsets.MF <- merge(GOsets.MF, MF3.digital, by="Term", all=T)
GOsets.MF <- merge(GOsets.MF, MF4.digital, by="Term", all=T)
GOsets.MF <- as.data.frame(ifelse(is.na(GOsets.MF[,2:5]), 0, 1), row.names=GOsets.MF$Term)

CC1.digital <- data.frame(Term=CC1$Term, PC1=1)
CC2.digital <- data.frame(Term=CC2$Term, PC2=1)
CC3.digital <- data.frame(Term=CC3$Term, PC3=1)

GOsets.CC <- merge(CC1.digital, CC2.digital, by="Term", all=T)
GOsets.CC <- merge(GOsets.CC, CC3.digital, by="Term", all=T)
GOsets.CC <- as.data.frame(ifelse(is.na(GOsets.CC[,2:4]), 0, 1), row.names=GOsets.CC$Term)

meUpset <- function(x, nPCs, ...){
	upset(x, 
	sets=paste("PC", 1:nPCs, sep=""), 
	keep.order=T,
	text.scale=3, ...)
}

GOsets.BP.plot <- meUpset(GOsets.BP, 4, point.size=10, line.size=5, matrix.dot.alpha=0)
GOsets.MF.plot <- meUpset(GOsets.MF, 4, point.size=10, line.size=5, matrix.dot.alpha=0)
GOsets.CC.plot <- meUpset(GOsets.CC, 3, point.size=10, line.size=5, matrix.dot.alpha=0)

GOsets.BP.plot
GOsets.MF.plot
GOsets.CC.plot

######################################################
##	does bacterial genus predict gene expression?	##
######################################################

	#create variable returning genus, or lowest known taxonomic unit
otuDecode$G2 <- with(otuDecode, 
	ifelse(S!="", 
	paste("g__", G, sep=""),
	OTU_lowest)
	)

bactOrder <- aggregate(as.matrix(OTUs) ~ G2, otuDecode, sum)
bactOrder <- droplevels(subset(bactOrder, bactOrder[,1]!=""))
rownames(bactOrder) <- bactOrder[,1]
bactOrder <- bactOrder[,2:7]
bactOrder <- t(bactOrder)
ggcorr(bactOrder)
pcaOrder <- PCA(log(bactOrder+1), scale.unit=T, graph=F)
summary(pcaOrder)
plot(pcaOrder)

pdf("PCA.percentage.explained.pdf")
par(las=2)
bp <- barplot(c(0,pcaOrder$eig[,2]), col=1, ylim=c(0,100), ylab="%age variance explained")
lines(x=bp, y=c(0,pcaOrder$eig[,3]), col="red", lwd=2)
dev.off()

pdf("PCA.pairs.genus.pdf")
pairs(pcaOrder$ind$coord, pch=16, col=brewer.pal(6, "Dark2"), cex=2, las=1, lower.panel=NULL)
dev.off()

pdf("PCA.legend.genus.pdf")
plot(1,1,col=0, axes=F)
legend("bottomleft", legend=rownames(thePca$ind$coord), col=brewer.pal(6, "Dark2"), pch=16)
dev.off()

pcaOrder_scaled <- apply(pcaOrder$ind$coord, 2, scale, center=T, scale=T)
colnames(pcaOrder_scaled) <- gsub("Dim.", "PC", colnames(pcaOrder_scaled))
pcaOrder_scaled <- data.frame(SampleName=rownames(pcaOrder$ind$coord), pcaOrder_scaled)

	# Making a deseq dataset with the data from tximport
deseq_dataset2 <- DESeqDataSetFromTximport(
	txi = count_data_PCs, 
	colData = pcaOrder_scaled, # colData is the metadata associated with the count data table
	design = ~ PC1 + PC2 + PC3 + PC4)

	#filter the dataset
notAllZero2 <- apply(assay(deseq_dataset2), 1, function(x){all(x>0)})
table(notAllZero2)

deseq_dataset2 <- deseq_dataset2[notAllZero2,]

	## PC1 ##
pc1Deseq2 <- DESeq(deseq_dataset2, test="LRT", reduced= ~ PC2 + PC3 + PC4)
respc1Deseq2 <- results(pc1Deseq2)
table(is.na(respc1Deseq2$padj))
table(is.na(p.adjust(respc1Deseq2$pvalue, method="fdr")))
respc1Deseq2$padj <- p.adjust(respc1Deseq2$pvalue, method="fdr")
sum(respc1Deseq2$padj <= pthresh, na.rm = TRUE)
pc1_DE2 <- giveGene(respc1Deseq2, pThresh)
pc1_DE2

	## PC2 ##
pc2Deseq2 <- DESeq(deseq_dataset2, test="LRT", reduced= ~ PC1 + PC3 + PC4)
respc2Deseq2 <- results(pc2Deseq2)
table(is.na(respc2Deseq2$padj))
table(is.na(p.adjust(respc2Deseq2$pvalue, method="fdr")))
respc2Deseq2$padj <- p.adjust(respc2Deseq2$pvalue, method="fdr")
sum(respc2Deseq2$padj <= pthresh, na.rm = TRUE)
pc2_DE2 <- giveGene(respc2Deseq2, pThresh)
pc2_DE2

	## PC3 ##
pc3Deseq2 <- DESeq(deseq_dataset2, test="LRT", reduced= ~ PC1 + PC2 + PC4)
respc3Deseq2 <- results(pc3Deseq2)
table(is.na(respc3Deseq2$padj))
table(is.na(p.adjust(respc3Deseq2$pvalue, method="fdr")))
respc3Deseq2$padj <- p.adjust(respc3Deseq2$pvalue, method="fdr")
sum(respc3Deseq2$padj <= pthresh, na.rm = TRUE)
pc3_DE2 <- giveGene(respc3Deseq2, pThresh)
pc3_DE2

	## PC4 ##
pc4Deseq2 <- DESeq(deseq_dataset2, test="LRT", reduced= ~ PC1 + PC2 + PC3)
respc4Deseq2 <- results(pc4Deseq2)
table(is.na(respc4Deseq2$padj))
table(is.na(p.adjust(respc4Deseq2$pvalue, method="fdr")))
respc4Deseq2$padj <- p.adjust(respc4Deseq2$pvalue, method="fdr")
sum(respc4Deseq2$padj <= pthresh, na.rm = TRUE)
pc4_DE2 <- giveGene(respc4Deseq2, pThresh)
pc4_DE2

	#mean fdr values
mean(subset(respc1Deseq2, padj <= pthresh)$padj)
mean(subset(respc2Deseq2, padj <= pthresh)$padj)
mean(subset(respc3Deseq2, padj <= pthresh)$padj)
mean(subset(respc4Deseq2, padj <= pthresh)$padj)

	#pull correlated transcripts
pc1genes2 <- assay(rlReads[rownames(rlReads) %in% pc1_DE2])
pc2genes2 <- assay(rlReads[rownames(rlReads) %in% pc2_DE2])
pc3genes2 <- assay(rlReads[rownames(rlReads) %in% pc3_DE2])
pc4genes2 <- assay(rlReads[rownames(rlReads) %in% pc4_DE2])

rownames(pc1genes2) <- genes$gene_name[match(rownames(pc1genes2), genes$gene_id)]
rownames(pc2genes2) <- genes$gene_name[match(rownames(pc2genes2), genes$gene_id)]
rownames(pc3genes2) <- genes$gene_name[match(rownames(pc3genes2), genes$gene_id)]
rownames(pc4genes2) <- genes$gene_name[match(rownames(pc4genes2), genes$gene_id)]
pc1genes2
pc2genes2
pc3genes2
pc4genes2

pc1genesScaled2 <- t(apply(pc1genes2, 1, scale, center=T, scale=T))
pc2genesScaled2 <- t(apply(pc2genes2, 1, scale, center=T, scale=T))
pc3genesScaled2 <- t(apply(pc3genes2, 1, scale, center=T, scale=T))
pc4genesScaled2 <- t(apply(pc4genes2, 1, scale, center=T, scale=T))

dimnames(pc1genesScaled2)[[2]] <- dimnames(pc2genesScaled2)[[2]] <- dimnames(pc3genesScaled2)[[2]] <- dimnames(pc4genesScaled2)[[2]] <- colnames(pc1genes2)

pdf("transcript_scaled_pc1.1.pdf")
superheat(
	pc1genesScaled2,
	yt=pcaOrder$ind$coord[,1], 
	order.cols=order(pcaOrder$ind$coord[,1]),
	yt.axis.name="PC1",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white",
	heat.lim=c(-2.5,2.5))
dev.off()

pdf("transcript_scaled_pc2.1.pdf")
superheat(
	pc2genesScaled2,
	yt=pcaOrder$ind$coord[,2], 
	order.cols=order(pcaOrder$ind$coord[,2]),
	yt.axis.name="PC2",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white",
	heat.lim=c(-2.5,2.5))
dev.off()

pdf("transcript_scaled_pc3.1.pdf")
superheat(
	pc3genesScaled2,
	yt=pcaOrder$ind$coord[,3], 
	order.cols=order(pcaOrder$ind$coord[,3]),
	yt.axis.name="PC3",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white",
	heat.lim=c(-2.5,2.5))
dev.off()

pdf("transcript_scaled_pc4.1.pdf")
superheat(
	pc4genesScaled2,
	yt=pcaOrder$ind$coord[,4], 
	order.cols=order(pcaOrder$ind$coord[,4]),
	yt.axis.name="PC4",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white",
	heat.lim=c(-2.5,2.5))
dev.off()

##############################
##	PC1 without sample 6	##
##############################

pc1genesScaled2.1 <- pc1genesScaled2[,1:5]
colnames(pc1genesScaled2.1) <- gsub("Gut", "", colnames(pc1genesScaled2.1))
pdf("transcript_scaled_pc1b_genus.pdf")
superheat(
	pc1genesScaled2[,1:5],
	yt=thePca$ind$coord[1:5,1], 
	order.cols=order(thePca$ind$coord[1:5,1]),
	yt.axis.name="PC1",
	yt.plot.type="line",
	yt.num.ticks=5,
	left.label.size=0.25,
	left.label.text.size=2.5,
	row.dendrogram = T,
	grid.hline.col="white",
	grid.vline.col="white",
	smooth.heat=T,
	grid.hline.size=0.1,
	grid.vline.size=0.1,
	bottom.label.col=brewer.pal(6, "Dark2"),
	bottom.label.text.col="white",
	heat.lim=c(-2,2))
dev.off()


######################################
##	GO for the genus-level genes		##
######################################

GOfunk3(datRef=universe, theGenes=pc1_DE2, dagFP=paste(dFP, "pc1.1",sep="/"), resFP=paste(tFP, "pc1.1",sep="/"), annotations= goAnnots, cutoff=CO, nodeSize=ns)
GOfunk3(datRef=universe, theGenes=pc2_DE2, dagFP=paste(dFP, "pc2.1",sep="/"), resFP=paste(tFP, "pc2.1",sep="/"), annotations= goAnnots, cutoff=CO, nodeSize=ns)
GOfunk3(datRef=universe, theGenes=pc3_DE2, dagFP=paste(dFP, "pc3.1",sep="/"), resFP=paste(tFP, "pc3.1",sep="/"), annotations= goAnnots, cutoff=CO, nodeSize=ns)
GOfunk3(datRef=universe, theGenes=pc4_DE2, dagFP=paste(dFP, "pc4.1",sep="/"), resFP=paste(tFP, "pc4.1", sep="/"), annotations= goAnnots, cutoff=CO, nodeSize=ns)

	#for the whole set of genes
GOfunk3(datRef=universe, theGenes=unique(c(pc1_DE2, pc2_DE2, pc3_DE2, pc4_DE2)), dagFP=paste(dFP, "all_genus",sep="/"), resFP=paste(tFP, "all_genus",sep="/"), annotations= goAnnots, cutoff=CO, nodeSize=ns)

	#replot the all-gene association
BP.1 <- read.table("GOresults/tables/all_genusBP.txt", header=T, sep="\t")
MF.1 <- read.table("GOresults/tables/all_genusMF.txt", header=T, sep="\t")
CC.1 <- read.table("GOresults/tables/all_genusCC.txt", header=T, sep="\t")

BP.1$OE <- with(BP.1, Significant/Expected)
BP.1$negLog <- -log10(BP.1$P)
MF.1$OE <- with(MF.1, Significant/Expected)
MF.1$negLog <- -log10(MF.1$P)
CC.1$OE <- with(CC.1, Significant/Expected)
CC.1$negLog <- -log10(CC.1$P)

BP.1 <- truncateLabels(BP.1)
MF.1 <- truncateLabels(MF.1)
CC.1 <- truncateLabels(CC.1)

BP.1fig <- replotGO(BP.1)
MF.1fig <- replotGO(MF.1)
CC.1fig <- replotGO(CC.1)

BPfig
MFfig
CCfig

ggsave("GOresults/graphs/allBP_dotPlot.pdf", plot = BPfig)
ggsave("GOresults/graphs/allMF_dotPlot.pdf", plot = MFfig)	
ggsave("GOresults/graphs/allCC_dotPlot.pdf", plot = CCfig)

	#GO with facets
	#replot the all-gene association
BP1.1 <- data.frame(read.table("GOresults/tables/pc1.1BP.txt", header=T, sep="\t"), PC=1, ontology="BP")
BP2.1 <- data.frame(read.table("GOresults/tables/pc2.1BP.txt", header=T, sep="\t"), PC=2, ontology="BP")
BP3.1 <- data.frame(read.table("GOresults/tables/pc3.1BP.txt", header=T, sep="\t"), PC=3, ontology="BP")
BP4.1 <- data.frame(read.table("GOresults/tables/pc4.1BP.txt", header=T, sep="\t"), PC=4, ontology="BP")
MF1.1 <- data.frame(read.table("GOresults/tables/pc1.1MF.txt", header=T, sep="\t"), PC=1, ontology="MF")
MF2.1 <- data.frame(read.table("GOresults/tables/pc2.1MF.txt", header=T, sep="\t"), PC=2, ontology="MF")
MF3.1 <- data.frame(read.table("GOresults/tables/pc3.1MF.txt", header=T, sep="\t"), PC=3, ontology="MF")
MF4.1 <- data.frame(read.table("GOresults/tables/pc4.1MF.txt", header=T, sep="\t"), PC=4, ontology="MF")
CC1.1 <- data.frame(read.table("GOresults/tables/pc1.1CC.txt", header=T, sep="\t"), PC=1, ontology="CC")
CC2.1 <- data.frame(read.table("GOresults/tables/pc2.1CC.txt", header=T, sep="\t"), PC=2, ontology="CC")
CC3.1 <- data.frame(read.table("GOresults/tables/pc3.1CC.txt", header=T, sep="\t"), PC=3, ontology="CC")
CC4.1 <- data.frame(read.table("GOresults/tables/pc4.1CC.txt", header=T, sep="\t"), PC=4, ontology="CC")

allGO.1 <- rbind(BP1.1, BP2.1, BP3.1, BP4.1, MF1.1, MF2.1, MF3.1, MF4.1, CC1.1, CC2.1, CC3.1
#, CC4.1)
)

allGO.1$OE <- with(allGO.1, Significant/Expected)
allGO.1$negLog <- -log10(allGO.1$P)
allGO.1 <- truncateLabels(allGO.1)

allGO.1plot <- ggplot(
	allGO.1, aes(x = negLog, y = reorder(Term, negLog))) + 
	geom_vline(xintercept = 1.3, col=2) +
	geom_point(aes(size=Significant, fill=OE),shape=21, alpha=0.75) +
	scale_y_discrete("GO categories") +
	scale_x_continuous(expression(paste("P-value (", -log[10],")"))) +
	scale_size_continuous("N. genes") +
	scale_fill_gradientn(colours = pal) + 
	labs(fill="Obs. / Exp.") +
	facet_grid(ontology ~ PC, scale="free_y", space="free") +
	theTheme2

ggsave("./GOresults/graphs/allGO.1_dotPlot_facets.pdf", 
	plot = allGO.1plot,
	units="mm",
	height=180,
	width=160)

	#look at the 5 most highly enriched terms per facet
allGO.1.1 <- rbind(
	BP1.1[1:5,], 
	BP2.1[1:5,], 
	BP3.1[1:5,], 
	BP4.1[1:5,], 
	MF1.1[1:5,], 
	MF2.1[1:5,], 
	MF3.1[1:5,], 
	MF4.1[1:5,], 
	CC1.1[1:5,], 
	CC2.1[1:5,], 
	CC3.1[1:5,]#, 
#	CC4.1[1:5,]
	)

	#remove NAs
nrow(allGO.1.1)
allGO.1.1 <- allGO.1.1[complete.cases(allGO.1.1),]
nrow(allGO.1.1)
	
allGO.1.1$OE <- with(allGO.1.1, Significant/Expected)
allGO.1.1$negLog <- -log10(allGO.1.1$P)
allGO.1.1 <- truncateLabels(allGO.1.1)

allGO.1.1plot <- ggplot(
	allGO.1.1, aes(x = negLog, y = reorder(Term, negLog))) + 
	geom_vline(xintercept = 1.3, col=2) +
	geom_point(aes(size=Significant, fill=OE),shape=21, alpha=0.75) +
	scale_y_discrete("GO categories") +
	scale_x_continuous(expression(paste("P-value (", -log[10],")"))) +
	scale_size_continuous("N. genes") +
	scale_fill_gradientn(colours = pal) + 
	labs(fill="Obs. / Exp.") +
	facet_grid(ontology ~ PC, scale="free_y", space="free") +
	theTheme +
	theme(legend.position="top")
allGO.1.1plot

ggsave("./GOresults/graphs/allGO.1.1_dotPlot_facets.pdf", 
	plot = allGO.1.1plot,
	units="mm",
	height=120,
	width=240)

	#look at co-occurence by upset
BP1.1.digital <- data.frame(Term=BP1.1$Term, PC1=1)
BP2.1.digital <- data.frame(Term=BP2.1$Term, PC2=1)
BP3.1.digital <- data.frame(Term=BP3.1$Term, PC3=1)
BP4.1.digital <- data.frame(Term=BP4.1$Term, PC4=1)

GOsets.BP.1 <- merge(BP1.1.digital, BP2.1.digital, by="Term", all=T)
GOsets.BP.1 <- merge(GOsets.BP.1, BP3.1.digital, by="Term", all=T)
GOsets.BP.1 <- merge(GOsets.BP.1, BP4.1.digital, by="Term", all=T)
GOsets.BP.1 <- as.data.frame(ifelse(is.na(GOsets.BP.1[,2:5]), 0, 1), row.names=GOsets.BP.1$Term)

MF1.1.digital <- data.frame(Term=MF1.1$Term, PC1=1)
MF2.1.digital <- data.frame(Term=MF2.1$Term, PC2=1)
MF3.1.digital <- data.frame(Term=MF3.1$Term, PC3=1)
MF4.1.digital <- data.frame(Term=MF4.1$Term, PC4=1)

GOsets.MF.1 <- merge(MF1.1.digital, MF2.1.digital, by="Term", all=T)
GOsets.MF.1 <- merge(GOsets.MF.1, MF3.1.digital, by="Term", all=T)
GOsets.MF.1 <- merge(GOsets.MF.1, MF4.1.digital, by="Term", all=T)
GOsets.MF.1 <- as.data.frame(ifelse(is.na(GOsets.MF.1[,2:5]), 0, 1), row.names=GOsets.MF.1$Term)

CC1.1.digital <- data.frame(Term=CC1.1$Term, PC1=1)
CC2.1.digital <- data.frame(Term=CC2.1$Term, PC2=1)
CC3.1.digital <- data.frame(Term=CC3.1$Term, PC3=1)
CC4.1.digital <- data.frame(Term=CC4.1$Term, PC4=1)

GOsets.CC.1 <- merge(CC1.1.digital, CC2.1.digital, by="Term", all=T)
GOsets.CC.1 <- merge(GOsets.CC.1, CC3.1.digital, by="Term", all=T)
GOsets.CC.1 <- merge(GOsets.CC.1, CC4.1.digital, by="Term", all=T)
GOsets.CC.1 <- as.data.frame(ifelse(is.na(GOsets.CC.1[,2:5]), 0, 1), row.names=GOsets.CC.1$Term)

GOsets.BP.1.plot <- meUpset(GOsets.BP.1, 4, point.size=10, line.size=5, matrix.dot.alpha=0)
GOsets.MF.1.plot <- meUpset(GOsets.MF.1, 4, point.size=10, line.size=5, matrix.dot.alpha=0)
GOsets.CC.1.plot <- meUpset(GOsets.CC.1, 3, point.size=10, line.size=5, matrix.dot.alpha=0)

GOsets.BP.1.plot
GOsets.MF.1.plot
GOsets.CC.1.plot

##############################
##	heatmaps of the genera	##
##############################

pcaOrderDesc <- dimdesc(pcaOrder, axes=1:5, proba=0.1)

str(pcaOrderDesc, max.level=2)
str(pcaOrderDesc$Dim.1$quanti)
pcaOrderDesc$Dim.1$quanti
pcaOrderDesc$Dim.2$quanti

	#transpose the genus matrix for heatmapping
tBactOrder <- t(bactOrder)

	#PC1 (by AD)
pdf("genus_pc1.pdf")
superheat(log(tBactOrder[rownames(tBactOrder) %in% rownames(pcaOrderDesc$Dim.1$quanti),]+1),
          yt=pcaOrder$ind$coord[,1], 
          order.cols=order(pcaOrder$ind$coord[,1]),
          yt.axis.name="PC1", 
          yt.plot.type="line",
          yt.num.ticks=5,
          row.dendrogram = T,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,12),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")
dev.off()          

	#PC2 (by AD)
pdf("genus_pc2.pdf")
superheat(log(tBactOrder[rownames(tBactOrder) %in% rownames(pcaOrderDesc$Dim.2$quanti),]+1),
          yt=pcaOrder$ind$coord[,2], 
          order.cols=order(pcaOrder$ind$coord[,2]),
          yt.axis.name="PC2", 
          yt.plot.type="line",
          yt.num.ticks=5,
          row.dendrogram = T,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,12),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")
dev.off()

	#PC3 (by AD)
pdf("genus_pc3.pdf")
superheat(log(tBactOrder[rownames(tBactOrder) %in% rownames(pcaOrderDesc$Dim.3$quanti),]+1),
          yt=pcaOrder$ind$coord[,3], 
          order.cols=order(pcaOrder$ind$coord[,3]),
          yt.axis.name="PC3", 
          yt.plot.type="line",
          yt.num.ticks=5,
          row.dendrogram = T,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,12),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")
dev.off()

	#PC4 (by AD)
pdf("genus_pc4.pdf")
superheat(log(tBactOrder[rownames(tBactOrder) %in% rownames(pcaOrderDesc$Dim.4$quanti),]+1),
          yt=pcaOrder$ind$coord[,4], 
          order.cols=order(pcaOrder$ind$coord[,4]),
          yt.axis.name="PC4", 
          yt.plot.type="line",
          yt.num.ticks=5,
          row.dendrogram = T,
          heat.col.scheme="red", 
          grid.hline.col="black", 
          grid.vline.col="black", 
          smooth.heat=T,
          heat.lim=c(1e-100,12),
          grid.hline.size=0.1,
          grid.vline.size=0.1,
          bottom.label.col=brewer.pal(6, "Dark2"),
          left.label.size=0.5,
          left.label.text.size=3,
          heat.na.col="lightgrey",
          bottom.label.text.col="white")
dev.off()

##############################
##	barplot of OTU measure	##
##############################

library("randomcoloR")

mahPal <- distinctColorPalette(k=length(unique(otuDecode$G2)))

otuG2long <- data.frame(taxon=rep(otuDecode$G2, 6), 
	sample=rep(1:6, each=nrow(otuDecode)),
	count=c(OTUs[,1],
	OTUs[,2],
	OTUs[,3],
	OTUs[,4],
	OTUs[,5],
	OTUs[,6]))

agg <- aggregate(count ~ taxon + sample, otuG2long, sum)

g1 <- ggplot(agg, aes(x = sample, y = count+1, fill = taxon))
g2 <- g1 + 
	geom_bar(stat = "identity", color="white", position="fill") +
	theme_classic() +
	theme(legend.position="none") + 
	scale_fill_viridis(option="H", discrete=T) +
	ylab("proportion") #+
#	scale_fill_manual(values=rev(mahPal)) #+ 
#	scale_y_continuous(trans="log2") #+ 
#	theme_bw()
g3 <- g2 + scale_x_discrete(breaks=1:6, limits=factor(1:6))
g3
pdf("bacterial diversity barplot.pdf", h=3, w=3)
print(g3)
dev.off() 