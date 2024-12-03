rm(list=ls())

library(networkD3)
library(tidyverse)

getwd()
setwd("~/Dobslab Dropbox/Adam Dobson/students/LlanwarneFrances/frances_manuscript")

sank <- read.table("OTUs_for_sankey.txt", sep="\t", header=T, stringsAsFactors=T)
dim(sank)

sankInd <- sank[,1]
sank <- sank[,c(2:ncol(sank))]

	#calculate %age variation attributed to each level of OTU hierarchy
kPerc <- aggregate(percReads ~ k, sank, sum)
pPerc <- aggregate(percReads ~ k/p, sank, sum)
cPerc <- aggregate(percReads ~ p/c, sank, sum)
oPerc <- aggregate(percReads ~ c/o, sank, sum)
fPerc <- aggregate(percReads ~ o/f, sank, sum)
gPerc <- aggregate(percReads ~ f/g, sank, sum)
sPerc <- aggregate(percReads ~ g/s, sank, sum)

myPerc <- function(x){
	y <- x[,2]
	z <- x[,1]
	rnded <- round(x$percReads, 1)
	rnded <- ifelse(rnded<0.1, "<0.1", rnded)
	x[,2] <- with(x, {
		paste(y, " (", rnded, ")", sep="")
	})
	rnded2 <- round(tapply(x$percReads, z, sum), 1)
	x[,1] <- with(x, {
		paste(z, " (", rnded2, ")", sep="")
	})
	x
}

	#remove null nodes
oPerc <- droplevels(subset(oPerc, c!="" & o!=""))
fPerc <- droplevels(subset(fPerc, o!="" & f!=""))
gPerc <- droplevels(subset(gPerc, f!="" & g!=""))
sPerc <- droplevels(subset(sPerc, g!="" & s!=""))

#kPerc$k <- "Bacteria (100)"
#pPerc <- myPerc(pPerc)
#cPerc <- myPerc(cPerc)
#oPerc <- myPerc(oPerc)
#fPerc <- myPerc(fPerc)
#gPerc <- myPerc(gPerc)
#sPerc <- myPerc(sPerc)

pPerc$k <- factor(paste("K:", pPerc$k, sep=""))
pPerc$p <- factor(paste("P:", pPerc$p, sep=""))

cPerc$p <- factor(paste("P:", cPerc$p, sep=""))
cPerc$c <- factor(paste("C:", cPerc$c, sep=""))

oPerc$c <- factor(paste("C:", oPerc$c, sep=""))
oPerc$o <- factor(paste("O:", oPerc$o, sep=""))

fPerc$o <- factor(paste("O:", fPerc$o, sep=""))
fPerc$f <- factor(paste("F:", fPerc$f, sep=""))

gPerc$f <- factor(paste("F:", gPerc$f, sep=""))
gPerc$g <- factor(paste("G:", gPerc$g, sep=""))

sPerc$g <- factor(paste("G:", sPerc$g, sep=""))
sPerc$s <- factor(paste("S:", sPerc$s, sep=""))


colnames(pPerc) <- colnames(cPerc) <- colnames(oPerc) <- colnames(fPerc) <- colnames(gPerc) <- colnames(sPerc) <- c("source", "target", "value")



percs <- rbind(pPerc, cPerc, oPerc, fPerc, gPerc, sPerc)

sankBank <- sank[,1:7]
unique(with(sankBank, paste(p,c)))

	#	A connection data frame is a list of flows with intensity for each flow
meLinks <- data.frame(source=percs$source,
			target=percs$target,
			value=percs$value)
dim(meLinks)			
meLinks <- subset(meLinks, as.character(source)!="")
nrow(meLinks)
meLinks <- subset(meLinks, as.character(target)!="")
nrow(meLinks)



# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(meLinks$source), 
  as.character(meLinks$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the meLinks dataframe.. So we need to reformat it.
meLinks$IDsource <- match(meLinks$source, nodes$name)-1 
meLinks$IDtarget <- match(meLinks$target, nodes$name)-1
 
# Make the Network
p <- sankeyNetwork(Links = meLinks, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE)
p