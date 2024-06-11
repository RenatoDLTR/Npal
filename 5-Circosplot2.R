library(circlize)
library(scales)

setwd("R:/2-NTNU-Holomuseomics/SharedAstbury/2023-12_NuclearGenome/SyntenyAndCircos/2024New")

###Data
Chrbed <- read.table("Chrs.bed", sep = "\t")
Chrbed <- Chrbed[c(1,7:14,2:6),]
#GCdata <- read.table("GCcontent100kbChronly.tsv", sep = "\t")
#GCdata <- read.table("GCcontent1kbChrOnly.tsv", sep = "\t")
GCdata <- read.table("GC10kb.tab", sep = "\t")
GCdata <- GCdata[,c(1:3,5)]
colnames(GCdata) <- c("chr","start","end","GCperc")
GCdata <- GCdata[-which(GCdata$GCperc>0.5),]

Genebed <- read.table("Gene.bed")
colnames(Genebed) <- c("chr","start","end","gene")
Repeatbed <- read.table("All_Repeat.bed")
colnames(Repeatbed) <- c("chr","start","end","repeat")


Linkdata <- read.table("PCG_Npal.PCG_Npal.anchors.simple")
Linkdata$V1 <- gsub(pattern = "_R[0-9]", replacement = "", x = Linkdata$V1)
Linkdata$V2 <- gsub(pattern = "_R[0-9]", replacement = "", x = Linkdata$V2)
Linkdata$V3 <- gsub(pattern = "_R[0-9]", replacement = "", x = Linkdata$V3)
Linkdata$V4 <- gsub(pattern = "_R[0-9]", replacement = "", x = Linkdata$V4)

View(Linkdata)
start1 <- Genebed[1,][-1,]
start2 <- Genebed[1,][-1,]
end1 <- Genebed[1,][-1,]
end2 <- Genebed[1,][-1,]


for(i in 1:nrow(Linkdata)){
  start1 <- rbind(start1, Genebed[grep(Linkdata[i,1],Genebed$gene),])
  start2 <- rbind(start2, Genebed[grep(Linkdata[i,2],Genebed$gene),])
  x <- end1
  end1 <- rbind(end1, Genebed[grep(Linkdata[i,3],Genebed$gene),])
  print(paste("For end1 number", i, tail(end1,1)[4]))
  ifelse(nrow(x)==nrow(end1),print("HEY"),print("NA"))
  end2 <- rbind(end2, Genebed[grep(Linkdata[i,4],Genebed$gene),])
}

nuc1 <- data.frame(cbind(start1$chr,start1$start,start2$end))
nuc2 <- data.frame(cbind(end1$chr,end1$start,end2$end))
colnames(nuc1) <- c("chr","start","end")
colnames(nuc2) <- c("chr","start","end")
nuc1$start <- as.numeric(nuc1$start)
nuc1$end <- as.numeric(nuc1$end)
nuc2$start <- as.numeric(nuc2$start)
nuc2$end <- as.numeric(nuc2$end)

nuc1 <- nuc1[-c(23,74),]
nuc2 <- nuc2[-c(23,74),]

rcols <- Linkdata$V6[-c(23,74)]
rcols <- ifelse(rcols=="+","darkred","darkred")
rcols <- scales::alpha(rcols, alpha=0.7)
#rcols <- scales::alpha(ifelse(sign(nuc1$start-nuc1$end) != sign(nuc2$start-nuc2$end), "#f46d43" , "#66c2a5"), alpha=0.4)

###Plot################################
set.seed(1)
circos.clear()
circos.par("track.height"=0.8, gap.degree=3, cell.padding=c(0,0,0.0,0.0), start.degree=87)
circos.initialize(sectors = Chrbed$V1, xlim = as.matrix(Chrbed[,2:3]))

set_track_gap(gap=0)

circos.track(ylim = c(0,1), panel.fun=function(x,y){
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.9, facing="bending.inside", niceFacing=TRUE, col="black")
}, bg.col="gray", bg.border = F, track.height=0.08)
brk <- c(0,6,12,18,24,30,36,42,48,54)*10^6
circos.track(track.index = get.current.track.index(), panel.fun=function(x,y){
  circos.axis(h="top", major.at=brk, labels=round(brk/10^6,1), labels.cex=0.7,
              col="black", labels.col="black",lwd=0.7,labels.facing="clockwise")
}, bg.border=F)

set_track_gap(gap=0.02)

circos.track(factors=GCdata$chr, x=GCdata$start, y=GCdata$GCperc, panel.fun=function(x,y){
  xlim=CELL_META$xlim
  circos.lines(x, y, col="black", lwd = 0.6)
  circos.segments(x0=0, x1=xlim, y0=0.25, y1=0.25, lwd=0.6, lty="11", col= "grey70")
  circos.segments(x0=0, x1=xlim, y0=0.3, y1=0.3, lwd=0.6, lty="11", col= "grey70")
  circos.segments(x0=0, x1=xlim, y0=0.35, y1=0.35, lwd=0.6, lty="11", col= "grey70")
  circos.segments(x0=0, x1=xlim, y0=0.4, y1=0.4, lwd=0.6, lty="11", col= "grey70")
  circos.segments(x0=0, x1=xlim, y0=0.45, y1=0.45, lwd=0.6, lty="11", col= "grey70")
}, ylim=c(0.2,0.45), track.height=0.1, bg.border=F)
circos.yaxis(at=c(0.3,0.4), side = "right",labels.cex=0.6, lwd=0, tick.length=0, labels.col="black", col="#FFFFFF")

set_track_gap(gap=0)

circos.genomicDensity(Repeatbed, window.size=5e5, count_by="perc", col=c("darkgreen"), track.height=0.1, ylim.force = T, bg.border=F)
circos.track(track.index = get.current.track.index(), panel.fun=function(x,y){
  xlim=CELL_META$xlim
  circos.segments(x0=0, x1=xlim, y0=0, y1=0, lwd=0.8)
  circos.segments(x0=0, x1=xlim, y0=1, y1=1, lwd=0.8)
}, bg.border=F)
circos.yaxis(at=c(0,1), side = "right",labels.cex=0.6, lwd=0.7, tick.length=0.7, labels.col="black", col="#FFFFFF")

set_track_gap(gap=0.05)

circos.genomicDensity(Genebed, window.size=5e5, count_by="perc", col = c("darkblue"), track.height=0.1, ylim.force = T, bg.border=F)
circos.track(track.index = get.current.track.index(), panel.fun=function(x,y){
  xlim=CELL_META$xlim
  circos.segments(x0=0, x1=xlim, y0=0, y1=0, lwd=0.8)
  circos.segments(x0=0, x1=xlim, y0=1, y1=1, lwd=0.8)
}, bg.border=F)
circos.yaxis(at=c(0,1), side = "right",labels.cex=0.6, lwd=0.7, tick.length=0.7, labels.col="black", col="#FFFFFF")

set_track_gap(gap = 0)

circos.genomicLink(nuc1,nuc2,col=rcols,border=NA)

