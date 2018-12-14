#Coursework 1

setwd("C:/Desktop/Bioinformatics/BioInf2018CourseWork-master")

#Exercise 1
#Question 1: 4719 genes have a padj value <0.05
DEtable<-read.delim("DE_Genes/GM12878_minus_HeLa_DEG.csv",sep=",",header=TRUE)
#Create significant genes table
significant<-DEtable[DEtable$padj < 0.05,]
numSignificant<-length(row.names(significant))

#Question 2:
library(ggplot2)
#Add -log10(Exp) to the significant genes table
significant$MinusLog10Pvalue<-log(significant$pvalue,base=10)
#Plot log2FoldChange v -log10(Exp) for significant genes
signifDEplot<-ggplot(significant,aes(x=log2FoldChange,y=MinusLog10Pvalue))+geom_point()+ylab("-log10 of Pvalue")+ggtitle("Volcano plot of GM12878 Minus HeLa")+theme_bw()+theme(panel.border = element_blank())
signifDEplot

#Question 3:
#Read in absExpTable
absExpTable<-read.delim("DE_Genes/Expression.csv",sep=",",header=TRUE)
#Add 1 to all expression values
absExpTable[,c(2:5)]<-absExpTable[,c(2:5)] + 1
#Add a Sample column for each tissue as a separate absExp table
absExp1<-data.frame(ID=absExpTable$ID,Sample="GM12878_1",Expression=absExpTable$GM12878_1)
absExp2<-data.frame(ID=absExpTable$ID,Sample="GM12878_2",Expression=absExpTable$GM12878_2)
absExp3<-data.frame(ID=absExpTable$ID,Sample="HeLa_1",Expression=absExpTable$HeLa_1)
absExp4<-data.frame(ID=absExpTable$ID,Sample="HeLa_2",Expression=absExpTable$HeLa_2)
#Combine absExp tables by ID
AbsExp12<-rbind(absExp1,absExp2)
AbsExp34<-rbind(absExp3,absExp4)
combAbsExpDF<-rbind(AbsExp12,AbsExp34)
#Add -log10(Exp) to combAbsExp table
combAbsExpDF$Expression<-log10(combAbsExpDF$Expression)
#Base graphics boxplot Expression v Sample
boxplot(Expression~Sample,data=combAbsExpDF,names=c("GM12878_1","GM12878_2","HeLa_1","HeLa_2"),range=0)

#Question 4:
#Subset significant genes (padj<0.05) by those with log2FoldChange>1
log2SignifDF<-significant[significant$log2FoldChange>1,]
#Subset log2significant genes table by those with padj<0.05
signifAbsExpTable<-merge.data.frame(combAbsExpDF,log2SignifDF,by="ID")
#Base graphics boxplot padj Exp v Sample < 0.05 and a log2FoldChange > 1
boxplot(Expression~Sample,data=signifAbsExpTable,names=c("GM12878_1","GM12878_2","HeLa_1","HeLa_2"))

#Question 5:
#Subset absExp table by those with expression in top 60%
expVect<-as.vector(combAbsExpDF$Expression)
maxExp<-max(expVect)
sixtyPercExp<-(maxExp*0.6)
sixtyPercExpDF<-combAbsExpDF[combAbsExpDF$Expression>=sixtyPercExp,]
#Subset DE genes table by those with at least 60 percent absExp
sixtyPercDEdf<-merge.data.frame(DEtable,sixtyPercExpDF,by="ID")
#Add log2(BaseMean)
sixtyPercDEdf$log2BaseMean<-log2(sixtyPercDEdf$baseMean)
#Add significance logical column
sigLogical<-sixtyPercDEdf$padj<0.05
sixtyPercDEdf$sigOrNot<-sigLogical
#Plot DE log2BaseMean v log2FoldChange of top 60% DE genes, color by significance
sixtyPercDEplot<-ggplot(sixtyPercDEdf,aes(x=log2BaseMean,y=log2FoldChange,color=sigOrNot))+geom_point()
sixtyPercDEplot
#Reproducing graph only possible by using all genes
DEtable$log2BaseMean<-log2(DEtable$baseMean)
#Add significance logical column
sigLogical<-DEtable$padj<0.05
DEtable$sigOrNot<-sigLogical
AllDEplot<-ggplot(DEtable,aes(x=log2BaseMean,y=log2FoldChange,color=sigOrNot))+geom_point()+theme_bw()+ggtitle("MA Plot of GM12878 Minus HeLa")+theme(panel.border=element_blank())
AllDEplot

#Exercise 2
#Question 1 Read in H3K27Ac_Limb_1.txt and report number of genomic locations there
#Read in H3K27Ac_Limb_1.txt
H3k27Table<-read.table("HOMER_peaks/H3K27Ac_Limb_1.txt",sep="\t",col.names = c("PeakID","chr","start","end","strand","Normalized Tag Count","region size","findPeaks Score","Total Tags","Control Tags (normalized to IP Experiment","Fold Change vs Control","p-value vs Control","Clonal Fold Change"))
numberLocations<-length(H3k27Table$PeakID)

#Question 2 Basegraphics histogram of log10(region size)
H3k27Table$Log10RegionSize<-log10(H3k27Table$region.size)
regionSizeHist<-hist(H3k27Table$Log10RegionSize,xlab="region size (log10)",main="Histogram log10 of region size")

#Question 3 ggplot density plot log10 regions sizes
H3k27DensityPlot<-ggplot(H3k27Table,aes(x=H3k27Table$Log10RegionSize))+geom_density(color="dark green",fill="dark green")+theme_bw()+xlab("region size (log10)")
H3k27DensityPlot

#Question 4 ggplot density plot by chr
H3k27ChrDensityPlot<-H3k27DensityPlot+facet_wrap(H3k27Table$chr)+theme(strip.background = element_rect(color="white",fill="white"),panel.border = element_blank())
H3k27ChrDensityPlot

#Question 5 boxplot of log10 findPeaks.Score for each chr
H3k27Table$log10findPeaks.Score<-log10(H3k27Table$findPeaks.Score)
H3k27ChrfindPeaks<-ggplot(H3k27Table,aes(x=H3k27Table$chr,y=H3k27Table$log10findPeaks.Score,fill=H3k27Table$chr))+geom_boxplot()+coord_flip()+theme_bw()+theme(panel.border = element_blank())
H3k27ChrfindPeaks

#Question 6 export HOMER genomic regions as a BED3 file
library(GenomicRanges)
library(rtracklayer)
IRange<-IRanges(start=H3k27Table$start,end=H3k27Table$end)
HomerGRange<-GRanges(seqnames=H3k27Table$chr,IRange)
export.bed(HomerGRange,con="HOMERGRange.bed")