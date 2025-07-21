
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

setwd("D:/SuperEnhancer/Fig1 SE overview/SuperEnhancer")
SE_files<-list.files(pattern="*SuperEnhancers.table.txt")

peak_all<-c()
for(i in 1:length(SE_files))
{
  data<-read.table(SE_files[i],header=TRUE)
  peak<-data[,c("CHROM","START","STOP","REGION_ID")]
  peak_all<-rbind(peak_all,peak)
}
dim(peak_all)
write.csv(peak_all,file="GBM_SE_peaks.csv")
#csv-txt delete first col (number)


SE_peaks=readPeakFile('GBM_SE_peaks2.txt')
#Pie plot
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno<-annotatePeak(SE_peaks, tssRegion = c(-3000, 3000), TxDb = txdb)
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno,
              title="Distribution of histone modification-binding loci relative to TSS")

setwd("D:/SuperEnhancer/Fig1 SE overview")
pdf("Fig 1A pieplot annotation.pdf",height=4,width=7)
plotAnnoPie(peakAnno)
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
SE_peaks_tagMatrix<-getTagMatrix(SE_peaks, windows=promoter)
pdf("Fig 1B read count frequency.pdf",width=5,height=4)
plotAvgProf(SE_peaks_tagMatrix, xlim=c(-3000, 3000), 
            xlab="Genomic Region (5' -> 3')", 
            ylab = "Read Count Frequency")
dev.off()


