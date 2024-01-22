
#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
bamfile=args[1]
outputname=args[2]
threads=args[3]

library(QDNAseq)
library(QDNAseq.hg38)

future::plan("multisession", workers=threads)

bins <- getBinAnnotations(binSize=10, genome="hg38")
readCounts <- binReadCounts(bins, bamfiles=bamfile)
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

#pdf(file=paste(outputname, "isobarplot.pdf", sep=""))
#par(mar=c(4.1, 4.4, 4.1, 1.0), xaxs="i", yaxs="i")
#isobarPlot(readCountsFiltered)
#dev.off()

readCountsFiltered <- estimateCorrection(readCountsFiltered)

#pdf(file=paste(outputname, "noiseplot.pdf", sep=""))
#par(mar=c(4.1, 4.4, 4.1, 1.0), xaxs="i", yaxs="i")
#noisePlot(readCountsFiltered)
#dev.off()

copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

#pdf(file=paste(outputname, "copy_numbers.pdf", sep=""))
#par(mar=c(4.1, 4.4, 4.1, 1.0), xaxs="i", yaxs="i")
#plot(copyNumbersSmooth)
#dev.off()

exportBins(copyNumbersSmooth, file=paste(outputname, ".cnv.bins.txt", sep=""))

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

#pdf(file=paste(outputname, "copy_numbers_segmented.pdf", sep=""))
#par(mar=c(4.1, 4.4, 4.1, 1.0), xaxs="i", yaxs="i")
#plot(copyNumbersSegmented)
#dev.off()

copyNumbersCalled <- callBins(copyNumbersSegmented)

pdf(file=paste(outputname, ".called_cnv.pdf", sep=""))
par(mar=c(4.1, 4.4, 4.1, 1.0), xaxs="i", yaxs="i")
plot(copyNumbersCalled)
dev.off()

exportBins(copyNumbersCalled, format="seg", file=paste(outputname, ".called_cnv.seg", sep=""))
exportBins(copyNumbersCalled, format="vcf", file=paste(outputname, ".called_cnv.vcf", sep=""))


