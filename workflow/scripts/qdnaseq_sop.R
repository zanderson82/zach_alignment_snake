#!/usr/bin/env Rscript --vanilla

suppressWarnings(suppressPackageStartupMessages(library(optparse)))

option_list <- list(make_option(c("-t", "--threads"), action="store", help="Threads to use", default=1),
                    make_option(c("-s", "--seed"), action="store", help="random seed", default=1234),
                    make_option(c("-b", "--binSize"), action="store", help="bin size: options are 1,5,10,15,30,50,100,500 or 1000kbp. Default=10", default=10)
                    )

parser <- OptionParser(usage="%prog [options] input output", option_list=option_list, description="\nRun QDNASeq with hg38 bins.\n
Example: ./qdnaseq_sop.R -b 15 -t 30 m12345.phased.bam m12345.qdnaseqoutput   -->     outputs m12345.qdnaseqoutput.called_cnv.seg, m12345.qdnaseqoutput.caled_cnv.vcf,\n m12345.qdnaseqoutput.called_cnv.pdf, and m12345.qdnaseqoutput.cnv.bins.txt with 10kb bins using 30 threads")

arguments <- parse_args(parser, positional_arguments=2)
opt <- arguments$options

bamfile <- arguments$args[1]
outputname <- arguments$args[2]
threads <- opt$threads
binSize <- opt$binSize

set.seed(opt$seed)

suppressWarnings(suppressPackageStartupMessages(library(QDNAseq)))
suppressWarnings(suppressPackageStartupMessages(library(QDNAseq.hg38)))

future::plan("multisession", workers=threads)

bins <- getBinAnnotations(binSize=binSize, genome="hg38")
readCounts <- binReadCounts(bins, bamfiles=bamfile)
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE, mappability=60)

#pdf(file=paste(outputname, "isobarplot.pdf", sep=""))
#par(mar=c(4.1, 4.4, 4.1, 1.0), xaxs="i", yaxs="i")
#isobarPlot(readCountsFiltered)
#dev.off()

readCountsFiltered <- estimateCorrection(readCountsFiltered)
readCountsFiltered <- applyFilters(readCountsFiltered, chromosomes="Y", residual=TRUE, blacklist=TRUE, mappability=60)

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

copyNumbersCalled <- callBins(copyNumbersSegmented, method="cutoff")

pdf(file=paste(outputname, ".called_cnv.pdf", sep=""))
par(mar=c(4.1, 4.4, 4.1, 1.0), xaxs="i", yaxs="i")
plot(copyNumbersCalled)
dev.off()

exportBins(copyNumbersCalled, format="seg", file=paste(outputname, ".called_cnv.seg", sep=""))
exportBins(copyNumbersCalled, format="vcf", file=paste(outputname, ".called_cnv.vcf", sep=""))


