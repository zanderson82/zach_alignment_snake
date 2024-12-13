#!/usr/bin/env Rscript --vanilla

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://cran.r-project.org")

# Set Bioconductor to version 3.17
BiocManager::install(version="3.17", ask=FALSE)

# Install QDNAseq
BiocManager::install("QDNAseq", version="3.17")

if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes", repos="https://cran.r-project.org")

if (!requireNamespace("QDNAseq.hg38", quietly=TRUE)) {
    remotes::install_github("asntech/QDNAseq.hg38")
}

suppressWarnings(suppressPackageStartupMessages(library(optparse)))

# parse command line arguments
option_list <- list(make_option(c("-t", "--threads"), action="store", help="Threads to use", default=1),
                    make_option(c("-s", "--seed"), action="store", help="random seed", default=1234),
                    make_option(c("-b", "--binSize"), action="store", help="bin size: options are 1,5,10,15,30,50,100,500 or 1000kbp. Default=10", default=10)
                    )

parser <- OptionParser(usage="%prog [options] input output", option_list=option_list, description="\nRun QDNASeq with hg38 bins.\n
Example: ./qdnaseq_sop.R -b 15 -t 30 m12345.phased.bam m12345.qdnaseqoutput   -->     outputs m12345.qdnaseqoutput.called_cnv.seg, m12345.qdnaseqoutput.caled_cnv.vcf,\n m12345.qdnaseqoutput.called_cnv.pdf, and m12345.qdnaseqoutput.cnv.bins.txt with 10kb bins using 30 threads")

arguments <- parse_args(parser, positional_arguments=2)
opt <- arguments$options

# extract command line arguments
bamfile <- arguments$args[1]
outputname <- arguments$args[2]
threads <- opt$threads
binSize <- opt$binSize

# set random seed for reproducibility
set.seed(opt$seed)

# load QDNAseq package
suppressWarnings(suppressPackageStartupMessages(library(QDNAseq)))
suppressWarnings(suppressPackageStartupMessages(library(QDNAseq.hg38)))

# set parallel processing plan
future::plan("multisession", workers=threads)

# get bin annotations from cnv bin size against hg38 genome
bins <- getBinAnnotations(binSize=binSize, genome="hg38")

# count reads in bins
readCounts <- binReadCounts(bins, bamfiles=bamfile)

# removes bins that are considered low quality or problematic
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE, mappability=60)
# residual=TRUE: remove bins with high residual noise
# blacklist=TRUE: remove bins in regions known to cause artifacts
# mappability=60: filters bins with low mappability (<60%), where reads cannot be reliably aligned.

# The isobarPlot step is commented out. It would visualize some properties of read counts. 
# The noisePlot step is also commented out. It would visualize noise characteristics.
### OPTIONAL PLOTS, COMMENT OUT LATER
pdf(file=paste(outputname, "isobarplot.pdf", sep=""))
isobarPlot(readCountsFiltered)
dev.off()
### OPTIONAL PLOTS, COMMENT OUT LATER

# computes correction factors for each bin to account for systematic biases (e.g., GC content, mappability differences, etc.).
readCountsFiltered <- estimateCorrection(readCountsFiltered)

# after correction estimation, you refilter the data, this time also excluding chromosome Y. Sometimes the Y chromosome introduces bias, especially in samples that may or may not have a Y chromosome.
readCountsFiltered <- applyFilters(readCountsFiltered, chromosomes="Y", residual=TRUE, blacklist=TRUE, mappability=60)

### OPTIONAL PLOTS, COMMENT OUT LATER
pdf(file=paste(outputname, "noiseplot.pdf", sep=""))
par(mar=c(4.1, 4.4, 4.1, 1.0), xaxs="i", yaxs="i")
noisePlot(readCountsFiltered)
dev.off()
### OPTIONAL PLOTS, COMMENT OUT LATER

# corrects for biases in read counts
copyNumbers <- correctBins(readCountsFiltered)

# normalizes the corrected read counts to a common scale
copyNumbersNormalized <- normalizeBins(copyNumbers)

# smooths out outliers in the data
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