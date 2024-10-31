#!/usr/bin/env Rscript --vanilla

args = commandArgs(trailingOnly=TRUE)
input=args[1]
output=args[2]
donestring=args[3]

outputparts=stringr::str_split(output, "\\.")[[1]]
mode=outputparts[length(outputparts)]

suppressWarnings(suppressPackageStartupMessages(library(karyoploteR)))
suppressWarnings(suppressPackageStartupMessages(library(regioneR)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))

df <- read.csv(input, header=TRUE)

segdf <- read.table("/n/dat/hg38/segDups.hg38.bed", header=FALSE, sep="\t")
colnames(segdf) <- c("Chromosome", "Start", "Stop", "Label", "Unknown", "Strand")
segrange <- toGRanges(segdf, format="BED")

repeatdf <- read.table("/n/dat/hg38/repeats.hg38.sorted.bed", sep="\t", header=FALSE)
colnames(repeatdf) <- c("Chromosome", "Start", "Stop", "Type")
repeatdf <- repeatdf %>% drop_na()
reprange <- toGRanges(repeatdf, format="BED")

# change this to svg after re-doing the conda environment

genes <- unique(df$Region)
for(gene in genes){
    minioutput=paste(output,".",gene,".",mode, sep="")
    if( mode == "png"){
        png(minioutput, width=4.5, height=5, unit="in", res=300)
    } else {
        svg(minioutput, width=4, height=4)
    }
    nomodf <- df %>% filter(Region==gene) ## replace with iteration and grid arrange
    if(nrow(nomodf) > 0){
        nomodf$RelativeDepth <- nomodf$Depth / max(nomodf$Depth)
        nomodf$RelativeHP1 <- nomodf$HP1 / max(nomodf$Depth)
        nomodf$RelativeHP2 <- nomodf$HP2 / max(nomodf$Depth)

        meanRelativeDepth <- mean(nomodf$RelativeDepth)
        plotstart <- min(nomodf$Position) + 1000
        plotend <- max(nomodf$Position) + 1000
        chrom <- nomodf$Chromosome[1]

        gdf <- data.frame(chrom, plotstart, plotend)
        gdf <- gdf %>% drop_na()
        if(nrow(gdf) > 4){
            regionbounds=toGRanges(data.frame(chrom, plotstart, plotend))
            print(paste("plotting karyotype for", chrom, plotstart, plotend))
            kp <- plotKaryotype(plot.type=2, chromosomes=c(chrom), zoom=regionbounds, main=gene, cex=0.6)

            kpPolygon(kp, chr=chrom, x=c(plotstart, nomodf$Position, plotend), y=c(0, nomodf$RelativeDepth, 0), col="gray", border=NA)
            kpPolygon(kp, chr=chrom, x=c(plotstart, nomodf$Position, plotend), y=c(0, nomodf$RelativeHP1+nomodf$RelativeHP2, 0), col="gold", border=NA)
            kpPolygon(kp, chr=chrom, x=c(plotstart, nomodf$Position, plotend), y=c(0, nomodf$RelativeHP1, 0), col="cornflowerblue", border=NA)

            kpAxis(kp, numticks=5, r0=0, r1=1, ymin=0, ymax=max(nomodf$Depth)+1, cex=0.5)
            kpAbline(kp, h=mean(nomodf$RelativeDepth), col="black", lty=2)

            kpRect(kp, data=segrange, y0=0, y1=1, col="lightskyblue1", border=NA, data.panel=2, r0=0.5, r1=0.75)
            kpPlotRegions(kp, data=repeatdf, border=NA, data.panel=2, col="#666666", r0=0.2, r1=0.5)

            kpAddLabels(kp, labels="Repeats", data.panel = 2, cex=0.5,col="#666666", r0=0.2, r1=0.5)
            kpAddLabels(kp, labels="Seg Dups", data.panel = 2, cex=0.5,col="dodgerblue3", r0=0.5, r1=0.75)

            plotrange=plotend-plotstart
            if(plotrange > 150000){
                tickdist=100000
                minortickdist=10000
                digits=2
            } else if(plotrange < 30000){
                tickdist=5000
                minortickdist=1000
                digits=3
            } else {
                tickdist=10000
                minortickdist=1000
                digits=2
            }

            kpAddBaseNumbers(kp, tick.dist = tickdist, cex=0.5, minor.tick.dist = minortickdist, minor.tick.col = "black", clipping=FALSE, digits=digits)
            dev.off()
        }
    else {
        print(paste("Could not make plot for", gene))
        dev.off()
    } }
}

write("complete", donestring)
