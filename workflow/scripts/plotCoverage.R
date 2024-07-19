#!/usr/bin/env Rscript --vanilla

args = commandArgs(trailingOnly=TRUE)
input=args[1]
output=args[2]

outputparts=stringr::str_split(output, "\\.")[[1]]
mode=outputparts[length(outputparts)]

suppressWarnings(suppressPackageStartupMessages(library(karyoploteR)))
suppressWarnings(suppressPackageStartupMessages(library(regioneR)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))

df <- read.table(input, header=FALSE, sep="\t")
colnames(df) <- c("Chromosome", "Position", "Depth", "QualString")

df$MeanQuality <- sapply(df$QualString, function(x) as.integer(mean(as.numeric(stringr::str_split(x, ",")[[1]]))))
df <- df %>% mutate(RelativeDepth=Depth/(mean(Depth)*2))

meanDepth <- mean(df$Depth)
meanRelativeDepth <- meanDepth/(mean(df$Depth)*2)
roundDepth <- round(meanDepth)

mapcolors <- brewer.pal(11, "PuOr")
mapcolors <- colorRampPalette(mapcolors)(60)
df$MapColor <- mapcolors[as.vector(df$MeanQuality)]

# change this to svg after re-doing the conda environment
if( mode == "png"){
    png(output, width=8.5, height=11, unit="in", res=300)
} else {
    svg(output, width=8.6, height=11)
}

kp <- plotKaryotype(plot.type=1)
kpSegments(kp, chr=df %>% pull(Chromosome), x0 = df %>% pull(Position), x1= df %>% pull(Position), y0=0, y1=df %>% pull(RelativeDepth), col=df %>% pull(MapColor))
kpText(kp, chr=seqlevels(kp$genome), y=meanRelativeDepth, x=0, data.panel=1, col="deeppink4", label=glue::glue("{roundDepth}X"), cex=0.4, pos=2)
kpAbline(kp, h=meanRelativeDepth, data.panel=1, col="deeppink4")
dev.off()