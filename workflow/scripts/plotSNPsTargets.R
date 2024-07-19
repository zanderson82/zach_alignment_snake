#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
input=args[1]
output=args[2]
multisnps=args[3]
highlights=args[4]
regionlist=args[5]

suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

df <- read.table(input, header=TRUE, sep="\t")
#segdf <- read.table("/n/dat/hg38/segDups.hg38.bed", header=FALSE, sep="\t")
#colnames(segdf) <- c("Chromosome", "Start", "Stop", "Label", "Unknown", "Strand")
#centrodf <- read.table("/n/dat/hg38/centromeres.txt", header=FALSE, sep="\t")
#colnames(centrodf) <- c("Bin", "Chromosome", "Start", "Stop", "Label")
exclusions <- read.table("/n/dat/hg38/all_exclusions.sorted.merged.bed", header=FALSE, sep="\t")
colnames(exclusions) <- c("Chromosome", "Start", "Stop")

highlightdf <- read.table(highlights, header=TRUE, sep="\t")

regions=stringr::str_split(regionlist, ";")[[1]]

theme_set(theme_classic(base_size=12))

chromosomes <- paste("chr",c(seq(1,22),"X","Y"), sep="")

df[which(df$Genotype=="1/1"),"GTQuality"] <- df[which(df$Genotype=="1/1"),"GTQuality"]*-1

allplots=list()
for(region in regions){
    #identify runs
    chr=stringr::str_split(region,":")[[1]][1]
    coords=stringr::str_split(region,":")[[1]][2]
    startcoord=as.integer(stringr::str_split(coords,"-")[[1]][1])
    stopcoord=as.integer(stringr::str_split(coords,"-")[[1]][2])
    dfsub <- df %>% filter(Chromosome==chr) %>% filter(Position >= startcoord) %>% filter(Position <= stopcoord)
    if(nrow(dfsub)>10000){
        dfsub <- dfsub %>% slice_sample(n=10000, replace=FALSE)
    }
    if(multisnps==0){
        g <- ggplot() + geom_rect(data=highlightdf %>% filter(Chromosome==chr) %>% filter(Stop >= startcoord) %>% filter(Start <= stopcoord), aes(xmin=Start, xmax=Stop, ymin=-60, ymax=60), fill="yellow") +
        geom_rect(data=exclusions %>% filter(Chromosome==chr) %>% filter(Stop >= startcoord) %>% filter(Start <= stopcoord), aes(xmin=Start, xmax=Stop, ymin=-60, ymax=60), fill="gray") + 
        geom_point(data=dfsub %>% filter(Genotype %in% c("1/1","0/1")), aes(x=Position, y=GTQuality), shape=19) + 
        scale_x_continuous(labels=scales::label_bytes("kB"), name="Coordinates") + ylab("Homozygous     |     Heterozygous") + geom_hline(yintercept=0) +
        coord_cartesian(ylim=c(-60,60)) + ggtitle(region)
    } else {
        g <- ggplot() + geom_rect(data=highlightdf %>% filter(Chromosome==chr) %>% filter(Stop >= startcoord) %>% filter(Start <= stopcoord), aes(xmin=Start, xmax=Stop, ymin=-60, ymax=60), fill="yellow") +
        geom_rect(data=exclusions %>% filter(Chromosome==chr) %>% filter(Stop >= startcoord) %>% filter(Start <= stopcoord), aes(xmin=Start, xmax=Stop, ymin=-60, ymax=60), fill="gray") + 
        geom_point(data=dfsub %>% filter(Genotype %in% c("1/1", "0/1", "1/2")), aes(x=Position, y=GTQuality, color=Genotype), shape=19) + 
        scale_x_continuous(labels=scales::label_bytes("kB"), name="coordinates") + ylab("Homozygous    |    Heterozygous") +
        geom_hline(yintercept=0) + coord_cartesian(ylim=c(-60,60)) + theme(legend.position="None") +
        scale_color_manual(values=c("black", "black", "blue")) + ggtitle(region)
    }
    allplots=c(allplots, list(g))
}

plotlist=list(plots=allplots, num=length(regions))
do.call(cowplot::plot_grid, c(plotlist$plots))

ggsave(output, dpi=300, width=8.5, height=10, units="in")


