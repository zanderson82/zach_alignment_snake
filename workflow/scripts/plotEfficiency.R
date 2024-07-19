#!/usr/bin/env Rscript --vanilla

args = commandArgs(trailingOnly=TRUE)
input=args[1]
output=args[2]
type=args[3]

suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

df <- read.table(input, header=TRUE, sep="\t")

df$Stage <- factor(df$Stage, levels=c("Sequencing", "Basecalling", "Alignment"))
theme_set(theme_classic(base_size=12, base_family="sans-serif"))

g <- ggplot(df %>% filter(Stat==type), aes(x=Stage, y=Value, fill=Stage)) + geom_bar(stat="identity") + ylab(type) + 
theme(legend.position="none") + scale_fill_manual(values=c("darkseagreen3", "darkslategray3", "dodgerblue3"))

if( type != "Reads"){
    g <- g + scale_y_continuous(labels=scales::label_bytes())
}

outputparts=stringr::str_split(output, "\\.")[[1]]
mode=outputparts[length(outputparts)]

if(mode == "png"){
    ggsave(output, plot=g, width=6, height=6, unit="in", dpi=300)
} else {
    ggsave(output, plot=g, width=6, height=6)
}
