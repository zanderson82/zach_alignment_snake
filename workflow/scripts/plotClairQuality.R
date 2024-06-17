#!/usr/bin/env Rscript --vanilla

args = commandArgs(trailingOnly=TRUE)
input=args[1]
output=args[2]
type=args[3]
plotcolor=args[4]

n=1000000

suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

theme_set(theme_classic(base_size=12))

df <- read.table(input, header=TRUE, sep="\t")

if(nrow(df) > n){
    df <- df %>% slice_sample(n=n, replace=FALSE)
}


typemean <- mean(df %>% filter(Type==type) %>% pull(Quality), na.rm=TRUE)
typemedian <- median(df %>% filter(Type==type) %>% pull(Quality), na.rm=TRUE)

g <- ggplot(df %>% filter(Type==type), aes(x=Quality)) + geom_density(alpha=0.5, fill=plotcolor, color=plotcolor, kernel="cosine", adjust=3) +
ggtitle(type) + xlab("Clair3 Quality Scores") + ylab("Frequency") + geom_vline(xintercept=typemean, linetype="dashed", color=plotcolor) +
geom_vline(xintercept=typemedian, linetype="solid", color=plotcolor)

legend <- stringr::str_glue("Mean: {signif(typemean, digits=3)} (dashed line)
    Median: {signif(typemedian, digits=3)} (solid line)")

p <- cowplot::add_sub(g, legend, x=0, hjust=0, size=10)

print(paste("saving plot to ", output, sep=""))

cowplot::ggsave2(output, plot=p, dpi=300, height=8, width=8, unit="in")