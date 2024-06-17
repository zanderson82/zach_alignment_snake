#!/usr/bin/env Rscript --vanilla

args = commandArgs(trailingOnly=TRUE)
input=args[1]
output=args[2]
plotcolor=args[3]

bins=150

suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

theme_set(theme_classic(base_size=12))

df <- read.table(input, header=TRUE, sep="\t")

halfgb <- sum(df$Length * df$Count)/2
sumdf <- df %>% arrange(desc(Length)) %>% dplyr::mutate("CumulativeKB" = cumsum(as.numeric(Length) * as.numeric(Count))) %>% summarise(minL=nth(Length, which.min(abs(CumulativeKB-halfgb))))
n50 <- as.numeric(sumdf[1])

bins <- seq(from=min(df$Length), to=max(df$Length), length.out=bins)
df$bincode <- .bincode(df$Length, bins, right=FALSE)

df.cumu <- df %>% dplyr::group_by(bincode) %>% dplyr::mutate("CumulativeKB" = cumsum(as.numeric(Length) * as.numeric(Count))) %>% dplyr::group_by(bincode) %>% dplyr::slice_max(CumulativeKB)
df.cumu <- data.frame(df.cumu)

names(bins) <- seq(1,length(bins))
df.cumu$bin.label <- bins[as.vector(df.cumu$bincode)]

print("making plot")

g <- ggplot(df.cumu, aes(x=bin.label, y=CumulativeKB))+ geom_bar(stat="identity", fill=plotcolor) + 
ggtitle("Read Lengths") + xlab("Read Length") + ylab("Bases") + geom_vline(xintercept=as.numeric(n50), linetype="dotted", linewidth=1.2)

g <- g + scale_x_continuous(label=scales::label_bytes(), breaks=waiver(), n.breaks=20) + coord_cartesian(xlim=c(1,165000)) + theme(axis.text.x=element_text(angle=45, hjust=1)) +
scale_y_continuous(label=scales::label_bytes(), breaks=waiver())


print(paste("saving plot to ", output, sep=""))
ggsave(output, plot=g, dpi=300, height=8, width=8, unit="in")

print("Done!")