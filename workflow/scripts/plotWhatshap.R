#!/usr/bin/env Rscript --vanilla

args = commandArgs(trailingOnly=TRUE)
input=args[1]
output=args[2]

suppressWarnings(suppressPackageStartupMessages(library(plyr)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))

theme_set(theme_classic(base_size=10))

df <- read.table(input, sep="\t", header=FALSE)
colnames <- c("Sample", "Chromosome", "Filename", "Variants", "Phased", "Unphased", "Singletons", "Blocks", "VariantsPerBlockMedian", "VariantsPerBlockAvg", "VariantsPerBlockMin", "VariantsPerBlockMax", "VariantsPerBlockSum", "BpPerBlocKMedian", "BpPerBlockAvg", "BpPerBlocKMin", "BpPerBlockMax", "BpPerBlocKSum", "HeterozygousVariants", "HeterozygousSNVs", "PhasedSNVs", "PhasedFraction", "PhasedSNVFraction","BlockN50")
colnames(df) <- colnames

df$NotPhasedSNVs <- df$HeterozygousSNVs - df$PhasedSNVs
df$HomozygousVariants <- df$Variants - df$HeterozygousVariants
df$HeterozygousINDELs <- df$HeterozygousVariants - df$HeterozygousSNVs

dfsub <- df %>% subset(select=c(HomozygousVariants, NotPhasedSNVs, PhasedSNVs, HeterozygousINDELs, Chromosome))
dfsub$Chromosome <- factor(dfsub$Chromosome, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))
dfmelt <- dfsub %>% melt(id.vars=c("Chromosome"))

colnames(dfmelt) <- c("Chromosome", "VariantType", "Count")

blockN50 <- df %>% filter(Chromosome == "ALL") %>% pull("BlockN50")
blocks <- df %>% filter(Chromosome == "ALL") %>% pull("Blocks")

caption <- stringr::str_glue("Block N50 (kb): {signif(blockN50/1000, digits=4)}
Avg number of blocks per Chr: {blocks}")

g <- ggplot(dfmelt %>% filter(Chromosome != "ALL"), aes(x=Chromosome, fill=VariantType, y=Count)) + geom_bar(position="stack", stat="identity") +
 facet_wrap(~Chromosome, scales="free_x") + theme(strip.background = element_blank(), strip.text.x = element_blank()) +
 scale_fill_manual(values=c("dimgray", "firebrick2", "darkseagreen3", "darkslateblue"))

p <- cowplot::add_sub(g, caption, x=0, hjust=0, size=10)

cowplot::ggsave2(output, plot=p, dpi=300, height=8, width=8, unit="in")
